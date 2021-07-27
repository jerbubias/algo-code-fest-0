#include "qap.h"
#include <float.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct state_size {
    int present;
};

struct sample_size {
    int add;
    int remove;
};

struct problem {
    int n;        /* Number of facilities and locations */
    int e;        /* Number of components */
    double *flow; /* Weight matrix */
    double *dist; /* Distance matrix */
};

struct solution {
    // Problem Abstraction
    struct problem *prob;
    int n_g;                      /* Number of components */
    int *groundSet;               /* Set of components */
    int *igroundSet;              /* Inverse permutation of the set of components */
    int n_p;                      /* Number of present components */
    int n_f;                      /* Number of forbidden components */
    int evalv;                    /* Flag indicating if the solution is evaluated */
    double objv;                  /* Objective value */
    int evalLB;                   /* Flag indicating if the lower bound is calculated */
    double objLB;                 /* Objective Lower Bound */
    struct state_size numEnumLim; /* Number of components left to enumerate */
    // Problem-Specific
    int n;                            /* Number of facilities and locations */
    int *facil;                       /* Permutation mapping the facilities */
    int *ifacil;                      /* Inverse permutation mapping the facilities */
    int *local;                       /* Permutation mapping the locations */
    int *ilocal;                      /* Inverse permutation mapping the locations */
    int n_a;                          /* Number of assignments */
    int *rndSample;                   /* Array used for sampling without replacement for moves of the subneighbourhood Add */
    struct sample_size sampleWORLim;  /* Number of unused moves for random move without replacement */
    struct sample_size sampleEnumLim; /* Number of unused moves for enumerate move */
    struct sample_size sampleSize;    /* Size of neighbourhood */
};

struct move {
    struct problem *prob;
    int data;       /* Component identifier */
    int facil;      /* Facility to be assigned to a location */
    int local;      /* Location to assign for next facility */
    int evalLBi[4]; /* Flag indicating if lower bound increment is evaluated for operations: { 0 - Add, 1 - Remove, 2 - Forbid, 3 - Permit } */
    double objLBi;  /* Lower bound increment */
};

extern gsl_rng *rng; /* The single rng instance used by the whole code */

/*********************************/
/* ----- Utility functions ----- */
/*********************************/

/*
 * Return random integer in the range 0..N
 */
static int randint(const int n_max)
{
    return gsl_rng_uniform_int(rng, n_max + 1);
}

/*
 * Uniquely encode an ordered pair (x, y) into a single natural number z 
 */
static int pairing(const int x, const int y, const int n)
{
    return x * n + y;
}

/*
 * Uniquely decode a single natural number z into an ordered pair (x, y)
 */
static void unpairing(int *x, int *y, const int z, const int n)
{
    *x = z / n;
    *y = z % n;
}

/*
 * Compute the sum of squares of the n consecutive natural numbers
 */
static int sumsquares(int n)
{
    return n * (2 * n + 1) * (n + 1) / 6;
}

/*
 * Exchange the values of the ith and jth elements of an array
 */
static void swap(int *data, const int i, const int j)
{
    if (i == j)
        return;
    int val = data[i];
    data[i] = data[j];
    data[j] = val;
}

/*
 * Exchange the values of the ith and jth elements of a permutation.
 * Update the inverse permutation.
 */
static void swap_i(int *data, int *idata, const int i, const int j)
{
    if (i == j)
        return;
    swap(data, i, j);
    swap(idata, *(data + i), *(data + j));
}

/*
 * Compute the number of components of a solution that are in each state
 */
static struct state_size st_size(const struct solution *s)
{
    struct state_size size;
    size.present = s->n_p;
    return size;
}

/*
 * Compute the number of neighbours of a solution that are in each
 * subneighbourhood
 */
static struct sample_size nh_size(const struct solution *s)
{
    struct sample_size size;
    size.add = (s->n - s->n_a) * (s->n - s->n_a);
    size.remove = s->n_p;
    return size;
}

/*
 * Compute lower bound
 */
static double lb(const struct solution *s)
{
    int i, j;
    double obj = 0.0;
    for (i = 0; i < s->n_a; ++i)
        for (j = 0; j < s->n_a; ++j)
            obj += s->prob->flow[s->facil[i] * s->prob->n + s->facil[j]] * s->prob->dist[s->local[i] * s->prob->n + s->local[j]];
    return obj;
}

/*
 * Compute lower bound increment after adding an assignment to a solution
 */
static double lbi_add(const struct move *v, const struct solution *s)
{
    int i;
    double obj = 0.0;
    obj += s->prob->flow[v->facil * s->prob->n + v->facil] * s->prob->dist[v->local * s->prob->n + v->local];
    for (i = 0; i < s->n_a; ++i)
        obj += s->prob->flow[v->facil * s->prob->n + s->facil[i]] * s->prob->dist[v->local * s->prob->n + s->local[i]] + s->prob->flow[s->facil[i] * s->prob->n + v->facil] * s->prob->dist[s->local[i] * s->prob->n + v->local];
    return obj;
}

/*
 * Compute lower bound increment after removing an assignment from a solution
 */
static double lbi_remove(const struct move *v, const struct solution *s)
{
    int i, n = s->ifacil[v->facil];
    double obj = 0.0;
    obj -= s->prob->flow[v->facil * s->prob->n + v->facil] * s->prob->dist[v->local * s->prob->n + v->local];
    for (i = 0; i < n; ++i)
        obj -= s->prob->flow[v->facil * s->prob->n + s->facil[i]] * s->prob->dist[v->local * s->prob->n + s->local[i]] + s->prob->flow[s->facil[i] * s->prob->n + v->facil] * s->prob->dist[s->local[i] * s->prob->n + v->local];
    for (i = n + 1; i < s->n_a; ++i)
        obj -= s->prob->flow[v->facil * s->prob->n + s->facil[i]] * s->prob->dist[v->local * s->prob->n + s->local[i]] + s->prob->flow[s->facil[i] * s->prob->n + v->facil] * s->prob->dist[s->local[i] * s->prob->n + v->local];
    return obj;
}

/*
 * Compute lower bound increment after forbidding an assignment in a solution
 */
static double lbi_forbid(const struct move *v, const struct solution *s)
{
    return 0.0;
}

/*
 * Compute lower bound increment after permitting an assignment in a solution
 */
static double lbi_permit(const struct move *v, const struct solution *s)
{
    return 0.0;
}

/*************************/
/* Problem instantiation */
/*************************/

/*
 * Quadratic assignment instantiation
 * Status: TENTATIVE
 * Notes:
 *   Needs further error checking
 */
struct problem *newProblem(const char *filename)
{
    int n, i, j;
    FILE *inFile;
    struct problem *p = NULL;
    inFile = fopen(filename, "r");
    if (inFile) {
        fscanf(inFile, "%d", &n);
        if (n > 2) {
            p = malloc(sizeof(struct problem));
            p->n = n;
            p->e = n * n;
            p->flow = malloc(sizeof(double) * p->e);
            p->dist = malloc(sizeof(double) * p->e);
            for (i = 0; i < n; ++i)
                for (j = 0; j < n; ++j)
                    fscanf(inFile, "%lf", &p->flow[i * n + j]);
            for (i = 0; i < n; ++i)
                for (j = 0; j < n; ++j)
                    fscanf(inFile, "%lf", &p->dist[i * n + j]);
        }
        else
            fprintf(stderr, "Invalid quadratic assignment instance %s.\n", filename);
        fclose(inFile);
    }
    else
        fprintf(stderr, "Cannot open file %s.\n", filename);
    return p;
}

/*
 * Return the number of possible assignments between all facilities and
 * locations
 * Status: FINAL
 */
long getNumComponents(const struct problem *p)
{
    return p->e;
}

/*
 * Return the number of assignments between all facilities and different
 * locations
 * Status: FINAL
 */
long getMaxSolutionSize(const struct problem *p)
{
    return p->n;
}

/*
 * Return the largest possible number of neighbours in a given subneighbourhood
 * Status: TENTATIVE
 * Notes:
 *   Redefine the concept of 'neighbourhood' and 'subneighbourhood'
 *   Should handle unimplemented exceptions
 */
long getMaxNeighbourhoodSize(const struct problem *p, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return p->e;
    case REMOVE:
        return p->n;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for getMaxNeighbourhoodSize().\n");
        break;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for getMaxNeighbourhoodSize().\n");
        break;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to getMaxNeighbourhoodSize().\n");
        break;
    }
    return -1;
}

/*********************/
/* Memory management */
/*********************/

/*
 * Allocate memory for a solution
 * Status: CHECK
 */
struct solution *allocSolution(struct problem *p)
{
    int i, j, a, b;
    struct solution *s = malloc(sizeof(struct solution));
    s->prob = p;
    s->n_g = p->e;
    s->groundSet = malloc(sizeof(int) * p->e);
    s->igroundSet = malloc(sizeof(int) * p->e);
    s->n = p->n;
    s->facil = malloc(sizeof(int) * p->n);
    s->ifacil = malloc(sizeof(int) * p->n);
    s->local = malloc(sizeof(int) * p->n);
    s->ilocal = malloc(sizeof(int) * p->n);
    /* enumMove() and randomMoveWOR() support */
    s->rndSample = malloc(sizeof(int) * sumsquares(p->n));
    for (i = 0; i < p->n; ++i) {
        a = sumsquares(i);
        b = sumsquares(i + 1);
        for (j = a; j < b; ++j)
            s->rndSample[j] = j - a;
    }
    return s;
}

/*
 * Allocate memory for a move
 * Status: CHECK
 */
struct move *allocMove(struct problem *p)
{
    struct move *v = malloc(sizeof(struct move));
    v->prob = p;
    return v;
}

/*
 * Free the memory used by a problem
 * Status: CHECK
 */
void freeProblem(struct problem *p)
{
    free(p->flow);
    free(p->dist);
    free(p);
}

/*
 * Free the memory used by a solution
 * Status: CHECK
 */
void freeSolution(struct solution *s)
{
    free(s->groundSet);
    free(s->igroundSet);
    free(s->local);
    free(s->ilocal);
    free(s->facil);
    free(s->ifacil);
    free(s->rndSample);
    free(s);
}

/*
 * Free the memory used by a move
 * Status: CHECK
 */
void freeMove(struct move *v)
{
    free(v);
}

/*************/
/* Reporting */
/*************/

/*
 * Print user-formatted representation of problem instance
 * Status: CHECK
 */
void printProblem(const struct problem *p)
{
    int i, j;
    printf("Quadratic Assignment (%d) problem\n  Flow:    ", p->n);
    for (i = 0; i < p->n; ++i) {
        for (j = 0; j < p->n; ++j)
            printf(" %*.0lf", 7, p->flow[i * p->n + j]);
        printf(i < p->n - 1 ? "\n\t   " : "");
    }
    printf("\n  Distance:");
    for (i = 0; i < p->n; ++i) {
        for (j = 0; j < p->n; ++j)
            printf(" %*.0lf", 7, p->dist[i * p->n + j]);
        printf(i < p->n - 1 ? "\n\t   " : "");
    }
    printf("\n\n");
}

/*
 * Print user-formatted representation of solution
 * Status: CHECK
 */
void printSolution(const struct solution *s)
{
    int i;
    printf("Quadratic Assignment (%d) solution\n  p:", s->n);
    for (i = 0; i < s->n; ++i)
        printf(" %d", s->ifacil[i] < s->n_a ? s->local[s->ifacil[i]] : -1);
#if __DEBUG__
    printf("\n  present  :");
    for (i = 0; i < s->n_p; ++i)
        printf(" %d", s->groundSet[i]);
    printf("\n  forbidden:");
    for (i = s->n_g - s->n_f; i < s->n_g; ++i)
        printf(" %d", s->groundSet[i]);
    printf("\n  undefined:");
    for (i = s->n_p; i < s->n_g - s->n_f; ++i)
        printf(" %d", s->groundSet[i]);
#endif
    if (s->n_p == s->n)
        printf("\n  solution is feasible and complete");
    else
        printf("\n  solution is unfeasible and partial");
    if (s->evalLB)
        printf("\n  lower bound: %*.1lf", 14, s->objLB);
    if (s->evalv)
        printf("\n  objective value: %*.1lf", 10, s->objv);
    printf("\n\n");
}

/*
 * Print user-formatted representation of move
 * Status: CHECK
 */
void printMove(const struct move *v)
{
    printf("Quadratic Assignment (%d) move: %d (%d, %d)\n\n", v->prob->n, v->data, v->facil, v->local);
}

/***********************/
/* Solution generation */
/***********************/

/*
 * Initialise empty solution
 * Status: CHECK
 */
struct solution *emptySolution(struct solution *s)
{
    /* solution s must have been allocated with allocSolution() */
    int i;
    for (i = 0; i < s->n_g; ++i)
        s->groundSet[i] = s->igroundSet[i] = i;
    s->n_p = s->n_f = 0;
    for (i = 0; i < s->n; ++i)
        s->facil[i] = s->ifacil[i] = s->local[i] = s->ilocal[i] = i;
    s->n_a = 0;
    /* solution not evaluated yet */
    s->evalv = s->evalLB = 0;
    s->numEnumLim = st_size(s);
    s->sampleWORLim = s->sampleEnumLim = s->sampleSize = nh_size(s);
    return s;
}

/*
 * Heuristically constructs a quadratic assignment solution using the suboptimal
 * procedure proposed on Gilmore, P. C. (1962). Optimal and Suboptimal
 * Algorithms for the Quadratic Assignment Problem. Journal of the Society for
 * Industrial and Applied Mathematics, 10(2), 305-313.
 * Status: INTERIM
 * Notes:
 *   Implement the suboptimal procedure proposed on Gilmore, P. C. (1962).
 * Optimal and Suboptimal Algorithms for the Quadratic Assignment Problem.
 * Journal of the Society for Industrial and Applied Mathematics, 10(2),
 * 305-313.
 */
struct solution *heuristicSolution(struct solution *s)
{
    if (s->n_a >= s->n) /* solution s is feasible, return it */
        return s;
    if (s->n_g - s->n_f < s->n) /* not enough edges to complete assignment, no feasible solution is possible */
        return NULL;
    fprintf(stderr, "heuristicSolution() not implemented yet for quadratic assignment problem.\n");
    /* update solution */
    s->evalv = s->evalLB = 0;
    s->numEnumLim = st_size(s);
    s->sampleWORLim = s->sampleEnumLim = s->sampleSize = nh_size(s);
    return s;
}

/***********************/
/* Solution inspection */
/***********************/

/*
 * Solution evaluation
 * Status: FINAL
 */
double *getObjectiveVector(double *objv, struct solution *s)
{
    if (s->n_a < s->n) /* solution is unfeasible, cannot evaluate it */
        return NULL;
    int i, j;
    double obj = 0.0;
    if (s->evalv) /* solution s is evaluated */
        *objv = s->objv;
    else { /* solution s is not evaluated */
        for (i = 0; i < s->n; ++i)
            for (j = 0; j < s->n; ++j)
                obj += s->prob->flow[s->facil[i] * s->prob->n + s->facil[j]] * s->prob->dist[s->local[i] * s->prob->n + s->local[j]];
        *objv = s->objv = obj;
        s->evalv = 1;
    }
    return objv;
}

/*
 * Lower bound evaluation
 * Status: INTERIM
 * Notes:
 *   Implement a tighter lower bound for quadratic assignment solution
 */
double *getObjectiveLB(double *objLB, struct solution *s)
{
    if (s->evalLB) /* solution s is evaluated */
        *objLB = s->objLB;
    else { /* solution s is not evaluated */
        *objLB = s->objLB = lb(s);
        s->evalLB = 1;
    }
    return objLB;
}

/*
 * Return the number of components of a solution that are in a given state
 * Status: FINAL
 */
long getNumSolutionComponents(const struct solution *s, const enum ComponentState st)
{
    switch (st) {
    case PRESENT:
        return s->n_p;
    case FORBIDDEN:
        return s->n_f;
    case UNDEFINED:
        return s->n_g - (s->n_p + s->n_f);
    default:
        fprintf(stderr, "Invalid state passed to getNumSolutionComponents().\n");
        break;
    }
    return -1;
}

/*
 * Return true if all facilities are assigned to different locations or false if
 * there are facilities that are not assigned to a location
 * Status: FINAL
 */
int isFeasible(struct solution *s)
{
    return s->n_p == s->n;
}

/*
 * Return true if no facilities need to be assigned to a location or false if
 * there are facilities that need to be assigned to a location
 * Status: FINAL
 */
int isComplete(struct solution *s)
{
    return s->n_p == s->n;
}

/*
 * Return the number of neighbours in a given subneighbouhood of a solution
 * Status: TENTATIVE
 * Notes:
 *   Redefine the concept of 'neighbourhood' and 'subneighbourhood'
 *   Should handle unimplemented exceptions
 */
long getNeighbourhoodSize(struct solution *s, const enum SubNeighbourhood nh)
{
    int uplim = s->n_g - s->n_f, n = s->n - s->n_a, size, i, j;
    switch (nh) {
    case ADD:
        for (i = s->n_a, size = n * n; i < s->n; ++i)
                for (j = s->n_a; j < s->n; ++j)
                    if (s->igroundSet[pairing(s->facil[i], s->local[j], s->n)] >= uplim) /* assignment is forbidden, decrement size */
                        size--;
        return size;
    case REMOVE:
        return s->n_a;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for getNeighbourhoodSize().\n");
        break;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for getNeighbourhoodSize().\n");
        break;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to getNeighbourhoodSize().\n");
        break;
    }
    return -1;
}

/*
 * Enumerate the components of a solution that are in a given state
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
long enumSolutionComponents(struct solution *s, const enum ComponentState st)
{
    switch (st) {
    case PRESENT:
        if (s->numEnumLim.present <= 0)
            return -1;
        return s->groundSet[s->n_p - s->numEnumLim.present--];
    case FORBIDDEN:
        fprintf(stderr, "Forbidden state not implemented for enumSolutionComponents().\n");
        break;
    case UNDEFINED:
        fprintf(stderr, "Undefined state not implemented for enumSolutionComponents().\n");
        break;
    default:
        fprintf(stderr, "Invalid state passed to enumSolutionComponents().\n");
        break;
    }
    return -1;
}

/*
 * Reset the enumeration of the components of a solution that are in a given
 * state
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
struct solution *resetEnumSolutionComponents(struct solution *s, const enum ComponentState st)
{
    switch (st) {
    case PRESENT:
        s->numEnumLim.present = s->n_p;
        return s;
    case FORBIDDEN:
        fprintf(stderr, "Forbidden state not implemented for resetEnumSolutionComponents().\n");
        break;
    case UNDEFINED:
        fprintf(stderr, "Undefined state not implemented for resetEnumSolutionComponents().\n");
        break;
    default:
        fprintf(stderr, "Invalid state passed to resetEnumSolutionComponents().\n");
        break;
    }
    return NULL;
}

/*******************/
/* Move generation */
/*******************/

/*
 * Reset the enumeration of a given subneighbourhood of a solution
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
struct solution *resetEnumMove(struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        s->sampleEnumLim.add = s->sampleSize.add;
        return s;
    case REMOVE:
        s->sampleEnumLim.remove = s->n_p;
        return s;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for resetEnumMove().\n");
        break;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for resetEnumMove().\n");
        break;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to resetEnumMove().\n");
        break;
    }
    return NULL;
}

/*
 * Reset the uniform random sampling without replacement of a given
 * subneighbourhood of a solution
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
struct solution *resetRandomMoveWOR(struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        s->sampleWORLim.add = s->sampleSize.add;
        return s;
    case REMOVE:
        s->sampleWORLim.remove = s->n_p;
        return s;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for resetRandomMoveWOR().\n");
        break;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for resetRandomMoveWOR().\n");
        break;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to resetRandomMoveWOR().\n");
        break;
    }
    return NULL;
}

static struct move *randomAddMove(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available (s->sampleSize.add is zero)
    //               - Complete solutions
    //               - All valid components are forbidden
    //           - Not enough components are available to complete solution
    if (!s->sampleSize.add || s->n_g - s->n_f < s->n)
        return NULL;

    // Sample a random valid component:
    //           - Cannot be forbidden
    //           - Available forbidden components cannot be resampled
    //           - Do not generate move if all available components are forbidden
    int i, j, k, x, y, z, idx, q = sumsquares(s->n - s->n_a - 1);

    for (;;) {
        // Sample new assignment
        idx = randint(s->sampleSize.add - 1);
        z = s->rndSample[q + idx];
        unpairing(&x, &y, z, s->n - s->n_a);
        i = s->facil[s->n_a + x];
        j = s->local[s->n_a + y];
        k = pairing(i, j, s->n);

        // Valid move if component is not forbidden
        if (s->igroundSet[k] < s->n_g - s->n_f)
            break;

        // Decrease number of valid moves available
        s->sampleSize.add--;

        // No valid moves are available
        if (!s->sampleSize.add)
            return NULL;

        // Exclude if component is forbidden
        swap(s->rndSample + q, idx, s->sampleSize.add);
    }
    v->facil = i;
    v->local = j;
    v->data = k;

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 4 * sizeof(int));
    return v;
}

static struct move *randomRemoveMove(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available (s->sampleSize.remove is zero)
    //               - Empty solutions
    if (!s->sampleSize.remove)
        return NULL;

    // Sample a random valid component
    int idx = randint(s->sampleSize.remove - 1);

    // Get facility, location, and component identifier
    v->facil = s->facil[idx];
    v->local = s->local[idx];
    v->data = pairing(v->facil, v->local, s->n);

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 4 * sizeof(int));
    return v;
}

struct move *randomMove(struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return randomAddMove(v, s);

    case REMOVE:
        return randomRemoveMove(v, s);

    case FORBID:
        return NULL; // TODO: Use unimplemented() function.

    case PERMIT:
        return NULL; // TODO: Use unimplemented() function.

    default:
        return NULL; // TODO: Use unimplemented() function.
    }
}

static struct move *randomAddMoveWOR(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available (s->sampleWORLim.add is zero)
    //               - Complete solutions
    //               - All valid components are forbidden or used
    //           - Not enough components are available to complete solution
    if (!s->sampleWORLim.add || s->n_g - s->n_f < s->n)
        return NULL;

    // Sample a random valid component:
    //           - Cannot be forbidden
    //           - Available forbidden and used components cannot be resampled
    //           - Do not generate move if all available components are forbidden
    int i, j, k, x, y, z, idx, q = sumsquares(s->n - s->n_a - 1);

    for (;;) {
        // Sample new assignment
        idx = (s->sampleSize.add - s->sampleWORLim.add) + randint(s->sampleWORLim.add - 1);
        z = s->rndSample[q + idx];
        unpairing(&x, &y, z, s->n - s->n_a);
        i = s->facil[s->n_a + x];
        j = s->local[s->n_a + y];
        k = pairing(i, j, s->n);

        // Valid move if component is not forbidden
        if (s->igroundSet[k] < s->n_g - s->n_f) {
            // Exclude component from resampling
            swap(s->rndSample + q, s->sampleSize.add - s->sampleWORLim.add, idx);

            // Decrease number of unused moves available
            s->sampleWORLim.add--;
            break;
        }

        // Decrease number of valid (and unused) moves available
        s->sampleSize.add--;
        s->sampleWORLim.add--;

        // No valid moves are available
        if (!s->sampleSize.add)
            return NULL;

        // Exclude if component is forbidden
        swap(s->rndSample + q, idx, s->sampleSize.add);
    }
    v->facil = i;
    v->local = j;
    v->data = k;

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 4 * sizeof(int));
    return v;
}

static struct move *randomRemoveMoveWOR(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available (s->sampleWORLim.remove is zero)
    //               - Empty solutions
    //               - All valid components are used
    if (!s->sampleWORLim.remove)
        return NULL;

    // Sample a random valid component:
    //           - Available used components cannot be resampled
    int idx = randint(s->sampleWORLim.remove - 1);

    // Get facility, location, and component identifier
    v->facil = s->facil[idx];
    v->local = s->local[idx];
    v->data = pairing(v->facil, v->local, s->n);

    // Exclude component from resampling
    swap_i(s->facil, s->ifacil, s->sampleWORLim.remove - 1, idx);
    swap_i(s->local, s->ilocal, s->sampleWORLim.remove - 1, idx);

    // Decrease number of unused moves available
    s->sampleWORLim.remove--;

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 4 * sizeof(int));
    return v;
}

struct move *randomMoveWOR(struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return randomAddMoveWOR(v, s);

    case REMOVE:
        return randomRemoveMoveWOR(v, s);

    case FORBID:
        return NULL; // TODO: Use unimplemented() function.

    case PERMIT:
        return NULL; // TODO: Use unimplemented() function.

    default:
        return NULL; // TODO: Use unimplemented() function.
    }
}

static struct move *enumAddMove(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available (s->sampleEnumLim.add is zero)
    //               - Complete solutions
    //               - All valid components are forbidden or used
    //           - Not enough components are available to complete solution
    if (!s->sampleEnumLim.add || s->n_g - s->n_f < s->n)
        return NULL;

    // Sample a random valid component:
    //           - Cannot be forbidden
    //           - Available forbidden and used components cannot be resampled
    //           - Do not generate move if all available components are forbidden
    int i, j, k, x, y, z, idx, q = sumsquares(s->n - s->n_a - 1);

    for (;;) {
        // Sample new assignment
        idx = s->sampleSize.add - s->sampleEnumLim.add;
        z = s->rndSample[q + idx];
        unpairing(&x, &y, z, s->n - s->n_a);
        i = s->facil[s->n_a + x];
        j = s->local[s->n_a + y];
        k = pairing(i, j, s->n);

        // Valid move if component is not forbidden
        if (s->igroundSet[k] < s->n_g - s->n_f) {
            // Decrease number of unused moves available
            s->sampleEnumLim.add--;
            break;
        }

        // Decrease number of valid (and unused) moves available
        s->sampleSize.add--;
        s->sampleEnumLim.add--;

        // No valid moves are available
        if (!s->sampleSize.add)
            return NULL;

        // Exclude if component is forbidden
        swap(s->rndSample + q, idx, s->sampleSize.add);
    }
    v->facil = i;
    v->local = j;
    v->data = k;

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 4 * sizeof(int));
    return v;
}

static struct move *enumRemoveMove(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available (s->sampleEnumLim.remove is zero)
    //               - Empty solutions
    //               - All valid components used
    if (!s->sampleEnumLim.remove)
        return NULL;

    // Sample a random valid component as for a given assignment:
    //           - Available used components cannot be resampled
    int idx = s->sampleSize.remove - s->sampleEnumLim.remove;

    // Get facility, location, and component identifier
    v->facil = s->facil[idx];
    v->local = s->local[idx];
    v->data = pairing(v->facil, v->local, s->n);

    // Decrease number of unused moves available
    s->sampleEnumLim.remove--;

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 4 * sizeof(int));
    return v;
}

struct move *enumMove(struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return enumAddMove(v, s);

    case REMOVE:
        return enumRemoveMove(v, s);

    case FORBID:
        return NULL; // TODO: Use unimplemented() function.

    case PERMIT:
        return NULL; // TODO: Use unimplemented() function.

    default:
        return NULL; // TODO: Use unimplemented() function.
    }
}

/**************************/
/* Operations on solutions*/
/**************************/

/*
 * Copy the arguments of the second argument to the first argument
 * Status: CHECK
 */
struct solution *copySolution(struct solution *dest, const struct solution *src)
{
    dest->prob = src->prob;
    dest->n_g = src->n_g;
    memcpy(dest->groundSet, src->groundSet, src->n_g * sizeof(int));
    memcpy(dest->igroundSet, src->igroundSet, src->n_g * sizeof(int));
    dest->n_p = src->n_p;
    dest->n_f = src->n_f;
    dest->n = src->n;
    memcpy(dest->facil, src->facil, src->n * sizeof(int));
    memcpy(dest->ifacil, src->ifacil, src->n * sizeof(int));
    memcpy(dest->local, src->local, src->n * sizeof(int));
    memcpy(dest->ilocal, src->ilocal, src->n * sizeof(int));
    dest->n_a = src->n_a;
    memcpy(dest->rndSample, src->rndSample, sumsquares(src->n) * sizeof(int));
    dest->numEnumLim = src->numEnumLim;
    dest->sampleSize = src->sampleSize;
    dest->sampleEnumLim = src->sampleEnumLim;
    dest->sampleWORLim = src->sampleWORLim;
    dest->evalv = src->evalv;
    dest->objv = src->objv;
    dest->evalLB = src->evalLB;
    dest->objLB = src->objLB;
    return dest;
}

static struct solution *applyAddMove(struct solution *s, const struct move *v)
{
    int k = s->igroundSet[v->data];
    int x = s->ifacil[v->facil];
    int y = s->ilocal[v->local];

    // Coverage: - Component is not undefined
    //           - Solution is complete
    //           - Facility x is already in use
    //           - Location y is already in use
    //
    // The || operator guarantees left-to-right evaluation.
    // If the first operand compares unequal to 0, the second operand is not evaluated.
#if __DEBUG__
    if (k < s->n_p || k >= s->n_g - s->n_f || s->n_p >= s->n || x < s->n_a || y < s->n_a) {
        fprintf(stderr, "Component %d is not undefined, solution is complete, facility %d or location %d is already in use.\nCould not add component.\n\n", s->groundSet[k], v->facil, v->local);
        return NULL;
    }
#endif

    // Update permutations mapping facilities and locations in regards to the new assignment
    swap_i(s->facil, s->ifacil, s->n_a, x);
    swap_i(s->local, s->ilocal, s->n_a, y);
    s->n_a++;

    // Add *undefined* component
    swap_i(s->groundSet, s->igroundSet, s->n_p, k);
    s->n_p++;

    // Update state of evaluation of objective value of incomplete (or incomplete) solution
    s->evalv = 0;

    // Update lower bound and position of the component with largest cost added in the lower bound calculation
    if (s->evalLB && v->evalLBi[0]) {
        s->objLB += v->objLBi;
    }
    else {
        s->evalLB = 0;
    }

    // Update number of components left to enumerate
    s->numEnumLim = st_size(s);

    // Update number of unused moves
    s->sampleWORLim = (s->sampleEnumLim = (s->sampleSize = nh_size(s)));

    return s;
}

static struct solution *applyForbidMove(struct solution *s, const struct move *v)
{
    int k = s->igroundSet[v->data];
    int lim = s->n_g - s->n_f - 1;

    // Coverage: - Component is not undefined
#if __DEBUG__
    if (k < s->n_p || k > lim) {
        fprintf(stderr, "Component %d is not undefined.\nCould not forbid component.\n\n", s->groundSet[k]);
        return NULL;
    }
#endif

    // Forbid *undefined* component
    swap_i(s->groundSet, s->igroundSet, k, lim);
    s->n_f++;

    // Update lower bound and position of the component with largest cost added in the lower bound calculation
    if (s->evalLB && v->evalLBi[2]) {
        s->objLB += v->objLBi;
    }
    else {
        s->evalLB = 0;
    }

    // Update number of components left to enumerate
    s->numEnumLim = st_size(s);

    // Update number of unused moves
    s->sampleWORLim = (s->sampleEnumLim = (s->sampleSize = nh_size(s)));

    return s;
}

static struct solution *applyRemoveMove(struct solution *s, const struct move *v)
{
    int k = s->igroundSet[v->data];
    int x = s->ifacil[v->facil];
    int y = s->ilocal[v->local];

    // Coverage: - Component is not present
    //           - Solution is empty
    //           - Facility x is not being used
    //           - Location y is not being used
    //
    // The || operator guarantees left-to-right evaluation.
    // If the first operand compares unequal to 0, the second operand is not evaluated.
#if __DEBUG__
    if (k >= s->n_p || s->n_p < 1 || x >= s->n_a || y >= s->n_a || x != y) {
        fprintf(stderr, "Component %d is not present, solution is empty, facility %d or location %d is not being used.\nCould not remove component.\n\n", s->groundSet[k], v->facil, v->local);
        return NULL;
    }
#endif

    // Update permutations mapping facilities and locations regarding the new assignment
    swap_i(s->facil, s->ifacil, x, s->n_a - 1);
    swap_i(s->local, s->ilocal, y, s->n_a - 1);
    s->n_a--;

    // Remove *present* component
    swap_i(s->groundSet, s->igroundSet, k, s->n_p - 1);
    s->n_p--;

    // Update state of evaluation of objective value of incomplete (or incomplete) solution
    s->evalv = 0;

    // Update lower bound and position of the component with largest cost added in the lower bound calculation
    if (s->evalLB && v->evalLBi[1]) {
        s->objLB += v->objLBi;
    }
    else {
        s->evalLB = 0;
    }

    // Update number of components left to enumerate
    s->numEnumLim = st_size(s);

    // Update number of unused moves
    s->sampleWORLim = (s->sampleEnumLim = (s->sampleSize = nh_size(s)));

    return s;
}

static struct solution *applyPermitMove(struct solution *s, const struct move *v)
{
    int k = s->igroundSet[v->data];
    int lim = s->n_g - s->n_f;

    // Coverage: - Component is not forbidden
#if __DEBUG__
    if (k < lim) {
        fprintf(stderr, "Component %d is not forbidden.\nCould not permit component.\n\n", s->groundSet[k]);
        return NULL;
    }
#endif

    // Permit *forbidden* component
    swap_i(s->groundSet, s->igroundSet, lim, k);
    s->n_f--;

    // Update lower bound and position of the component with largest cost added in the lower bound calculation
    if (s->evalLB && v->evalLBi[3]) {
        s->objLB += v->objLBi;
    }
    else {
        s->evalLB = 0;
    }

    // Update number of components left to enumerate
    s->numEnumLim = st_size(s);

    // Update number of unused moves
    s->sampleWORLim = (s->sampleEnumLim = (s->sampleSize = nh_size(s)));

    return s;
}

struct solution *applyMove(struct solution *s, const struct move *v, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return applyAddMove(s, v);

    case REMOVE:
        return applyRemoveMove(s, v);

    case FORBID:
        return applyForbidMove(s, v);

    case PERMIT:
        return applyPermitMove(s, v);

    default:
        return NULL; // TODO: Use unimplemented() function.
    }
}

/***********************/
/* Operations on moves */
/***********************/

/*
 * Copy the contents of the second argument to the first argument
 * Status: CHECK
 */
struct move *copyMove(struct move *dest, const struct move *src)
{
    dest->prob = src->prob;
    dest->data = src->data;
    dest->facil = src->facil;
    dest->local = src->local;
    memcpy(dest->evalLBi, src->evalLBi, 4 * sizeof(int));
    dest->objLBi = src->objLBi;
    return dest;
}

/*******************/
/* Move inspection */
/*******************/

/*
 * Return the unique component identifier with respect to a given move
 * Status: FINAL
 */
long getComponentFromMove(const struct move *v)
{
    return v->data;
}

/*
 * Move evaluation
 * Status: INTERIM
 * Notes:
 *   Redefine the concept of 'subneighbourhood'
 *   Implement a tighter lower bound for quadratic assignment solution
 */
double *getObjectiveLBIncrement(double *obji, struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    int i;
    switch (nh) {
    case ADD:
#if __DEBUG__
        if (s->ifacil[v->facil] < s->n_a || s->ilocal[v->local] < s->n_a || s->igroundSet[v->data] >= s->n_g - s->n_f) { /* assignment does not belong to Add neighbourhood */
            fprintf(stderr, "Invalid move passed to getObjectiveLBIncrement(). It does not belong to Add neighbourhood.\n");
            return NULL;
        }
#endif
        i = 0;
        break;
    case REMOVE:
#if __DEBUG__
        if (s->ifacil[v->facil] >= s->n_a || s->ilocal[v->local] >= s->n_a || s->ifacil[v->facil] != s->ilocal[v->local]) { /* assignment does not belong to Remove neighbourhood */
            fprintf(stderr, "Invalid move passed to getObjectiveLBIncrement(). It does not belong to Remove neighbourhood.\n");
            return NULL;
        }
#endif
        i = 1;
        break;
    case FORBID:
#if __DEBUG__
        if (s->igroundSet[v->data] < s->n_p || s->igroundSet[v->data] >= s->n_g - s->n_f) { /* item is not undefined */
            fprintf(stderr, "Invalid move passed to getObjectiveLBIncrement(). It does not belong to Forbid neighbourhood.\n");
            return NULL;
        }
#endif
        i = 2;
        break;
    case PERMIT:
#if __DEBUG__
        if (s->igroundSet[v->data] < s->n_g - s->n_f) { /* item is not forbidden */
            fprintf(stderr, "Invalid move passed to getObjectiveLBIncrement(). It does not belong to Permit neighbourhood.\n");
            return NULL;
        }
#endif
        i = 3;
        break;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to getObjectiveLBIncrement().\n");
        return NULL;
    }
    if (v->evalLBi[i]) /* move v is evaluated */
        *obji = v->objLBi;
    else { /* move v is not evaluated */
        memset(v->evalLBi, 0, sizeof(int) * 4);
        switch (nh) {
        case ADD:
            *obji = v->objLBi = lbi_add(v, s);
            break;
        case REMOVE:
            *obji = v->objLBi = lbi_remove(v, s);
            break;
        case FORBID:
            *obji = v->objLBi = lbi_forbid(v, s);
            break;
        case PERMIT:
            *obji = v->objLBi = lbi_permit(v, s);
            break;
        default:
            *obji = v->objLBi = DBL_MAX;
            break;
        }
        v->evalLBi[i] = 1;
    }
    return obji;
}
