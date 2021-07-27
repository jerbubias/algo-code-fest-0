#include "qmkp.h"
#include <float.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
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
    int n;            /* Number of items */
    int m;            /* Number of constraints */
    double *profit;   /* Profit matrix */
    double sum;       /* Sum of (nonnegative) profit of items */
    double *cost;     /* Cost array */
    double *scost;    /* Cost array regarding the surrogate constraint */
    double *capacity; /* Capacity */
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
    double *amount;                        /* Amount occupied by the items in the knapsack */
    int n_v;                               /* Number of valid components */
    struct sample_size sampleWORLim;       /* Number of unused moves for random move without replacement */
    struct sample_size sampleEnumLim;      /* Number of unused moves for enumerate move */
    struct sample_size sampleHeuristicLim; /* Number of unused moves for heuristic move */
    // Heuristic Construction
    double *densityLB; /* Density array used in the lower bound calculation */
};

struct move {
    struct problem *prob;
    int data;       /* Component identifier */
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
    size.add = s->n_v;
    size.remove = s->n_p;
    return size;
}

/*
 * Compute upper bound
 */
static double ub(struct solution *s)
{
    int lim = s->n_p + s->n_v, i, j;
    double obj = 0.0, p;
    /* compute profit of solution s */
    if (s->evalv) /* solution s is evaluated */
        obj += s->prob->sum - s->objv;
    else { /* solution s is not evaluated */
        for (i = 0; i < s->n_p; ++i)
            for (j = 0; j < s->n_p; ++j)
                obj += s->prob->profit[s->groundSet[i] * s->prob->n + s->groundSet[j]];
        s->objv = s->prob->sum - obj;
        s->evalv = 1;
    }
    /* compute the remaining upper bound */
    for (i = s->n_p; i < lim; ++i) {
        /* add positive profit of valid item */
        p = s->prob->profit[s->groundSet[i] * s->prob->n + s->groundSet[i]];
        obj += p > 0.0 ? p : 0.0;
        for (j = 0; j < i; ++j) {
            p = s->prob->profit[s->groundSet[i] * s->prob->n + s->groundSet[j]];
            obj += p > 0.0 ? p : 0.0;
            p = s->prob->profit[s->groundSet[j] * s->prob->n + s->groundSet[i]];
            obj += p > 0.0 ? p : 0.0;
        }
    }
    return obj;
}

/*
 * Compute lower bound increment after adding an item to a solution
 */
static double lbi_add(const struct move *v, const struct solution *s)
{
    int lim = s->n_p + s->n_v, n = s->igroundSet[v->data], i, j, k;
    double obj = 0.0, p;
    /* add positive profit of item */
    p = s->prob->profit[v->data * s->prob->n + v->data];
    obj -= p < 0.0 ? p : 0.0;
    for (i = 0; i < s->n_p; ++i) {
        p = s->prob->profit[s->groundSet[i] * s->prob->n + v->data];
        obj -= p < 0.0 ? p : 0.0;
        p = s->prob->profit[v->data * s->prob->n + s->groundSet[i]];
        obj -= p < 0.0 ? p : 0.0;
    }
    swap_i(s->groundSet, s->igroundSet, s->n_p, n);
    for (i = s->n_p + 1; i < lim; ++i)
        for (j = 0; j < s->prob->m; ++j)
            if (s->amount[j] + s->prob->cost[j * s->prob->n + v->data] + s->prob->cost[j * s->prob->n + s->groundSet[i]] > s->prob->capacity[j]) { /* item is invalid */
                /* remove positive profit of valid item */
                p = s->prob->profit[s->groundSet[i] * s->prob->n + s->groundSet[i]];
                obj += p > 0.0 ? p : 0.0;
                for (k = 0; k < i; ++k) {
                    p = s->prob->profit[s->groundSet[k] * s->prob->n + s->groundSet[i]];
                    obj += p > 0.0 ? p : 0.0;
                    p = s->prob->profit[s->groundSet[i] * s->prob->n + s->groundSet[k]];
                    obj += p > 0.0 ? p : 0.0;
                }
                for (k = i + 1; k < lim; ++k) {
                    p = s->prob->profit[s->groundSet[k] * s->prob->n + s->groundSet[i]];
                    obj += p > 0.0 ? p : 0.0;
                    p = s->prob->profit[s->groundSet[i] * s->prob->n + s->groundSet[k]];
                    obj += p > 0.0 ? p : 0.0;
                }
                swap_i(s->groundSet, s->igroundSet, i--, --lim);
                break;
            }
    return obj;
}

/*
 * Compute lower bound increment after removing an item from a solution
 */
static double lbi_remove(const struct move *v, const struct solution *s)
{
    int lim1 = s->n_p + s->n_v, lim2 = s->n_g - s->n_f, n = s->igroundSet[v->data], i, j, k;
    double obj = 0.0, p;
    /* remove positive profit of item */
    p = s->prob->profit[v->data * s->prob->n + v->data];
    obj += p < 0.0 ? p : 0.0;
    for (i = 0; i < n; ++i) {
        p = s->prob->profit[s->groundSet[i] * s->prob->n + v->data];
        obj += p < 0.0 ? p : 0.0;
        p = s->prob->profit[v->data * s->prob->n + s->groundSet[i]];
        obj += p < 0.0 ? p : 0.0;
    }
    for (i = n + 1; i < s->n_p; ++i) {
        p = s->prob->profit[s->groundSet[i] * s->prob->n + v->data];
        obj += p < 0.0 ? p : 0.0;
        p = s->prob->profit[v->data * s->prob->n + s->groundSet[i]];
        obj += p < 0.0 ? p : 0.0;
    }
    for (i = lim1; i < lim2; ++i) {
        for (j = 0; j < s->prob->m; ++j)
            if (s->amount[j] - s->prob->cost[j * s->prob->n + v->data] + s->prob->cost[j * s->prob->n + s->groundSet[i]] > s->prob->capacity[j])
                break;
        if (j == s->prob->m) { /* item is valid */
            /* add positive profit of valid item */
            p = s->prob->profit[s->groundSet[i] * s->prob->n + s->groundSet[i]];
            obj -= p > 0.0 ? p : 0.0;
            for (k = 0; k < lim1; ++k) {
                p = s->prob->profit[s->groundSet[k] * s->prob->n + s->groundSet[i]];
                obj -= p > 0.0 ? p : 0.0;
                p = s->prob->profit[s->groundSet[i] * s->prob->n + s->groundSet[k]];
                obj -= p > 0.0 ? p : 0.0;
            }
            swap_i(s->groundSet, s->igroundSet, i, lim1++);
        }
    }
    return obj;
}

/*
 * Compute lower bound increment after forbidding an item in a solution
 */
static double lbi_forbid(const struct move *v, const struct solution *s)
{
    int lim = s->n_p + s->n_v, n = s->igroundSet[v->data], i;
    double obj = 0.0, p;
    if (n < lim) { /* item is valid */
        /* remove positive profit of item */
        p = s->prob->profit[v->data * s->prob->n + v->data];
        obj += p > 0.0 ? p : 0.0;
        for (i = 0; i < n; ++i) {
            p = s->prob->profit[v->data * s->prob->n + s->groundSet[i]];
            obj += p > 0.0 ? p : 0.0;
            p = s->prob->profit[s->groundSet[i] * s->prob->n + v->data];
            obj += p > 0.0 ? p : 0.0;
        }
        for (i = n + 1; i < lim; ++i) {
            p = s->prob->profit[v->data * s->prob->n + s->groundSet[i]];
            obj += p > 0.0 ? p : 0.0;
            p = s->prob->profit[s->groundSet[i] * s->prob->n + v->data];
            obj += p > 0.0 ? p : 0.0;
        }
    }
    return obj;
}

/*
 * Compute lower bound increment after permitting an item in a solution
 */
static double lbi_permit(const struct move *v, const struct solution *s)
{
    int lim = s->n_p + s->n_v, i, j;
    double obj = 0.0, p;
    for (i = 0; i < s->prob->m; ++i)
        if (s->amount[i] + s->prob->cost[i * s->prob->n + v->data] > s->prob->capacity[i])
            break;
    if (i == s->prob->m) { /* item is valid */
        /* add positive profit of item */
        p = s->prob->profit[v->data * s->prob->n + v->data];
        obj -= p > 0.0 ? p : 0.0;
        for (j = 0; j < lim; ++j) {
            p = s->prob->profit[v->data * s->prob->n + s->groundSet[j]];
            obj -= p > 0.0 ? p : 0.0;
            p = s->prob->profit[s->groundSet[j] * s->prob->n + v->data];
            obj -= p > 0.0 ? p : 0.0;
        }
    }
    return obj;
}

/*************************/
/* Problem instantiation */
/*************************/

/*
 * Quadratic multidimensional knapsack instantiation
 * Status: TENTATIVE
 * Notes:
 *   Needs further error checking
 */
struct problem *newProblem(const char *filename)
{
    int n, m, i, j;
    FILE *inFile;
    struct problem *p = NULL;
    inFile = fopen(filename, "r");
    if (inFile) {
        fscanf(inFile, "%d %d", &n, &m);
        if (n > 0 && m > 0) {
            p = malloc(sizeof(struct problem));
            p->n = n;
            p->m = m;
            p->profit = malloc(sizeof(double) * n * n);
            p->cost = malloc(sizeof(double) * m * n);
            p->scost = malloc(sizeof(double) * n);
            p->capacity = malloc(sizeof(double) * m);
            for (i = 0, p->sum = 0.0; i < n; ++i)
                for (j = 0; j < n; ++j) {
                    fscanf(inFile, "%lf", &p->profit[i * n + j]);
                    p->sum += p->profit[i * n + j] > 0.0 ? p->profit[i * n + j] : 0.0;
                }
            for (i = 0; i < m; ++i)
                for (j = 0; j < n; ++j)
                    fscanf(inFile, "%lf", &p->cost[i * n + j]);
            for (i = 0; i < n; ++i) {
                p->scost[i] = 0.0;
                for (j = 0; j < m; ++j)
                    p->scost[i] += p->cost[j * n + i];
            }
            for (i = 0; i < m; ++i)
                fscanf(inFile, "%lf", &p->capacity[i]);
        }
        else
            fprintf(stderr, "Invalid quadratic multidimensional knapsack instance %s.\n", filename);
        fclose(inFile);
    }
    else
        fprintf(stderr, "Cannot open file %s.\n", filename);
    return p;
}

/*
 * Return the number of items
 * Status: FINAL
 */
long getNumComponents(const struct problem *p)
{
    return p->n;
}

/*
 * Return the number of items
 * Status: FINAL
 */
long getMaxSolutionSize(const struct problem *p)
{
    return p->n;
}

/*
 * Return the number of items as the largest possible number of neighbours in a
 * given subneighbourhood
 * Status: TENTATIVE
 * Notes:
 *   Redefine the concept of 'neighbourhood' and 'subneighbourhood'
 *   Should handle unimplemented exceptions
 */
long getMaxNeighbourhoodSize(const struct problem *p, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return p->n;
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
    struct solution *s = malloc(sizeof(struct solution));
    s->prob = p;
    s->n_g = p->n;
    s->groundSet = malloc(sizeof(int) * p->n);
    s->igroundSet = malloc(sizeof(int) * p->n);
    s->amount = malloc(sizeof(double) * p->m);
    /* heuristicSolution() support */
    s->densityLB = malloc(sizeof(double) * p->n);
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
    free(p->profit);
    free(p->cost);
    free(p->scost);
    free(p->capacity);
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
    free(s->amount);
    free(s->densityLB);
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
    printf("Quadratic Multidimensional Knapsack (%d) problem\n  Profit:  ", p->n);
    for (i = 0; i < p->n; ++i) {
        for (j = 0; j < p->n; ++j)
            printf(" %*.1lf", 9, p->profit[i * p->n + j]);
        printf(i < p->n - 1 ? "\n\t   " : "");
    }
    printf("\n  Weight:  ");
    for (i = 0; i < p->m; ++i) {
        for (j = 0; j < p->n; ++j)
            printf(" %*.1lf", 9, p->cost[i * p->n + j]);
        printf(i < p->m - 1 ? " %*.1lf\n\t   " : " %*.1lf", 9, p->capacity[i]);
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
    printf("Quadratic Multidimensional Knapsack (%d) solution\n  x:", s->prob->n);
    for (i = 0; i < s->n_p; ++i)
        printf(" %d", s->groundSet[i]);
#if __DEBUG__
    printf("\n  present  :");
    for (i = 0; i < s->n_p; ++i)
        printf(" %d", s->groundSet[i]);
    printf("\n  forbidden:");
    for (i = s->n_g - s->n_f; i < s->n_g; ++i)
        printf(" %d", s->groundSet[i]);
    printf("\n  valid    :");
    for (i = s->n_p; i < s->n_p + s->n_v; ++i)
        printf(" %d", s->groundSet[i]);
    printf("\n  invalid  :");
    for (i = s->n_p + s->n_v; i < s->n_g - s->n_f; ++i)
        printf(" %d", s->groundSet[i]);
#endif
    if (!s->n_v)
        printf("\n  solution is feasible and complete");
    else
        printf("\n  solution is feasible and partial");
    if (s->evalLB)
        printf("\n  upper bound: %*.1lf", 14, s->prob->sum - s->objLB);
    if (s->evalv)
        printf("\n  objective value: %*.1lf", 10, s->prob->sum - s->objv);
#if __DEBUG__
    printf("\n  constraint(s):");
    for (i = 0; i < s->prob->m; ++i)
        printf(i < s->prob->m - 1 ? "    %*.1lf <= %*.1lf\n\t\t" : "    %*.1lf <= %*.1lf", 9, s->amount[i], 9, s->prob->capacity[i]);
#endif
    printf("\n\n");
}

/*
 * Print user-formatted representation of move
 * Status: CHECK
 */
void printMove(const struct move *v)
{
    printf("Quadratic Multidimensional Knapsack (%d) move: %d\n\n", v->prob->n, v->data);
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
    memset(s->amount, 0, s->prob->m * sizeof(double));
    s->n_v = s->n_g;
    /* solution not evaluated yet */
    s->evalv = s->evalLB = 0;
    s->numEnumLim = st_size(s);
    s->sampleWORLim = s->sampleEnumLim = nh_size(s);
    return s;
}

/*
 * Heuristically constructs a quadratic multidimensional knapsack solution using
 * the algorithm proposed in Billionnet, A., & Calmels, F. (1996). Linear
 * Programming for the 0–1 Quadratic Knapsack Problem. European Journal of
 * Operational Research, 92(2), 310-325.
 * Status: TENTATIVE
 * Notes:
 *   Implement fill-up and exchange procedure as proposed in Gallo, G., Hammer,
 * P. L., & Simeone, B. (1980). Quadratic Knapsack Problems. In Combinatorial
 * Optimization (pp. 132-149). Springer, Berlin, Heidelberg. This procedure is
 * used in the second stage of the algorithm
 */
struct solution *heuristicSolution(struct solution *s)
{
    if (!s->n_v) /* solution s is complete, return it */
        return s;
    int oldlim = s->n_p + s->n_v, newlim = oldlim, z, item = -1, i, j;
    double mindensity = DBL_MAX;
    /* Heuristic method proposed in Chaillou, P., Hansen, P., & Mahieu, Y. (1989).
     * Best Network Flow Bounds for the Quadratic Knapsack Problem. In Combinatorial
     * Optimization (pp. 225-235). Springer, Berlin, Heidelberg. This method is used
     * in the first stage of the algorithm
     */
    /* Add valid items into knapsack. Compute density of each item regarding all
     * items in the knapsack.
     */
    for (i = s->n_p; i < oldlim; ++i) {
        z = s->groundSet[i]; /* undefined item */
        /* add item into knapsack */
        for (j = 0; j < s->prob->m; ++j)
            s->amount[j] += s->prob->cost[j * s->prob->n + z];
        /* compute density of item */
        s->densityLB[z] = s->prob->profit[z * s->prob->n + z];
        for (j = 0; j < i; ++j)
            s->densityLB[z] += s->prob->profit[s->groundSet[j] * s->prob->n + z] + s->prob->profit[z * s->prob->n + s->groundSet[j]];
        for (j = i + 1; j < oldlim; ++j)
            s->densityLB[z] += s->prob->profit[s->groundSet[j] * s->prob->n + z] + s->prob->profit[z * s->prob->n + s->groundSet[j]];
        s->densityLB[z] = s->densityLB[z] / s->prob->scost[z];
        if (s->densityLB[z] < mindensity) {
            item = z; /* item with lowest density */
            mindensity = s->densityLB[z];
        }
    }
    /* Remove item with lowest density until a feasible solution is obtained and the
     * lowest density is nonnegative.
     */
    for (i = 0; i < s->prob->m; ++i)
        if (s->amount[i] > s->prob->capacity[i]) /* solution is unfeasible */
            break;
    while (i < s->prob->m || mindensity < 0.0) {
        /* remove item with lowest density */
        swap_i(s->groundSet, s->igroundSet, s->igroundSet[item], --newlim);
        for (j = 0; j < s->prob->m; ++j)
            s->amount[j] -= s->prob->cost[j * s->prob->n + item];
        for (j = s->n_p, mindensity = DBL_MAX; j < newlim; ++j) {
            z = s->groundSet[j]; /* undefined item */
            /* update density of item */
            s->densityLB[z] -= (s->prob->profit[item * s->prob->n + z] + s->prob->profit[z * s->prob->n + item]) / s->prob->scost[z];
            if (s->densityLB[z] < mindensity) {
                item = z; /* new item with lowest density */
                mindensity = s->densityLB[z];
            }
        }
        for (i = 0; i < s->prob->m; ++i)
            if (s->amount[i] > s->prob->capacity[i])
                break;
    }
    /* update solution */
    s->n_p = newlim;
    for (i = oldlim - 1; i >= newlim; --i) {
        z = s->groundSet[i]; /* undefined item */
        for (j = 0; j < s->prob->m; ++j)
            if (s->amount[j] + s->prob->cost[j * s->prob->n + z] > s->prob->capacity[j])
                break;
        if (j < s->prob->m) /* item is invalid, exclude it */
            swap_i(s->groundSet, s->igroundSet, i, --oldlim);
    }
    s->n_v = oldlim - newlim;
    /* solution not evaluated yet */
    s->evalv = s->evalLB = 0;
    s->numEnumLim = st_size(s);
    s->sampleWORLim = s->sampleEnumLim = nh_size(s);
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
    /* solution s is feasible */
    int i, j;
    double obj = 0.0;
    if (s->evalv) /* solution s is evaluated */
        *objv = s->objv;
    else { /* solution s is not evaluated */
        for (i = 0; i < s->n_p; ++i)
            for (j = 0; j < s->n_p; ++j)
                obj += s->prob->profit[s->groundSet[i] * s->prob->n + s->groundSet[j]];
        *objv = s->objv = s->prob->sum - obj;
        s->evalv = 1;
    }
    return objv;
}

/*
 * Lower bound evaluation
 * Status: INTERIM
 * Notes:
 *   Implement a tighter upper bound for quadratic multidimensional knapsack
 *   solution
 */
double *getObjectiveLB(double *objLB, struct solution *s)
{
    if (s->evalLB) /* solution s is evaluated */
        *objLB = s->objLB;
    else { /* solution s is not evaluated */
        *objLB = s->objLB = s->prob->sum - ub(s);
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
 * Return true if the items do not exceed the capacity of the multidimensional
 * knapsack or false if they exceed its capacity
 * Status: FINAL
 */
int isFeasible(struct solution *s)
{
    return 1;
}

/*
 * Return true if no items can be added to the multidimensional knapsack or
 * false if there are items that can be added to it
 * Status: FINAL
 */
int isComplete(struct solution *s)
{
    return !s->n_v;
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
    switch (nh) {
    case ADD:
        return s->n_v;
    case REMOVE:
        return s->n_p;
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
        s->sampleEnumLim.add = s->n_v;
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
        s->sampleWORLim.add = s->n_v;
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
    // Coverage: - No valid moves available (s->n_v is zero)
    //               - Complete solutions
    if (!s->n_v)
        return NULL;

    // Sample a random valid component
    v->data = s->groundSet[s->n_p + randint(s->n_v - 1)];

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 4 * sizeof(int));
    return v;
}

static struct move *randomRemoveMove(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available (s->n_p is zero)
    //               - Solution is empty
    if (!s->n_p)
        return NULL;

    // Sample a random valid component
    v->data = s->groundSet[randint(s->n_p - 1)];

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
    //               - All valid components are used
    if (!s->sampleWORLim.add)
        return NULL;

    // Sample a random valid component:
    //           - Available used components cannot be resampled
    // Decrease number of unused moves available
    v->data = s->groundSet[s->n_p + randint(--s->sampleWORLim.add)];

    // Exclude component
    swap_i(s->groundSet, s->igroundSet, s->igroundSet[v->data], s->n_p + s->sampleWORLim.add);

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 4 * sizeof(int));
    return v;
}

static struct move *randomRemoveMoveWOR(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available (s->sampleWORLim.remove is zero)
    //               - Solution is empty
    //               - All valid moves are used
    if (!s->sampleWORLim.remove)
        return NULL;

    // Sample a random valid component:
    //           - Available used components cannot be resampled
    // Decrease number of unused moves available
    v->data = s->groundSet[randint(--s->sampleWORLim.remove)];

    // Exclude component
    swap_i(s->groundSet, s->igroundSet, s->igroundSet[v->data], s->sampleWORLim.remove);

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
    //               - All valid moves are used
    if (!s->sampleEnumLim.add)
        return NULL;

    // Sample a valid component:
    //           - Available used components cannot be resampled
    // Decrease number of unused moves available
    v->data = s->groundSet[s->n_p + s->n_v - s->sampleEnumLim.add--];

    
    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 4 * sizeof(int));
    return v;
}

static struct move *enumRemoveMove(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available (s->sampleEnumLim.remove is zero)
    //               - Solution is empty
    //               - All valid moves are used
    if (!s->sampleEnumLim.remove)
        return NULL;

    // Sample a valid component:
    //           - Available used components cannot be resampled
    // Decrease number of unused moves available
    v->data = s->groundSet[s->n_p - s->sampleEnumLim.remove--];

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
 * Copy the contents of the second argument to the first argument
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
    dest->n_v = src->n_v;
    memcpy(dest->amount, src->amount, src->prob->m * sizeof(double));
    dest->numEnumLim = src->numEnumLim;
    dest->sampleEnumLim = src->sampleEnumLim;
    dest->sampleHeuristicLim = src->sampleHeuristicLim;
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

    // Coverage: - Component is not valid
    //           - Solution is complete
#if __DEBUG__
    if (k < s->n_p || k > s->n_p + s->n_v - 1) {
        fprintf(stderr, "Component %d is not valid or solution is complete.\nCould not add component.\n\n", s->groundSet[k]);
        return NULL;
    }
#endif

    int i, j;

    // Add *undefined* component
    swap_i(s->groundSet, s->igroundSet, s->n_p, k);
    s->n_p++;
    s->n_v--;

    // Update amount occupied by the items in the knapsack
    for (i = 0; i < s->prob->m; ++i)
        s->amount[i] += s->prob->cost[i * s->prob->n + v->data];

    // Update invalid components
    for (i = s->n_p + s->n_v - 1; i > s->n_p - 1; --i)
        for (j = 0; j < s->prob->m; ++j)
            if (s->amount[j] + s->prob->cost[j * s->prob->n + s->groundSet[i]] > s->prob->capacity[j]) {
                swap_i(s->groundSet, s->igroundSet, i, s->n_p + --s->n_v);
                break;
            }

    s->evalv = 0;

    // Update lower bound
    if (s->evalLB && v->evalLBi[0]) {
        s->objLB += v->objLBi;
    }
    else {
        s->evalLB = 0;
    }

    // Update number of components left to enumerate
    s->numEnumLim = st_size(s);

    // Update number of unused moves
    s->sampleWORLim = (s->sampleEnumLim = nh_size(s));
    // s->sampleHeuristicLim = ; // FIX ME

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

    // Update invalid components
    if (k < s->n_p + s->n_v)
        s->n_v--;

    swap_i(s->groundSet, s->igroundSet, k, s->n_p + s->n_v);

    // Forbid *undefined* component
    swap_i(s->groundSet, s->igroundSet, s->n_p + s->n_v, lim);
    s->n_f++;

    s->evalv = 0;

    // Update lower bound
    if (s->evalLB && v->evalLBi[2]) {
        s->objLB += v->objLBi;
    }
    else {
        s->evalLB = 0;
    }

    // Update number of components left to enumerate
    s->numEnumLim = st_size(s);

    // Update number of unused moves
    s->sampleWORLim = (s->sampleEnumLim = nh_size(s));
    // s->sampleHeuristicLim = ; // FIX ME

    return s;
}

static struct solution *applyRemoveMove(struct solution *s, const struct move *v)
{
    int k = s->igroundSet[v->data];
    int lim = s->n_p - 1;

    // Coverage: - Component is not present
    //           - Solution is empty
#if __DEBUG__
    if (k > lim || s->n_p < 1) {
        fprintf(stderr, "Component %d is not present or solution is empty.\nCould not remove component.\n\n", s->groundSet[k]);
        return NULL;
    }
#endif

    int i, j;

    // Remove *present* component
    swap_i(s->groundSet, s->igroundSet, k, lim);
    s->n_p--;
    s->n_v++;

    // Update amount occupied by the items in the knapsack
    for (i = 0; i < s->prob->m; ++i)
        s->amount[i] -= s->prob->cost[i * s->prob->n + v->data];

    // Update valid components
    for (i = s->n_p + s->n_v; i < s->n_g - s->n_f; ++i) {
        for (j = 0; j < s->prob->m; ++j)
            if (s->amount[j] + s->prob->cost[j * s->prob->n + s->groundSet[i]] > s->prob->capacity[j])
                break;
        if (j == s->prob->m)
            swap_i(s->groundSet, s->igroundSet, i, s->n_p + s->n_v++);
    }

    s->evalv = 0;

    // Update lower bound
    if (s->evalLB && v->evalLBi[1]) {
        s->objLB += v->objLBi;
    }
    else {
        s->evalLB = 0;
    }

    // Update number of components left to enumerate
    s->numEnumLim = st_size(s);

    // Update number of unused moves
    s->sampleWORLim = (s->sampleEnumLim = nh_size(s));
    // s->sampleHeuristicLim = ; // FIX ME

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

    int i;

    // Permit *forbidden* component
    swap_i(s->groundSet, s->igroundSet, lim, k);
    s->n_f--;

    // Update valid components
    for (i = 0; i < s->prob->m; ++i)
        if (s->amount[i] + s->prob->cost[i * s->prob->n + v->data] > s->prob->capacity[i])
            break;
    if (i == s->prob->m)
        swap_i(s->groundSet, s->igroundSet, s->n_p + s->n_v++, lim);

    s->evalv = 0;

    // Update lower bound
    if (s->evalLB && v->evalLBi[3]) {
        s->objLB += v->objLBi;
    }
    else {
        s->evalLB = 0;
    }

    // Update number of components left to enumerate
    s->numEnumLim = st_size(s);

    // Update number of unused moves
    s->sampleWORLim = (s->sampleEnumLim = nh_size(s));
    // s->sampleHeuristicLim = ; // FIX ME

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
 *   Implement a tighter upper bound for quadratic multidimensional knapsack
 *   solution
 */
double *getObjectiveLBIncrement(double *obji, struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    int i;
    switch (nh) {
    case ADD:
#if __DEBUG__
        if (s->igroundSet[v->data] < s->n_p || s->igroundSet[v->data] >= s->n_p + s->n_v) { /* item is not valid */
            fprintf(stderr, "Invalid move passed to getObjectiveLBIncrement(). It does not belong to Add neighbourhood.\n");
            return NULL;
        }
#endif
        i = 0;
        break;
    case REMOVE:
#if __DEBUG__
        if (s->igroundSet[v->data] >= s->n_p) { /* item is not present */
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
