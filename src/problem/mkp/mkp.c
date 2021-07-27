#include "mkp.h"
#include "sort_r.h"
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
    double *profit;   /* Profit array */
    double sum;       /* Sum of profit of items */
    double *mu;       /* Multiplier array used for the surrogate constraint */
    double *cost;     /* Cost array */
    double *scost;    /* Cost array regarding the surrogate constraint */
    double *density;  /* Density array */
    int *ordered;     /* Components sorted in decreasing order based on density regarding the surrogate constraint */
    int *iordered;    /* Inverse permutation of the sorted array of components */
    double *capacity; /* Capacity */
    double scapacity; /* Capacity regarding the surrogate constraint */
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
    double samount;                        /* Amount occupied by the items regarding the surrogate constraint */
    int n_v;                               /* Number of valid components */
    struct sample_size sampleWORLim;       /* Number of unused moves for random move without replacement */
    struct sample_size sampleEnumLim;      /* Number of unused moves for enumerate move */
    struct sample_size sampleHeuristicLim; /* Number of unused moves for heuristic move */
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
 * Compare the values of the ath and bth elements of an array.
 * This function is called repeatedly by sort_r().
 */
static int compar_r(const void *_a, const void *_b, void *_arg)
{
    const int a = *(const int *)_a;
    const int b = *(const int *)_b;
    const double *arg = (const double *)_arg;
    return *(arg + a) < *(arg + b);
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
 * Compute the number of heuristically sorted neighbours of a solution that are
 * in each neighbourhood
 */
static struct sample_size br_size(const struct solution *s)
{
    struct sample_size size;
    size.add = size.remove = s->prob->n;
    return size;
}

/*
 * Compute upper bound using a greedy approach that consists of adding the items
 * with highest density and a fraction of the last item that exceeds the
 * remanining capacity of the knapsack
 */
static double ub(struct solution *s)
{
    int lim = s->n_p + s->n_v, z, i;
    double obj = 0.0, cap;
    /* compute profit of solution s */
    if (s->evalv) /* solution s is evaluated */
        obj += s->prob->sum - s->objv;
    else { /* solution s is not evaluated */
        for (i = 0; i < s->n_p; ++i)
            obj += s->prob->profit[s->groundSet[i]];
        s->objv = s->prob->sum - obj;
        s->evalv = 1;
    }
    /* compute the remaining upper bound */
    cap = s->prob->scapacity - s->samount; /* remaining capacity */
    for (i = 0; i < s->prob->n; ++i) {
        z = s->prob->ordered[i]; /* item */
        if (s->igroundSet[z] >= s->n_p && s->igroundSet[z] < lim) { /* item is valid */
            /* add item */
            obj += cap >= s->prob->scost[z] ? s->prob->profit[z] : s->prob->profit[z] * cap / s->prob->scost[z];
            cap -= s->prob->scost[z];
        }
        if (cap <= 0.0) /* no remaining capacity */
            break;
    }
    return obj;
}

/*
 * Compute upper bound after adding an item to a solution
 */
static double ub_add(const struct move *v, struct solution *s)
{
    int lim = s->n_p + s->n_v, z, i, j;
    double obj = 0.0, cap;
    /* compute profit of solution s */
    if (s->evalv) /* solution s is evaluated */
        obj += s->prob->sum - s->objv;
    else { /* solution s is not evaluated */
        for (i = 0; i < s->n_p; ++i)
            obj += s->prob->profit[s->groundSet[i]];
        s->objv = s->prob->sum - obj;
        s->evalv = 1;
    }
    /* add profit of item in move v */
    obj += s->prob->profit[v->data];
    /* compute the remaining upper bound */
    cap = s->prob->scapacity - s->prob->scost[v->data] - s->samount; /* remaining capacity */
    for (i = 0; i < s->prob->n; ++i) {
        z = s->prob->ordered[i]; /* item */
        if (s->igroundSet[z] >= s->n_p && s->igroundSet[z] < lim && z != v->data) { /* item is valid and not in move v */
            for (j = 0; j < s->prob->m; ++j)
                if (s->amount[j] + s->prob->cost[j * s->prob->n + v->data] + s->prob->cost[j * s->prob->n + z] > s->prob->capacity[j])
                    break;
            if (j == s->prob->m) { /* item is valid after adding item in move v */
                /* add item */
                obj += cap >= s->prob->scost[z] ? s->prob->profit[z] : s->prob->profit[z] * cap / s->prob->scost[z];
                cap -= s->prob->scost[z];
            }
        }
        if (cap <= 0.0) /* no remaining capacity */
            break;
    }
    return obj;
}

/*
 * Compute upper bound after removing an item from a solution
 */
static double ub_remove(const struct move *v, struct solution *s)
{
    int lim1 = s->n_p + s->n_v, lim2 = s->n_g - s->n_f, z, i, j;
    double obj = 0.0, cap;
    /* compute profit of solution s */
    if (s->evalv) /* solution s is evaluated */
        obj += s->prob->sum - s->objv;
    else { /* solution s is not evaluated */
        for (i = 0; i < s->n_p; ++i)
            obj += s->prob->profit[s->groundSet[i]];
        s->objv = s->prob->sum - obj;
        s->evalv = 1;
    }
    /* remove profit of item in move v */
    obj -= s->prob->profit[v->data];
    /* compute the remaining upper bound */
    cap = s->prob->scapacity + s->prob->scost[v->data] - s->samount; /* remaining capacity */
    for (i = 0; i < s->prob->n; ++i) {
        z = s->prob->ordered[i]; /* item */
        if ((s->igroundSet[z] >= s->n_p && s->igroundSet[z] < lim1) || z == v->data) { /* item is valid or in move v */
            /* add item */
            obj += cap >= s->prob->scost[z] ? s->prob->profit[z] : s->prob->profit[z] * cap / s->prob->scost[z];
            cap -= s->prob->scost[z];
        }
        else if (s->igroundSet[z] >= lim1 && s->igroundSet[z] < lim2) { /* item is invalid */
            for (j = 0; j < s->prob->m; ++j)
                if (s->amount[j] - s->prob->cost[j * s->prob->n + v->data] + s->prob->cost[j * s->prob->n + z] > s->prob->capacity[j])
                    break;
            if (j == s->prob->m) { /* item is valid after removing item in move v */
                /* add item */
                obj += cap >= s->prob->scost[z] ? s->prob->profit[z] : s->prob->profit[z] * cap / s->prob->scost[z];
                cap -= s->prob->scost[z];
            }
        }
        if (cap <= 0.0) /* no remaining capacity */
            break;
    }
    return obj;
}

/*
 * Compute upper bound after forbidding an item in a solution
 */
static double ub_forbid(const struct move *v, struct solution *s)
{
    int lim = s->n_p + s->n_v, z, i;
    double obj = 0.0, cap;
    /* compute profit of solution s */
    if (s->evalv) /* solution s is evaluated */
        obj += s->prob->sum - s->objv;
    else { /* solution s is not evaluated */
        for (i = 0; i < s->n_p; ++i)
            obj += s->prob->profit[s->groundSet[i]];
        s->objv = s->prob->sum - obj;
        s->evalv = 1;
    }
    /* compute the remaining upper bound */
    cap = s->prob->scapacity - s->samount; /* remaining capacity */
    for (i = 0; i < s->prob->n; ++i) {
        z = s->prob->ordered[i]; /* item */
        if (s->igroundSet[z] >= s->n_p && s->igroundSet[z] < lim && z != v->data) { /* item is valid and not in move v */
            /* add item */
            obj += cap >= s->prob->scost[z] ? s->prob->profit[z] : s->prob->profit[z] * cap / s->prob->scost[z];
            cap -= s->prob->scost[z];
        }
        if (cap <= 0.0) /* no remaining capacity */
            break;
    }
    return obj;
}

/*
 * Compute upper bound after permitting an item in a solution
 */
static double ub_permit(const struct move *v, struct solution *s)
{
    int lim = s->n_p + s->n_v, z, i, j;
    double obj = 0.0, cap;
    /* compute profit of solution s */
    if (s->evalv) /* solution s is evaluated */
        obj += s->prob->sum - s->objv;
    else { /* solution s is not evaluated */
        for (i = 0; i < s->n_p; ++i)
            obj += s->prob->profit[s->groundSet[i]];
        s->objv = s->prob->sum - obj;
        s->evalv = 1;
    }
    /* compute the remaining upper bound */
    cap = s->prob->scapacity - s->samount; /* remaining capacity */
    for (i = 0; i < s->prob->n; ++i) {
        z = s->prob->ordered[i]; /* item */
        if (s->igroundSet[z] >= s->n_p && s->igroundSet[z] < lim) { /* item is valid */
            /* add item */
            obj += cap >= s->prob->scost[z] ? s->prob->profit[z] : s->prob->profit[z] * cap / s->prob->scost[z];
            cap -= s->prob->scost[z];
        }
        else if (z == v->data) { /* item is in move v */
            for (j = 0; j < s->prob->m; ++j)
                if (s->amount[j] + s->prob->cost[j * s->prob->n + z] > s->prob->capacity[j])
                    break;
            if (j == s->prob->m) { /* item is valid after removing item in move v */
                /* add item */
                obj += cap >= s->prob->scost[z] ? s->prob->profit[z] : s->prob->profit[z] * cap / s->prob->scost[z];
                cap -= s->prob->scost[z];
            }
        }
        if (cap <= 0.0) /* no remaining capacity */
            break;
    }
    return obj;
}

/*************************/
/* Problem instantiation */
/*************************/

/*
 * Multidimensional knapsack instantiation
 * Status: TENTATIVE
 * Notes:
 *   Should compute the surrogate dual using linear programming
 *   Needs further error checking
 */
struct problem *newProblem(const char *filename)
{
    int n, m, i, j;
    double aux = 0.0;
    FILE *inFile;
    struct problem *p = NULL;
    inFile = fopen(filename, "r");
    if (inFile) {
        fscanf(inFile, "%d %d %lf", &n, &m, &aux);
        if (n > 0 && m > 0) {
            p = malloc(sizeof(struct problem));
            p->n = n;
            p->m = m;
            p->profit = malloc(sizeof(double) * n);
            p->mu = malloc(sizeof(double) * m);
            p->cost = malloc(sizeof(double) * m * n);
            p->scost = malloc(sizeof(double) * n);
            p->density = malloc(sizeof(double) * n);
            p->capacity = malloc(sizeof(double) * m);
            for (i = 0, p->sum = 0.0; i < n; ++i) {
                fscanf(inFile, "%lf", &p->profit[i]);
                p->sum += p->profit[i];
            }
            for (i = 0; i < m; p->mu[i++] = 1.0)
                for (j = 0; j < n; ++j)
                    fscanf(inFile, "%lf", &p->cost[i * n + j]);
            for (i = 0; i < n; ++i) {
                p->scost[i] = 0.0;
                for (j = 0; j < m; ++j)
                    p->scost[i] += p->mu[j] * p->cost[j * n + i];
                if (p->scost[i] > 0.0 && p->profit[i] / p->scost[i] > aux)
                    aux = p->profit[i] / p->scost[i];
            }
            for (i = 0; i < n; ++i)
                p->density[i] = p->scost[i] > 0.0 ? p->profit[i] / p->scost[i] : aux + p->profit[i];
            for (i = 0, p->scapacity = 0.0; i < m; ++i) {
                fscanf(inFile, "%lf", &p->capacity[i]);
                p->scapacity += p->mu[i] * p->capacity[i];
            }
            p->ordered = malloc(sizeof(int) * n);
            for (i = 0; i < n; ++i)
                p->ordered[i] = i;
            sort_r(p->ordered, n, sizeof(int), compar_r, p->density); /* portable qsort_r implementation by Isaac Turner in https://github.com/noporpoise/sort_r */
            p->iordered = malloc(sizeof(int) * n);
            for (i = 0; i < n; ++i)
                p->iordered[p->ordered[i]] = i;
        }
        else
            fprintf(stderr, "Invalid multidimensional knapsack instance %s.\n", filename);
        fclose(inFile);
    }
    else
        fprintf(stderr, "Cannot open file %s.\n", filename);
    return p;
}

/**********************/
/* Problem inspection */
/**********************/

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
    free(p->mu);
    free(p->cost);
    free(p->scost);
    free(p->density);
    free(p->capacity);
    free(p->ordered);
    free(p->iordered);
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
    printf("Multidimensional Knapsack (%d) problem\n  Profit:  ", p->n);
    for (i = 0; i < p->n; ++i)
        printf(" %*.1lf", 9, p->profit[i]);
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
    printf("Multidimensional Knapsack (%d) solution\n  x:", s->prob->n);
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
    printf("Multidimensional Knapsack (%d) move: %d\n\n", v->prob->n, v->data);
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
    s->n_v = s->n_g;
    memset(s->amount, 0, s->prob->m * sizeof(double));
    s->samount = 0.0;
    /* solution not evaluated yet */
    s->evalv = s->evalLB = 0;
    s->numEnumLim = st_size(s);
    s->sampleWORLim = s->sampleEnumLim = nh_size(s);
    s->sampleHeuristicLim = br_size(s);
    return s;
}

/*
 * Heuristically constructs a multidimensional knapsack solution using the
 * algorithm proposed in Akçay, Y., Li, H., & Xu, S. H. (2007). Greedy Algorithm
 * for the General Multidimensional Knapsack Problem. Annals of Operations
 * Research, 150(1), 17-29.
 * Status: CHECK
 */
struct solution *heuristicSolution(struct solution *s)
{
    if (!s->n_v) /* solution s is complete, return it */
        return s;
    int z, item, load, minload, i, j;
    double value, maxvalue;
    /* Primal Effective Capacity Heuristic (PECH) */
    do {
        for (i = s->n_p + s->n_v - 1, item = -1, maxvalue = 0.0; i >= s->n_p; --i) {
            z = s->groundSet[i]; /* undefined item */
            /* compute maximum number of copies of item that fit into the knapsack */
            for (j = 0, minload = INT_MAX; j < s->prob->m; ++j) {
                load = s->prob->cost[j * s->prob->n + z] > 0 ? (s->prob->capacity[j] - s->amount[j]) / s->prob->cost[j * s->prob->n + z] : INT_MAX;
                if (load < minload)
                    minload = load;
                if (!minload)
                    break;
            }
            if (!minload) /* item is invalid, exclude it */
                swap_i(s->groundSet, s->igroundSet, i, s->n_p + --s->n_v);
            else {
                /* compute maximum reward of item */
                value = s->prob->profit[z] * minload;
                if (value > maxvalue) {
                    item = z; /* item with largest attainable reward */
                    maxvalue = value;
                }
            }
        }
        if (item > -1) {
            /* update solution */
            swap_i(s->groundSet, s->igroundSet, s->n_p++, s->igroundSet[item]);
            s->n_v--;
            for (i = 0; i < s->prob->m; ++i)
                s->amount[i] += s->prob->cost[i * s->prob->n + item];
            s->samount += s->prob->scost[item];
        }
    } while (s->n_v);
    /* solution not evaluated yet */
    s->evalv = s->evalLB = 0;
    s->numEnumLim = st_size(s);
    s->sampleWORLim = s->sampleEnumLim = nh_size(s);
    s->sampleHeuristicLim = br_size(s);
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
    int i;
    double obj = 0.0;
    if (s->evalv) /* solution s is evaluated */
        *objv = s->objv;
    else { /* solution s is not evaluated */
        for (i = 0; i < s->n_p; ++i)
            obj += s->prob->profit[s->groundSet[i]];
        *objv = s->objv = s->prob->sum - obj;
        s->evalv = 1;
    }
    return objv;
}

/*
 * Lower bound evaluation
 * Status: FINAL
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

static struct move *heuristicAddMove(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available (s->n_v is zero)
    //               - Complete solutions
    //               - All moves are present, forbidden or invalid
    if (!s->n_v)
        return NULL;

    int i, k, idx, pres_lim = s->n_p - 1, valid_lim = s->n_p + s->n_v;

    // Sample a valid component:
    //           - Cannot be present, forbidden or invalid
    for (i = 0; i < s->prob->n; ++i) {
        // Get next item k in sequence of components in decreasing order based on density
        k = s->prob->ordered[i];
        // Get position of components in component set
        idx = s->igroundSet[k];
        // Component is valid
        if (idx > pres_lim && idx < valid_lim) {
            v->data = k;
            break;
        }
    }

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 4 * sizeof(int));
    return v;
}

static struct move *heuristicRemoveMove(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available (s->n_p is zero)
    //               - Solution is empty
    //               - All moves are undefined or forbidden
    if (!s->n_p)
        return NULL;

    int i, k, idx, pres_lim = s->n_p;

    // Sample a valid component:
    //           - Cannot be undefined or forbidden
    for (i = s->prob->n - 1; i > -1; --i) {
        // Get previous item k in sequence of components in decreasing order based on density
        k = s->prob->ordered[i];
        // Get position of components in component set
        idx = s->igroundSet[k];
        // Component is present
        if (idx < pres_lim) {
            v->data = k;
            break;
        }
    }

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 4 * sizeof(int));
    return v;
}

struct move *heuristicMove(struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return heuristicAddMove(v, s);

    case REMOVE:
        return heuristicRemoveMove(v, s);

    case FORBID:
        return NULL; // TODO: Use unimplemented() function.

    case PERMIT:
        return NULL; // TODO: Use unimplemented() function.

    default:
        return NULL; // TODO: Use unimplemented() function.
    }
}

static struct move *heuristicAddMoveWOR(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available (s->n_v and s->sampleHeuristicLim.add are zero)
    //               - Complete solutions
    //               - All moves are present, forbidden, invalid or used
    if (!s->n_v || !s->sampleHeuristicLim.add)
        return NULL;

    int i, k, idx, pres_lim = s->n_p - 1, valid_lim = s->n_p + s->n_v;

    // Sample a valid component:
    //           - Cannot be present, forbidden or invalid
    //           - Available used components cannot be resampled
    for (i = s->prob->n - s->sampleHeuristicLim.add; i < s->prob->n; ++i) {
        // Get next item k in sequence of components in decreasing order based on density
        k = s->prob->ordered[i];
        // Decrease number of valid (and unused) moves available
        s->sampleHeuristicLim.add--;
        // Get position of components in component set
        idx = s->igroundSet[k];
        // Component is valid
        if (idx > pres_lim && idx < valid_lim) {
            v->data = k;
            break;
        }
    }

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 4 * sizeof(int));
    return v;
}

static struct move *heuristicRemoveMoveWOR(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available (s->sampleHeuristicLim.remove is zero)
    //               - Solution is empty
    //               - All moves are undefined, forbidden or used
    if (!s->n_p || !s->sampleHeuristicLim.remove)
        return NULL;

    int i, k, idx, pres_lim = s->n_p;

    // Sample a valid component:
    //           - Cannot be undefined or forbidden
    //           - Available used components cannot be resampled
    for (i = s->sampleHeuristicLim.remove - 1; i > -1; --i) {
        // Get previous item k in sequence of components in decreasing order based on density
        k = s->prob->ordered[i];
        // Decrease number of valid (and unused) moves available
        s->sampleHeuristicLim.remove--;
        // Get position of components in component set
        idx = s->igroundSet[k];
        // Component is present
        if (idx < pres_lim) {
            v->data = k;
            break;
        }
    }

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 4 * sizeof(int));
    return v;
}

struct move *heuristicMoveWOR(struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return heuristicAddMoveWOR(v, s);

    case REMOVE:
        return heuristicRemoveMoveWOR(v, s);

    case FORBID:
        return NULL; // TODO: Use unimplemented() function.

    case PERMIT:
        return NULL; // TODO: Use unimplemented() function.

    default:
        return NULL; // TODO: Use unimplemented() function.
    }
}

static struct solution *resetHeuristicAddMoveWOR(struct solution *s)
{
    s->sampleHeuristicLim.add = s->prob->n;
    return s;
}

static struct solution *resetHeuristicRemoveMoveWOR(struct solution *s)
{
    s->sampleHeuristicLim.remove = s->prob->n;
    return s;
}

struct solution *resetHeuristicMoveWOR(struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return resetHeuristicAddMoveWOR(s);

    case REMOVE:
        return resetHeuristicRemoveMoveWOR(s);

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
    dest->samount = src->samount;
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

    // Update amount occupied by the items regarding the surrogate constraint
    s->samount += s->prob->scost[v->data];

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
    s->sampleHeuristicLim = br_size(s);

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
    s->sampleHeuristicLim = br_size(s);

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

    // Update amount occupied by the items regarding the surrogate constraint
    s->samount -= s->prob->scost[v->data];

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
    s->sampleHeuristicLim = br_size(s);

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
    s->sampleHeuristicLim = br_size(s);

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
 * Status: TENTATIVE
 * Notes:
 *   Redefine the concept of 'subneighbourhood'
 */
double *getObjectiveLBIncrement(double *obji, struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    int i;
    double obj1, obj2;
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
        if (s->evalLB) /* solution s is evaluated */
            obj1 = s->objLB;
        else { /* solution s is not evaluated */
            obj1 = s->objLB = s->prob->sum - ub(s);
            s->evalLB = 1;
        }
        switch (nh) {
        case ADD:
            obj2 = s->prob->sum - ub_add(v, s);
            break;
        case REMOVE:
            obj2 = s->prob->sum - ub_remove(v, s);
            break;
        case FORBID:
            obj2 = s->prob->sum - ub_forbid(v, s);
            break;
        case PERMIT:
            obj2 = s->prob->sum - ub_permit(v, s);
            break;
        default:
            obj2 = -DBL_MAX;
            break;
        }
        *obji = v->objLBi = obj2 - obj1;
        v->evalLBi[i] = 1;
    }
    return obji;
}
