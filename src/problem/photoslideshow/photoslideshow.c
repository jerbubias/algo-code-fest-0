#include "photoslideshow.h"
#include <float.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_TAGS 100

struct state_size {
    int present; /* Number of present components */
};

struct sample_size {
    int add;    /* Size of subneighborhood resulting by adding a component to the solution */
    int remove; /* Size of subneighborhood resulting by removing a component from the solution */
};

struct photo {
    int id;   /* Photo identifier */
    char o;   /* Photo orientation */
    int oid;  /* Photo identifier regarding orientation */
    int m;    /* Number of tags in the photo */
    int *tgs; /* Tags in the photo */
};

struct problem {
    int n;                  /* Number of photos in the collection */
    int u;                  /* Number of unique tags in the problem instance */
    int h;                  /* Number of horizontal photos in the collection */
    int v;                  /* Number of vertical photos in the collection */
    int p;                  /* Maximum number of photos in the photo slideshow */
    int s;                  /* Maximum number of slides in the photo slideshow */
    long e;                 /* Number of components in the problem instance */
    struct photo *photos;   /* Photos in the collection */
    struct photo **hphotos; /* Horizontal photos in the collection */
    struct photo **vphotos; /* Vertical photos in the collection */
    double maxUB;           /* Maximum upper bound value for the problem instance */
    double maxUBi;          /* Maximum upper bound value for a photo transition */
};

struct solution {
    // Problem Abstraction
    struct problem *prob;
    int evalv;                    /* Flag indicating if the solution is evaluated */
    double objv;                  /* Objective value */
    int evalLB;                   /* Flag indicating if the lower bound is calculated */
    double objLB;                 /* Objective Lower Bound */
    struct state_size numEnumLim; /* Number of components left to enumerate */
    // Problem-Specific
    int *pss;                         /* Sequence of photos in the slideshow */
    int *ipss;                        /* Inverse permutation of sequence of photos in the slideshow */
    int i_f;                          /* Index of the first photo in the slideshow */
    int n_p;                          /* Number of photos added to the photo slideshow */
    int side;                         /* Flag referring to the side where last photo was placed */
    int *vss;                         /* Vertical photos in the slideshow */
    int *ivss;                        /* Inverse permutation of vertical photos in the slideshow */
    int n_v;                          /* Number of vertical photos added to the slideshow */
    struct sample_size sampleWORLim;  /* Number of unused moves for random move without replacement */
    struct sample_size sampleEnumLim; /* Number of unused moves for enumerate move */
    int *tgs;
};

struct move {
    struct problem *prob;
    long data;      /* Component identifier */
    int photos[2];  /* Photos */
    int side;       /* Flag referring to the side where photo is placed */
    int evalLBi[2]; /* Flag indicating if lower bound increment is evaluated for operations: { 0 - Add, 1 - Remove } */
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
 * Compare the values of a and b.
 * This function is called repeatedly by qsort().
 */
static int compar(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

/*
 * Uniquely encode two distinct natural numbers x and y into a single natural
 * number z
 */
static long pairing(const int x, const int y)
{
    return x > y ? (long)x * (x - 1) / 2 + y : (long)y * (y - 1) / 2 + x;
}

/*
 * Uniquely decode a single natural number z into two distinct natural numbers x
 * and y
 */
static void unpairing(int *x, int *y, const long z)
{
    int w = (-1 + sqrt(1 + 8 * z)) / 2;
    *x = w + 1;
    *y = z - (long)w * *x / 2;
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
    size.present = s->n_p - 1;
    return size;
}

/*
 * Compute the number of neighbours of a solution that are in each
 * subneighbourhood
 */
static struct sample_size nh_size(struct solution *s)
{
    struct sample_size size;
    int a = s->n_p > 1 ? 2 : 1, b = s->n_p < s->prob->p;
    size.add = s->n_v % 2 ? s->prob->v - s->n_v : a * (s->prob->n - s->n_p) * b;
    size.remove = s->n_p > 2 ? a : a - 1;
    return size;
}

/*
 * Return the smallest of a, b and c
 */
static int min(const int a, const int b, const int c)
{
    return a < b ? a < c ? a : c : b < c ? b : c;
}

/*
 * Combines the elements in the sorted arrays tag1 and tag2, into a new array
 * tag with all its elements sorted
 */
static int merge(int *tag, const int *tag1, const int n1, const int *tag2, const int n2)
{
    int result, n = 0, i = 0, j = 0;
    while (i < n1 && j < n2) {
        result = tag1[i] - tag2[j];
        if (result < 0)
            tag[n++] = tag1[i++];
        else if (result > 0)
            tag[n++] = tag2[j++];
        else {
            tag[n++] = tag1[i++];
            j++;
        }
    }
    while (i < n1)
        tag[n++] = tag1[i++];
    while (j < n2)
        tag[n++] = tag2[j++];
    return n;
}

/*
 * Return the number elements that are present in both sorted arrays tag1 and
 * tag2
 */
static int common(const int *tag1, const int n1, const int *tag2, const int n2)
{
    int result, n = 0, i = 0, j = 0;
    while (i < n1 && j < n2) {
        result = tag1[i] - tag2[j];
        if (result < 0)
            i++;
        else if (result > 0)
            j++;
        else {
            n++;
            i++;
            j++;
        }
    }
    return n;
}

/*
 * Compute upper bound
 */
static double ub(struct solution *s)
{
    int p1, p2, *tag1, n1, *tag2, n2, n3, bit, n = s->i_f + s->n_p - (s->n_v % 2 && s->side), i = s->i_f + (s->n_v % 2 && !s->side);
    double obj = 0.0;
    if (s->evalv) /* solution s is evaluated */
        obj += s->prob->maxUB - s->objv;
    else { /* solution s is not evaluated */
        /* get first valid slide */
        p1 = s->pss[i % s->prob->n];
        p2 = s->pss[(i + 1) % s->prob->n];
        if (s->prob->photos[p1].o == 'V') {
            n1 = merge(s->tgs, s->prob->photos[p1].tgs, s->prob->photos[p1].m, s->prob->photos[p2].tgs, s->prob->photos[p2].m);
            tag1 = s->tgs;
            i++;
        }
        else {
            n1 = s->prob->photos[p1].m;
            tag1 = s->prob->photos[p1].tgs;
        }
        /* compute interest factor for each valid transition */
        for (i = i + 1, bit = 1; i < n; ++i, bit = 1 - bit) {
            /* get next valid slide */
            p1 = s->pss[i % s->prob->n];
            p2 = s->pss[(i + 1) % s->prob->n];
            if (s->prob->photos[p1].o == 'V') {
                n2 = merge(s->tgs + MAX_TAGS * bit * 2, s->prob->photos[p1].tgs, s->prob->photos[p1].m, s->prob->photos[p2].tgs, s->prob->photos[p2].m);
                tag2 = s->tgs + MAX_TAGS * bit * 2;
                i++;
            }
            else {
                n2 = s->prob->photos[p1].m;
                tag2 = s->prob->photos[p1].tgs;
            }
            n3 = common(tag1, n1, tag2, n2); /* number of common tags */
            obj += min(n1 - n3, n3, n2 - n3);
            /* update previous slide */
            n1 = n2;
            tag1 = tag2;
        }
    }
    /* compute upper bound for invalid transition */
    if (s->n_v % 2) { /* solution s is unfeasible */
        // TODO: Compute upper bound for invalid transition
    }
    else { /* solution s is feasible */
        /* update objective value */
        s->objv = s->prob->maxUB - obj;
        s->evalv = 1;
    }
    /* compute upper bound for remaining transitions */
    // TODO: Compute upper bound for a endpoint transition
    // TODO: Compute upper bound for other transitions
    int t = s->n_p - (s->n_v + 1) / 2 - 1;
    obj += s->prob->maxUB - (t >= 0 ? t : 0) * s->prob->maxUBi; // FIX ME: Compute a tighter upper bound for invalid and remaining transitions
    return obj;
}

/*
 * Compute lower bound increment after adding a photo to a solution
 */
static double lbi_add(const struct move *v, const struct solution *s)
{
    int p1, p2, *tag1, n1, *tag2, n2, n3, i;
    double obj = 0.0;
    /* compute interest factor for transition */
    if (s->prob->photos[v->photos[1]].o == 'H') { /* photo is horizontal */
        /* get first slide */
        n1 = s->prob->photos[v->photos[1]].m;
        tag1 = s->prob->photos[v->photos[1]].tgs;
        /* get second slide */
        i = s->i_f + (s->n_p - 1 - (s->prob->photos[s->pss[(s->i_f + s->n_p - 1) % s->prob->n]].o == 'V')) * v->side;
        p1 = s->pss[i % s->prob->n];
        p2 = s->pss[(i + 1) % s->prob->n];
        if (s->prob->photos[p1].o == 'V') {
            n2 = merge(s->tgs, s->prob->photos[p1].tgs, s->prob->photos[p1].m, s->prob->photos[p2].tgs, s->prob->photos[p2].m);
            tag2 = s->tgs;
        }
        else {
            n2 = s->prob->photos[p1].m;
            tag2 = s->prob->photos[p1].tgs;
        }
        n3 = common(tag1, n1, tag2, n2); /* number of common tags */
        obj -= min(n1 - n3, n3, n2 - n3) - s->prob->maxUBi;
    }
    else if (s->n_p > 1 && s->n_v % 2) { /* photo is vertical */
        /* get first slide */
        p1 = s->pss[(s->i_f + (s->n_p - 1) * v->side) % s->prob->n];
        n1 = merge(s->tgs, s->prob->photos[p1].tgs, s->prob->photos[p1].m, s->prob->photos[v->photos[1]].tgs, s->prob->photos[v->photos[1]].m);
        tag1 = s->tgs;
        /* get second slide */
        i = s->i_f + 1 + (s->n_p - 3 - (s->prob->photos[s->pss[(s->i_f + s->n_p - 2) % s->prob->n]].o == 'V')) * v->side;
        p1 = s->pss[i % s->prob->n];
        p2 = s->pss[(i + 1) % s->prob->n];
        if (s->prob->photos[p1].o == 'V') {
            n2 = merge(s->tgs + MAX_TAGS * 2, s->prob->photos[p1].tgs, s->prob->photos[p1].m, s->prob->photos[p2].tgs, s->prob->photos[p2].m);
            tag2 = s->tgs + MAX_TAGS * 2;
        }
        else {
            n2 = s->prob->photos[p1].m;
            tag2 = s->prob->photos[p1].tgs;
        }
        n3 = common(tag1, n1, tag2, n2); /* number of common tags */
        obj -= min(n1 - n3, n3, n2 - n3) - s->prob->maxUBi;
    }
    /* compute upper bound increment for invalid transition */
    if (s->prob->photos[v->photos[1]].o == 'V' && !(s->n_v % 2)) {
        // TODO: Compute upper bound increment for invalid transition
    }
    else if (s->n_v % 2) {
    }
    /* compute upper bound increment for remaining transitions */
    // TODO: Compute upper bound increment for a endpoint transition
    // TODO: Compute upper bound increment for other transitions
    return obj;
}

/*
 * Compute lower bound increment after removing a photo from a solution
 */
static double lbi_remove(const struct move *v, const struct solution *s)
{
    int p1, p2, *tag1, n1, *tag2, n2, n3, i;
    double obj = 0.0;
    /* compute interest factor for transition */
    if (s->prob->photos[v->photos[1]].o == 'H') { /* photo is horizontal */
        /* get first slide */
        n1 = s->prob->photos[v->photos[1]].m;
        tag1 = s->prob->photos[v->photos[1]].tgs;
        /* get second slide */
        i = s->i_f + 1 + (s->n_p - 3 - (s->prob->photos[s->pss[(s->i_f + s->n_p - 2) % s->prob->n]].o == 'V')) * v->side;
        p1 = s->pss[i % s->prob->n];
        p2 = s->pss[(i + 1) % s->prob->n];
        if (s->prob->photos[p1].o == 'V') {
            n2 = merge(s->tgs, s->prob->photos[p1].tgs, s->prob->photos[p1].m, s->prob->photos[p2].tgs, s->prob->photos[p2].m);
            tag2 = s->tgs;
        }
        else {
            n2 = s->prob->photos[p1].m;
            tag2 = s->prob->photos[p1].tgs;
        }
        n3 = common(tag1, n1, tag2, n2); /* number of common tags */
        obj += min(n1 - n3, n3, n2 - n3) - s->prob->maxUBi;
    }
    else if (s->n_p > 2 && !(s->n_v % 2)) { /* photo is vertical */
        /* get first slide */
        p2 = s->pss[(s->i_f + 1 + (s->n_p - 3) * v->side) % s->prob->n];
        n1 = merge(s->tgs, s->prob->photos[v->photos[1]].tgs, s->prob->photos[v->photos[1]].m, s->prob->photos[p2].tgs, s->prob->photos[p2].m);
        tag1 = s->tgs;
        /* get second slide */
        i = s->i_f + 2 + (s->n_p - 5 - (s->prob->photos[s->pss[(s->i_f + s->n_p - 3) % s->prob->n]].o == 'V')) * v->side;
        p1 = s->pss[i % s->prob->n];
        p2 = s->pss[(i + 1) % s->prob->n];
        if (s->prob->photos[p1].o == 'V') {
            n2 = merge(s->tgs + MAX_TAGS * 2, s->prob->photos[p1].tgs, s->prob->photos[p1].m, s->prob->photos[p2].tgs, s->prob->photos[p2].m);
            tag2 = s->tgs + MAX_TAGS * 2;
        }
        else {
            n2 = s->prob->photos[p1].m;
            tag2 = s->prob->photos[p1].tgs;
        }
        n3 = common(tag1, n1, tag2, n2); /* number of common tags */
        obj += min(n1 - n3, n3, n2 - n3) - s->prob->maxUBi;
    }
    /* compute upper bound increment for invalid transition */
    if (s->prob->photos[v->photos[1]].o == 'V' && s->n_v % 2) {
        // TODO: Compute upper bound increment for invalid transition
    }
    else if (!(s->n_v % 2)) {
    }
    /* compute upper bound increment for remaining transitions */
    // TODO: Compute upper bound increment for a endpoint transition
    // TODO: Compute upper bound increment for other transitions
    return obj;
}

/*************************/
/* Problem instantiation */
/*************************/

/*
 * Photo slideshow instantiation
 * Status: TENTATIVE
 * Notes:
 *   Needs further error checking
 */
struct problem *newProblem(const char *filename)
{
    int n, u, h = 0, v = 0, i, j;
    FILE *infile;
    struct problem *p = NULL;
    infile = fopen(filename, "r");
    if (infile) {
        fscanf(infile, "%d %d\n", &n, &u);
        if (n > 2) {
            p = malloc(sizeof(struct problem));
            p->n = n;
            p->u = u;
            p->photos = malloc(sizeof(struct photo) * n);
            for (i = 0; i < n; ++i) {
                fscanf(infile, "%c %d", &p->photos[i].o, &p->photos[i].m);
                p->photos[i].tgs = malloc(sizeof(int) * p->photos[i].m);
                for (j = 0; j < p->photos[i].m; ++j)
                    fscanf(infile, " %d", &p->photos[i].tgs[j]);
                fscanf(infile, "\n");
                qsort(p->photos[i].tgs, p->photos[i].m, sizeof(int), compar);
                if (p->photos[i].o == 'H')
                    p->photos[i].oid = h++;
                else if (p->photos[i].o == 'V')
                    p->photos[i].oid = v++;
                p->photos[i].id = i;
            }
            p->h = h;
            p->v = v;
            p->s = h + v / 2;
            p->p = n - v % 2;
            p->e = (long)n * (n - 1) / 2;
            p->hphotos = malloc(sizeof(struct photo *) * h);
            p->vphotos = malloc(sizeof(struct photo *) * v);
            for (i = 0, h = 0, v = 0; i < n; ++i)
                if (p->photos[i].o == 'H')
                    p->hphotos[h++] = p->photos + i;
                else if (p->photos[i].o == 'V')
                    p->vphotos[v++] = p->photos + i;
            p->maxUBi = (double)MAX_TAGS; // FIX ME: Compute maximum upper bound value for a photo transition
            p->maxUB = (p->s - 1.0) * p->maxUBi;
        }
        else
            fprintf(stderr, "Invalid photo slideshow instance %s.\n", filename);
        fclose(infile);
    }
    else
        fprintf(stderr, "Cannot open file %s.\n", filename);
    return p;
}

/**********************/
/* Problem inspection */
/**********************/

/*
 * Return the number of possible trasitions between each pair of photos
 * Status: FINAL
 */
long getNumComponents(const struct problem *p)
{
    return p->e;
}

/*
 * Return the maximum number of transitions in a photo slideshow
 * Status: FINAL
 */
long getMaxSolutionSize(const struct problem *p)
{
    return p->n - 1 - (p->v % 2);
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
        return (p->n - 2) * 2;
    case REMOVE:
        return 2;
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
    s->pss = malloc(sizeof(int) * p->n);
    s->ipss = malloc(sizeof(int) * p->n);
    s->vss = malloc(sizeof(int) * p->v);
    s->ivss = malloc(sizeof(int) * p->v);
    s->tgs = malloc(sizeof(int) * MAX_TAGS * 4);
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
    int i;
    for (i = 0; i < p->n; ++i)
        free(p->photos[i].tgs);
    free(p->vphotos);
    free(p->hphotos);
    free(p->photos);
    free(p);
}

/*
 * Free the memory used by a solution
 * Status: CHECK
 */
void freeSolution(struct solution *s)
{
    free(s->pss);
    free(s->ipss);
    free(s->vss);
    free(s->ivss);
    free(s->tgs);
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
    printf("Photo SlideShow (%d) problem\n  The collection has %d photos", p->n, p->n);
    for (i = 0; i < p->n; ++i) {
        if (p->photos[i].o == 'H')
            printf("\n  Photo %d is horizontal and has tags [", i);
        else if (p->photos[i].o == 'V')
            printf("\n  Photo %d is vertical and has tags [", i);
        for (j = 0; j < p->photos[i].m; ++j)
            printf(j < p->photos[i].m - 1 ? "%d, " : "%d]", p->photos[i].tgs[j]);
    }
    printf("\n\n");
}

/*
 * Print user-formatted representation of solution
 * Status: CHECK
 */
void printSolution(const struct solution *s)
{
    int p, cnt, bit = !s->side && s->n_v % 2 ? -1 : 0, n = s->i_f + s->n_p, i;
    printf("Photo SlideShow (%d) solution\n  The slideshow has %d slides", s->prob->n, s->n_p - s->n_v / 2);
    for (i = s->i_f, cnt = 0; i < n; ++i) {
        p = s->pss[i % s->prob->n];
        if (s->prob->photos[p].o == 'H')
            printf("\n  Slide %d contains photo %d", cnt++, p);
        else if (bit > 0)
            printf(" %d", p), bit = 0;
        else if (s->prob->photos[p].o == 'V')
            printf("\n  Slide %d contains photos %d and", cnt++, p), bit++;
    }
    if (s->n_p == s->prob->p)
        printf("\n  solution is feasible and complete");
    else if (s->n_v % 2)
        printf("\n  solution is unfeasible and partial");
    else
        printf("\n  solution is feasible and partial");
    if (s->evalLB)
        printf("\n  upper bound: %*.1lf", 14, s->prob->maxUB - s->objLB);
    if (s->evalv)
        printf("\n  objective value: %*.1lf", 10, s->prob->maxUB - s->objv);
    printf("\n\n");
}

/*
 * Print user-formatted representation of move
 * Status: CHECK
 */
void printMove(const struct move *v)
{
    printf("Photo SlideShow (%d) move: %lu (%d, %d)\n\n", v->prob->n, v->data, v->photos[0], v->photos[1]);
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
    int k, i;
    for (i = 0; i < s->prob->n; ++i)
        s->pss[i] = s->ipss[i] = i;
    for (i = 0; i < s->prob->v; ++i)
        s->vss[i] = s->ivss[i] = i;
    s->n_p = s->i_f = s->n_v = 0;
    s->side = 1;
    /* add photo uniformly at random, prioritizing horizontal photos */
    if (s->prob->h) {
        k = randint(s->prob->h - 1);
        swap_i(s->pss, s->ipss, s->ipss[s->prob->hphotos[k]->id], s->n_p++);
    }
    else {
        k = randint(s->prob->v - 1);
        swap_i(s->pss, s->ipss, s->ipss[s->prob->vphotos[k]->id], s->n_p++);
        swap_i(s->vss, s->ivss, s->ivss[s->prob->vphotos[k]->oid], s->n_v++);
    }
    /* solution not evaluated yet */
    s->evalv = s->evalLB = 0;
    s->numEnumLim = st_size(s);
    s->sampleWORLim = s->sampleEnumLim = nh_size(s);
    return s;
}

/*
 * Heuristically constructs a photo slideshow solution
 * Status: INTERIM
 * Notes:
 *   Implement heuristic solution construction
 */
struct solution *heuristicSolution(struct solution *s)
{
    if (s->n_p >= s->prob->p) /* solution s is complete, return it */
        return s;
    fprintf(stderr, "heuristicSolution() not implemented yet for photo slideshow problem.\n");
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
    if (s->n_v % 2) /* solution is unfeasible, cannot evaluate it */
        return NULL;
    int p1, p2, *tag1, n1, *tag2, n2, n3, bit, n = s->i_f + s->n_p, i = s->i_f;
    double obj = 0.0;
    if (s->evalv) /* solution s is evaluated */
        *objv = s->objv;
    else { /* solution s is not evaluated */
        /* get first slide */
        p1 = s->pss[i % s->prob->n];
        p2 = s->pss[(i + 1) % s->prob->n];
        if (s->prob->photos[p1].o == 'V') {
            n1 = merge(s->tgs, s->prob->photos[p1].tgs, s->prob->photos[p1].m, s->prob->photos[p2].tgs, s->prob->photos[p2].m);
            tag1 = s->tgs;
            i++;
        }
        else {
            n1 = s->prob->photos[p1].m;
            tag1 = s->prob->photos[p1].tgs;
        }
        /* compute interest factor for each transition */
        for (i = i + 1, bit = 1; i < n; ++i, bit = 1 - bit) {
            /* get next slide */
            p1 = s->pss[i % s->prob->n];
            p2 = s->pss[(i + 1) % s->prob->n];
            if (s->prob->photos[p1].o == 'V') {
                n2 = merge(s->tgs + MAX_TAGS * bit * 2, s->prob->photos[p1].tgs, s->prob->photos[p1].m, s->prob->photos[p2].tgs, s->prob->photos[p2].m);
                tag2 = s->tgs + MAX_TAGS * bit * 2;
                i++;
            }
            else {
                n2 = s->prob->photos[p1].m;
                tag2 = s->prob->photos[p1].tgs;
            }
            n3 = common(tag1, n1, tag2, n2); /* number of common tags */
            obj += min(n1 - n3, n3, n2 - n3);
            /* update previous slide */
            n1 = n2;
            tag1 = tag2;
        }
        *objv = s->objv = s->prob->maxUB - obj;
        s->evalv = 1;
    }
    return objv;
}

/*
 * Lower bound evaluation
 * Status: INTERIM
 * Notes:
 *   Implement a tighter upper bound for photo slideshow solution
 */
double *getObjectiveLB(double *objLB, struct solution *s)
{
    if (s->evalLB) /* solution s is evaluated */
        *objLB = s->objLB;
    else { /* solution s is not evaluated */
        *objLB = s->objLB = s->prob->maxUB - ub(s);
        s->evalLB = 1;
    }
    return objLB;
}

/*
 * Return the number of components of a solution that are in a given state
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
long getNumSolutionComponents(const struct solution *s, const enum ComponentState st)
{
    switch (st) {
    case PRESENT:
        return s->n_p - 1 > 0 ? s->n_p - 1 : 0;
    case FORBIDDEN:
        fprintf(stderr, "Forbidden state not implemented for getNumSolutionComponents().\n");
        break;
    case UNDEFINED:
        return s->n_p - 1 > 0 ? s->prob->e - (s->n_p - 1) : s->prob->e;
    default:
        fprintf(stderr, "Invalid state passed to getNumSolutionComponents().\n");
        break;
    }
    return -1;
}

/*
 * Return true if all slides in the photo slideshow are complete or false if
 * there is a slide with only one vertical photo
 * Status: FINAL
 */
int isFeasible(struct solution *s)
{
    return !(s->n_v % 2);
}

/*
 * Return true if no photos can be added to the photo slideshow or false if
 * there are photos that can be added to it
 * Status: FINAL
 */
int isComplete(struct solution *s)
{
    return s->n_p == s->prob->p;
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
    int a, b;
    switch (nh) {
    case ADD:
        /* Compute the number of vertical photos not included in the photo slideshow if
         * the number of vertical photos in the solution is odd or the number of photos
         * not included in the photo slideshow that may be added to both extremes if the
         * number of vertical photos in the solution is even.
         */
        a = s->n_p > 1 ? 2 : 1;
        b = s->n_p < s->prob->p;
        return s->n_v % 2 ? s->prob->v - s->n_v : a * (s->prob->n - s->n_p) * b;
    case REMOVE:
        return s->n_p > 2 ? 2 : s->n_p > 1;
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
    int aux;
    switch (st) {
    case PRESENT:
        if (s->numEnumLim.present <= 0)
            return -1;
        aux = s->i_f + s->n_p - s->numEnumLim.present--;
        return pairing(s->pss[aux % s->prob->n], s->pss[(aux - 1) % s->prob->n]);
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
        s->numEnumLim.present = s->n_p - 1 > 0 ? s->n_p - 1 : 0;
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
    int a, b;
    switch (nh) {
    case ADD:
        /* Compute the number of vertical photos not included in the photo slideshow if
         * the number of vertical photos in the solution is odd or the number of photos
         * not included in the photo slideshow that may be added to both extremes if the
         * number of vertical photos in the solution is even.
         */
        a = s->n_p > 1 ? 2 : 1;
        b = s->n_p < s->prob->p;
        s->sampleEnumLim.add = s->n_v % 2 ? s->prob->v - s->n_v : a * (s->prob->n - s->n_p) * b;
        return s;
    case REMOVE:
        s->sampleEnumLim.remove = s->n_p > 2 ? 2 : s->n_p > 1;
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
    int a, b;
    switch (nh) {
    case ADD:
        /* Compute the number of vertical photos not included in the photo slideshow if
         * the number of vertical photos in the solution is odd or the number of photos
         * not included in the photo slideshow that may be added to both extremes if the
         * number of vertical photos in the solution is even.
         */
        a = s->n_p > 1 ? 2 : 1;
        b = s->n_p < s->prob->p;
        s->sampleWORLim.add = s->n_v % 2 ? s->prob->v - s->n_v : a * (s->prob->n - s->n_p) * b;
        return s;
    case REMOVE:
        s->sampleWORLim.remove = s->n_p > 2 ? 2 : s->n_p > 1;
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
    // Coverage: - No valid moves available
    //               - Complete solutions
    if (s->n_p == s->prob->p)
        return NULL;

    // Sample a random valid component
    if (s->n_v % 2) { // Number of vertical photos added to the slideshow is odd [Slide with a vertical photo needs to be completed]

        // Photo is added to complete slide with one vertical photo
        v->side = s->side;
        v->photos[0] = s->pss[(s->i_f + (s->n_p - 1) * v->side) % s->prob->n];
        v->photos[1] = s->prob->vphotos[s->vss[s->n_v + randint(s->prob->v - s->n_v - 1)]]->id;
    }
    else { // Number of vertical photos added to the slideshow is even [No slide with vertical photos needs to be completed]

        int n, i, a = s->n_p > 1 ? 2 : 1, b;

        // Number of the remaining photos not included in the photo slideshow
        n = s->prob->n - s->n_p;

        // Random number to select a photo not included in the photo slideshow that may be added to both ends of the sequence of photos
        i = randint(a * n - 1);

        b = i < n; // Photo is added to { True - Left, False - Right } of the sequence

        v->side = 1 - b;
        v->photos[0] = s->pss[(s->i_f + (s->n_p - 1) * v->side) % s->prob->n];
        v->photos[1] = s->pss[((s->i_f + s->n_p + i - n * v->side) % s->prob->n)];
    }

    // Encode photo transition to single natural number
    v->data = pairing(v->photos[0], v->photos[1]);

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 2 * sizeof(int));
    return v;
}

static struct move *randomRemoveMove(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available
    //               - Empty solutions
    if (s->n_p < 2)
        return NULL;

    // Sample a random valid component
    if (s->n_v % 2) { // Number of vertical photos added to the slideshow is odd [Slide with a vertical photo needs to be removed]

        // Photo is removed from slide
        v->side = s->side;
    }
    else { // Number of vertical photos added to the slideshow is even [No slide with vertical photos needs to be removed]

        int i, a = s->n_p > 2 ? 2 : 1, b;

        // Random number to select a photo included in the photo slideshow that may be removed from both ends of the sequence of photos
        i = randint(a - 1);

        b = i < 1; // Photo is removed from the { True - Left, False - Right } of the sequence

        v->side = 1 - b;
    }

    // Photo is removed
    v->photos[1 - v->side] = s->pss[(s->i_f + (s->n_p - 2) * v->side) % s->prob->n];
    v->photos[v->side] = s->pss[(s->i_f + (s->n_p - 2) * v->side + 1) % s->prob->n];

    // Encode photo transition to single natural number
    v->data = pairing(v->photos[0], v->photos[1]);

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 2 * sizeof(int));
    return v;
}

struct move *randomMove(struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return randomAddMove(v, s);

    case REMOVE:
        return randomRemoveMove(v, s);

    default:
        return NULL; // TODO: Use unimplemented() function.
    }
}

static struct move *enumAddMove(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available
    //               - Complete solutions
    //               - All valid moves are used
    if (!s->sampleEnumLim.add)
        return NULL;

    // Sample a random valid component
    if (s->n_v % 2) { // Number of vertical photos added to the slideshow is odd [Slide with a vertical photo needs to be completed]

        // Photo is added to complete slide with one vertical photo
        v->side = s->side;
        v->photos[0] = s->pss[(s->i_f + (s->n_p - 1) * v->side) % s->prob->n];
        v->photos[1] = s->prob->vphotos[s->vss[s->n_v + s->prob->v - s->n_v - s->sampleEnumLim.add--]]->id;
    }
    else { // Number of vertical photos added to the slideshow is even [No slide with vertical photos needs to be completed]

        int n, i, a = s->n_p > 1 ? 2 : 1, b;

        // Number of the remaining photos not included in the photo slideshow
        n = s->prob->n - s->n_p;

        // Random number to select a photo not included in the photo slideshow that may be added to both ends of the sequence of photos
        i = a * n - s->sampleEnumLim.add--;

        b = i < n; // Photo is added to { True - Left, False - Right } of the sequence

        v->side = 1 - b;
        v->photos[0] = s->pss[(s->i_f + (s->n_p - 1) * v->side) % s->prob->n];
        v->photos[1] = s->pss[((s->i_f + s->n_p + i - n * v->side) % s->prob->n)];
    }

    // Encode photo transition to single natural number
    v->data = pairing(v->photos[0], v->photos[1]);

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 2 * sizeof(int));
    return v;
}

static struct move *enumRemoveMove(struct move *v, struct solution *s)
{
    // Coverage: - No valid moves available
    //               - Empty solutions
    //               - All valid moves are used
    if (!s->sampleEnumLim.remove)
        return NULL;

    // Sample a random valid component
    if (s->n_v % 2) { // Number of vertical photos added to the slideshow is odd [Slide with a vertical photo needs to be removed]

        // Photo is removed from slide
        v->side = s->side;
    }
    else { // Number of vertical photos added to the slideshow is even [No slide with vertical photos needs to be removed]

        int i, a = s->n_p > 2 ? 2 : 1, b;

        // Random number to select a photo included in the photo slideshow that may be removed from both ends of the sequence of photos
        i = a - s->sampleEnumLim.remove--;

        b = i < 1; // Photo is removed from the { True - Left, False - Right } of the sequence

        v->side = 1 - b;
    }

    // Photo is removed
    v->photos[1 - v->side] = s->pss[(s->i_f + (s->n_p - 2) * v->side) % s->prob->n];
    v->photos[v->side] = s->pss[(s->i_f + (s->n_p - 2) * v->side + 1) % s->prob->n];

    // Encode photo transition to single natural number
    v->data = pairing(v->photos[0], v->photos[1]);

    // Update state of evaluation of lower bound increment
    memset(v->evalLBi, 0, 2 * sizeof(int));
    return v;
}

struct move *enumMove(struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return enumAddMove(v, s);

    case REMOVE:
        return enumRemoveMove(v, s);

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
    dest->prob = dest->prob;
    memcpy(dest->pss, src->pss, src->prob->n * sizeof(int));
    memcpy(dest->ipss, src->ipss, src->prob->n * sizeof(int));
    dest->i_f = src->i_f;
    dest->n_p = src->n_p;
    dest->side = src->side;
    memcpy(dest->vss, src->vss, src->prob->v * sizeof(int));
    memcpy(dest->ivss, src->ivss, src->prob->v * sizeof(int));
    dest->n_v = src->n_v;
    dest->numEnumLim = src->numEnumLim;
    dest->sampleEnumLim = src->sampleEnumLim;
    dest->sampleWORLim = src->sampleWORLim;
    dest->evalv = src->evalv;
    dest->objv = src->objv;
    dest->evalLB = src->evalLB;
    dest->objLB = src->objLB;
    return dest;
}

struct solution *applyAddMove(struct solution *s, const struct move *v)
{
    // Coverage: - Solution is complete
    //           - Photo transition is not incident to an end of the sequence of photos in the slideshow
    //           - New photo is already included in photo slideshow
#if __DEBUG__
    if (s->n_p == s->prob->p || v->photos[0] != s->pss[(s->i_f + (s->n_p - 1) * v->side) % s->prob->n] || (s->ipss[v->photos[1]] + s->prob->n - s->i_f) % s->prob->n < (s->n_p + s->prob->n) % s->prob->n) {
        fprintf(stderr, "Component %lu is not undefined, solution is complete, or photo transition |%c(%d)-%c(%d)| is not incident to an end of the sequence of photos in the slideshow.\n"
                        "Could not add component.\n\n",
                v->data, s->prob->photos[v->photos[0]].o, s->prob->photos[v->photos[0]].id, s->prob->photos[v->photos[1]].o, s->prob->photos[v->photos[1]].id);
        return NULL;
    }
#endif

    // Update vertical photos in the slideshow if one was added
    if (s->prob->photos[v->photos[1]].o == 'V') {
#if __DEBUG__
        if (s->ivss[s->prob->photos[v->photos[1]].oid] < s->n_v) {
            fprintf(stderr, "Vertical photo %d is already in the slideshow.\n"
                            "Could not add component.\n\n",
                    v->photos[1]);
            return NULL;
        }
#endif
        swap_i(s->vss, s->ivss, s->ivss[s->prob->photos[v->photos[1]].oid], s->n_v++);
    }

    // Update sequence of photos in the slideshow
    if (v->side) {
        swap_i(s->pss, s->ipss, s->ipss[v->photos[1]], (s->i_f + s->n_p) % s->prob->n);
    }
    else {
        s->i_f = (s->i_f - 1 + s->prob->n) % s->prob->n;
        swap_i(s->pss, s->ipss, s->ipss[v->photos[1]], s->i_f % s->prob->n);
    }

    // Update number of photos added to the photo slideshow
    s->n_p++;

    // Update flag referring to the side where last photo was placed
    s->side = v->side;

    // Update state of evaluation of objective value of incomplete (or incomplete) solution
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
    s->sampleWORLim = s->sampleEnumLim = nh_size(s);

    return s;
}

struct solution *applyRemoveMove(struct solution *s, const struct move *v)
{
    // Coverage: - Solution is empty
    //           - Photo transition is not incident to an end of the sequence of photos in the slideshow
#if __DEBUG__
    if (s->n_p < 2 || v->photos[1 - v->side] != s->pss[(s->i_f + (s->n_p - 2) * v->side) % s->prob->n] || v->photos[v->side] != s->pss[(s->i_f + (s->n_p - 2) * v->side + 1) % s->prob->n]) {
        fprintf(stderr, "Component %lu is not present, solution is empty or photo transition |%c(%d)-%c(%d)| is not incident to an end of the sequence of photos in the slideshow.\n"
                        "Could not remove component.\n\n",
                v->data, s->prob->photos[v->photos[0]].o, s->prob->photos[v->photos[0]].id, s->prob->photos[v->photos[1]].o, s->prob->photos[v->photos[1]].id);
        return NULL;
    }
#endif

    // Update vertical photos in the slideshow if one was removed
    if (s->prob->photos[v->photos[1]].o == 'V') {
#if __DEBUG__
        if (s->ivss[s->prob->photos[v->photos[1]].oid] > s->n_v - 1) {
            fprintf(stderr, "Vertical photo %d is not in the slideshow.\n"
                            "Could not remove component.\n\n",
                    v->photos[1]);
            return NULL;
        }
#endif
        swap_i(s->vss, s->ivss, s->ivss[s->prob->photos[v->photos[1]].oid], --s->n_v);
    }

    // Update sequence of photos in the slideshow
    if (v->side) {
        swap_i(s->pss, s->ipss, s->ipss[v->photos[1]], (s->i_f + s->n_p - 1) % s->prob->n);
    }
    else {
        swap_i(s->pss, s->ipss, s->ipss[v->photos[1]], s->i_f % s->prob->n);
        s->i_f = (s->i_f + 1 + s->prob->n) % s->prob->n;
    }

    // Update number of photos added to the photo slideshow
    s->n_p--;

    // Update flag referring to the side where last photo was removed
    s->side = v->side;

    // Update state of evaluation of objective value of incomplete (or incomplete) solution
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
    s->sampleWORLim = s->sampleEnumLim = nh_size(s);

    return s;
}

struct solution *applyMove(struct solution *s, const struct move *v, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return applyAddMove(s, v);

    case REMOVE:
        return applyRemoveMove(s, v);

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
    memcpy(dest->photos, src->photos, 2 * sizeof(int));
    dest->side = src->side;
    memcpy(dest->evalLBi, src->evalLBi, 2 * sizeof(int));
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
 *   Implement a tighter upper bound for photo slideshow solution
 */
double *getObjectiveLBIncrement(double *obji, struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    int i;
    switch (nh) {
    case ADD:
#if __DEBUG__
        if (s->n_p >= s->prob->p || v->photos[0] != s->pss[(s->i_f + (s->n_p - 1) * v->side) % s->prob->n] || (s->ipss[v->photos[1]] - s->i_f + s->prob->n) % s->prob->n < s->n_p) {
            fprintf(stderr, "Invalid move passed to getObjectiveLBIncrement(). It does not belong to Add neighbourhood.\n");
            return NULL;
        }
#endif
        i = 0;
        break;
    case REMOVE:
#if __DEBUG__
        int aux = s->i_f + (s->n_p - 2) * v->side;
        if (s->n_p < 2 || v->photos[1 - v->side] != s->pss[aux % s->prob->n] || v->photos[v->side] != s->pss[(aux + 1) % s->prob->n]) {
            fprintf(stderr, "Invalid move passed to getObjectiveLBIncrement(). It does not belong to Remove neighbourhood.\n");
            return NULL;
        }
#endif
        i = 1;
        break;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for getObjectiveLBIncrement().\n");
        return NULL;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for getObjectiveLBIncrement().\n");
        return NULL;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to getObjectiveLBIncrement().\n");
        return NULL;
    }
    if (v->evalLBi[i]) /* move v is evaluated */
        *obji = v->objLBi;
    else { /* move v is not evaluated */
        memset(v->evalLBi, 0, sizeof(int) * 2);
        switch (nh) {
        case ADD:
            *obji = v->objLBi = lbi_add(v, s);
            break;
        case REMOVE:
            *obji = v->objLBi = lbi_remove(v, s);
            break;
        default:
            *obji = v->objLBi = DBL_MAX;
            break;
        }
        v->evalLBi[i] = 1;
    }
    return obji;
}
