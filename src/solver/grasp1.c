/* grasp1.c
 *
 * (C) 2021 Samuel B. Outeiro <souteiro@student.dei.uc.pt>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License, version 3, as
 * published by the Free Software Foundation.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Greedy Randomized Adaptive Search Procedure */

/* Inspired on the pseudocode (Algorithm 108) in Essentials of Metaheuristics by Sean Luke
 * https://cs.gmu.edu/~sean/book/metaheuristics/Essentials.pdf
 */

#include "problem.h"
#include "sort_r.h"
#include <float.h>
#include <gsl/gsl_rng.h>
#include <time.h>

gsl_rng *rng; /* The single rng instance used by the whole code */

static int compar_r(const void *_a, const void *_b, void *_arg)
{
    const int a = *(const int *)_a;
    const int b = *(const int *)_b;
    const double *arg = (const double *)_arg;
    return *(arg + a) < *(arg + b);
}

int main(int argc, char **argv)
{
    struct problem *p;
    struct solution *s0, *s1, *best, *tmp;
    struct move **v;
    int max_iter, n, *enumOrder, i, j, k, l;
    double pr, cost, auxcost, mincost = DBL_MAX, *enumCost;

    if (argc < 4) {
        fprintf(stderr, "Usage: %s <file name> <pr> <max iter>\n", argv[0]);
        return 0;
    }

    /* Input arguments */
    pr = atof(argv[2]);
    max_iter = atoi(argv[3]);

    /* Input error handling */
    if (pr < 0.0 || pr > 1.0) {
        fprintf(stderr, "Invalid value for the proportion of components for restricted candidate list. Should be included in the unit interval [0,1].\n"
                        "Usage: %s <file name> <pr> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (max_iter < 1) {
        fprintf(stderr, "Invalid number of iterations. Should be a positive value.\n"
                        "Usage: %s <file name> <pr> <max iter>\n",
                argv[0]);
        return 0;
    }

    /* Set up random number generation */
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(0));

    /* Problem instantiation */
    p = newProblem(argv[1]);
    if (p != NULL) {
        printProblem(p);

        /* Memory allocation */
        n = getMaxNeighbourhoodSize(p, ADD);
        s0 = allocSolution(p);
        s1 = allocSolution(p);
        best = allocSolution(p);
        v = malloc((n + 1) * sizeof(struct move *));
        for (i = 0; i <= n; ++i)
            v[i] = allocMove(p);
        enumCost = malloc((n + 1) * sizeof(double));
        enumOrder = malloc((n + 1) * sizeof(int));

        /* Run */
        for (i = 0; i < max_iter; ++i) {
#if __DEBUG__
            printf("----------------- Iteration %d -----------------\n\n", i + 1);
#endif

            /* Construct greedy randomized solution */
            emptySolution(s1);
#if __DEBUG__
            printSolution(s1);
#endif
            auxcost = DBL_MAX;
            for (;;) {
                j = 0;
                while (enumMove(v[j], s1, ADD) != NULL) {
                    getObjectiveLBIncrement(&cost, v[j], s1, ADD);
                    enumOrder[j] = j;
                    enumCost[j++] = cost;
                }
                if (!j && isFeasible(s1))
                    break;
                else if (!j) {
                    emptySolution(s1);
#if __DEBUG__
                    printSolution(s1);
#endif
                    auxcost = DBL_MAX;
                    continue;
                }
                else {
#if __DEBUG__
                    printf("Enumerated Components:\nComponent\t");
                    for (k = 0; k < j; ++k)
                        printf("%ld\t", getComponentFromMove(v[k]));
                    printf("\nCost\t\t");
                    for (k = 0; k < j; ++k) {
                        getObjectiveLBIncrement(&cost, v[k], s1, ADD);
                        printf("%.1lf\t", cost);
                    }
                    printf("\n\n");
#endif
                    sort_r(enumOrder, j, sizeof(int), compar_r, enumCost);
                    k = pr * j > 1 ? pr * j + 0.5 : 1;
                    for (; k < j && enumCost[enumOrder[k - 1]] == enumCost[enumOrder[k]]; ++k);
                    l = gsl_rng_uniform_int(rng, k);
                    if (applyMove(s1, v[l], ADD) == NULL)
                        exit(EXIT_FAILURE);
#if __DEBUG__
                    printSolution(s1);
#endif
                    if (getObjectiveVector(&cost, s1) != NULL && cost < auxcost) {
                        copySolution(s0, s1);
                        auxcost = cost;
                    }
                }
            }

            getObjectiveVector(&cost, s0);
            if (cost < mincost) { /* avoid copying */
                tmp = best;
                best = s0;
                s0 = tmp;
                mincost = cost;
                printSolution(best);
            }
        }

        /* Report result */
        printSolution(best);

        /* Clean up */
        freeSolution(s0);
        freeSolution(s1);
        freeSolution(best);
        for (i = 0; i <= n; ++i)
            freeMove(v[i]);
        free(v);
        free(enumCost);
        free(enumOrder);
        freeProblem(p);
    }
    gsl_rng_free(rng);
    return 0;
}