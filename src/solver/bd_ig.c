/* bd_ig.c
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

/* Iterated Greedy With Biased Destruction */

/* Inspired on the pseudocode (Figure 3) in Iterated Greedy (Handbook of Heuristics) by
 * Rubén Ruiz and Thomas Stützle
 */

#include "problem.h"
#include <gsl/gsl_rng.h>
#include <math.h>
#include <time.h>

gsl_rng *rng; /* The single rng instance used by the whole code */

int rouletteWheelSelection(double *val, int size, double sum) /* TODO: API for selection methods */
{
    int i;
    double r, cumsum = 0.0;

    r = gsl_rng_uniform(rng);
    for (i = 0; i < size; ++i) {
        cumsum += val[i] / sum;
        if (r < cumsum)
            return i;
    }
    return i - 1;
}

int main(int argc, char **argv)
{
    struct problem *p;
    struct solution *s0, *s1, *best, *tmp;
    struct move **v;
    int d, max_iter, n, i, j, k, l;
    double T, cost, auxcost, mincost, *enumCost, sum;

    if (argc < 5) {
        fprintf(stderr, "Usage: %s <file name> <d> <T> <max iter>\n", argv[0]);
        return 0;
    }

    /* Input arguments */
    d = atoi(argv[2]);
    T = atof(argv[3]);
    max_iter = atoi(argv[4]);

    /* Input error handling */
    if (d < 1) {
        fprintf(stderr, "Invalid destruction size. Should be a positive value.\n"
                        "Usage: %s <file name> <d> <T> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (T <= 0.0) {
        fprintf(stderr, "Invalid temperature value. Should be positive.\n"
                        "Usage: %s <file name> <d> <T> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (max_iter < 1) {
        fprintf(stderr, "Invalid number of iterations. Should be a positive value.\n"
                        "Usage: %s <file name> <d> <T> <max iter>\n",
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
        n = getMaxNeighbourhoodSize(p, REMOVE);
        s0 = allocSolution(p);
        s1 = allocSolution(p);
        best = allocSolution(p);
        v = malloc((n + 1) * sizeof(struct move *));
        for (i = 0; i <= n; ++i)
            v[i] = allocMove(p);
        enumCost = malloc(n * sizeof(double));

        /* Run */
        emptySolution(s0);

        /* Heuristic search */
        if (heuristicSolution(s0) != NULL) {
            getObjectiveVector(&auxcost, s0);
            copySolution(best, s0);
            mincost = auxcost;
            printSolution(best);

            for (i = 0; i < max_iter; ++i) {
#if __DEBUG__
                printf("----------------- Iteration %d -----------------\n\n", i + 1);
#endif
                copySolution(s1, s0);
#if __DEBUG__
                printSolution(s1);
#endif

                /* Biased Destruction */
                for (j = 0; j < d; ++j) {
                    k = 0;
                    sum = 0.0;
                    while (enumMove(v[k], s1, REMOVE) != NULL) {
                        getObjectiveLBIncrement(&cost, v[k], s1, REMOVE);
                        enumCost[k] = -cost;
                        sum += enumCost[k++];
                    }
#if __DEBUG__
                    printf("Enumerated Components:\nComponent\t");
                    for (l = 0; l < k; ++l)
                        printf("%ld\t", getComponentFromMove(v[l]));
                    printf("\nCost\t\t");
                    for (l = 0; l < k; ++l)
                        printf("%.1lf\t", enumCost[l]);
                    printf("\nProbability\t");
                    for (l = 0; l < k; ++l)
                        printf("%3.2lf%%\t", sum > 0.0 ? enumCost[l] * 100.0 / sum : 100.0 / k);
                    printf("\n\n");
#endif
                    if (k > 0) {
                        l = sum > 0.0 ? rouletteWheelSelection(enumCost, k, sum) : (int)gsl_rng_uniform_int(rng, k);
#if __DEBUG__
                        printMove(v[l]);
#endif
                        if (applyMove(s1, v[l], REMOVE) == NULL)
                            exit(EXIT_FAILURE);
#if __DEBUG__
                        printSolution(s1);
#endif
                    }
                    else
                        break;
                }

                /* Construction */
                if (heuristicSolution(s1) != NULL) {

                    /* Acceptance criterion */
                    getObjectiveVector(&cost, s1);
#if __DEBUG__
                    printSolution(s1);
#endif
                    if (cost < auxcost) { /* avoid copying */
                        tmp = s0;
                        s0 = s1;
                        s1 = tmp;
                        auxcost = cost;
                        if (auxcost < mincost) {
                            copySolution(best, s0);
                            mincost = auxcost;
                            printSolution(best);
                        }
                    }
                    else if (gsl_rng_uniform(rng) < exp((auxcost - cost) / T)) { /* avoid copying */
#if __DEBUG__
                        printf("Accepted Worse Solution With a Probability of %3.2lf%%\n\n", exp((auxcost - cost) / T) * 100.0);
#endif
                        tmp = s0;
                        s0 = s1;
                        s1 = tmp;
                        auxcost = cost;
                    }
                }
            }

            /* Report result */
            printSolution(best);
        }
        else
            fprintf(stderr, "No feasible solution found for %s.\n", argv[0]);

        /* Clean up */
        freeSolution(s0);
        freeSolution(s1);
        freeSolution(best);
        for (i = 0; i <= n; ++i)
            freeMove(v[i]);
        free(v);
        free(enumCost);
        freeProblem(p);
    }
    gsl_rng_free(rng);
    return 0;
}