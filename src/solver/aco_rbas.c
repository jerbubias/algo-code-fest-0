/* aco_rbas.c
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

/* Ant Colony Optimisation - Rank-Based Ant System */

#include "minmax_heap.h"
#include "problem.h"
#include <float.h>
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
    struct solution *s0, *s1, **ant, *best, *tmp;
    struct move **v;
    struct heap *h;
    int pop_size, w, max_iter, n, m, cnt, i, j, k, l;
    double ph_val, ph_wt, cost_wt, evp_rate, cost, auxcost, mincost = DBL_MAX, *ph, *enumDes, *antCost, sum;

    if (argc < 9) {
        fprintf(stderr, "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <evp rate> <w> <max iter>\n", argv[0]);
        return 0;
    }

    /* Input arguments */
    pop_size = atoi(argv[2]);
    ph_val = atof(argv[3]);
    ph_wt = atof(argv[4]);
    cost_wt = atof(argv[5]);
    evp_rate = atof(argv[6]);
    w = atoi(argv[7]);
    max_iter = atoi(argv[8]);

    /* Input error handling */
    if (pop_size < 1) {
        fprintf(stderr, "Invalid population size. Should be a positive value.\n"
                        "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <evp rate> <w> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (ph_val < 0.0) {
        fprintf(stderr, "Invalid initial pheromone value. Should be nonnegative.\n"
                        "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <evp rate> <w> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (evp_rate <= 0.0 || evp_rate > 1.0) {
        fprintf(stderr, "Invalid evaporation rate value. Should be included in the unit interval (0,1].\n"
                        "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <evp rate> <w> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (w < 1 || w > pop_size + 1) {
        fprintf(stderr, "Invalid number of ants that lay pheromones. Should be a positive value, no greater than the size of the population with one additional ant.\n"
                        "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <evp rate> <w> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (max_iter < 1) {
        fprintf(stderr, "Invalid number of iterations. Should be a positive value.\n"
                        "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <evp rate> <w> <max iter>\n",
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
        m = getNumComponents(p);
        s0 = allocSolution(p);
        s1 = allocSolution(p);
        ant = malloc((w - 1) * sizeof(struct solution *));
        for (i = 0; i < w - 1; ++i)
            ant[i] = allocSolution(p);
        best = allocSolution(p);
        v = malloc((n + 1) * sizeof(struct move *));
        for (i = 0; i <= n; ++i)
            v[i] = allocMove(p);
        h = newHeap(w - 1);
        ph = malloc(m * sizeof(double));
        enumDes = malloc(n * sizeof(double));
        antCost = malloc((w - 1) * sizeof(double));

        /* Run */

        /* Pheromone initialisation */
        for (i = 0; i < m; ++i)
            ph[i] = ph_val;

        for (i = 0; i < max_iter; ++i) {
#if __DEBUG__
            printf("Pheromone Array:\nComponent\t");
            for (j = 0; j < m; ++j)
                printf("%d\t", j);
            printf("\nPheromone Value\t");
            for (j = 0; j < m; ++j)
                printf("%.4lf\t", ph[j]);
            printf("\n\n");
#endif

            for (j = 0, cnt = 0; j < pop_size; ++j) {
#if __DEBUG__
                printf("-------- Iteration %d | Individual %d --------\n\n", i + 1, j);
#endif

                /* Ant based solution construction */
                emptySolution(s1);
#if __DEBUG__
                printSolution(s1);
#endif
                auxcost = DBL_MAX;
                for (;;) {
                    k = 0;
                    sum = 0.0;
                    while (enumMove(v[k], s1, ADD) != NULL) {
                        getObjectiveLBIncrement(&cost, v[k], s1, ADD);
                        enumDes[k] = pow(ph[getComponentFromMove(v[k])], ph_wt) * pow((1.0 / (cost + 1.0)), cost_wt);
                        sum += enumDes[k++];
                    }
                    if (!k && isFeasible(s1))
                        break;
                    else if (!k) {
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
                        for (l = 0; l < k; ++l)
                            printf("%ld\t", getComponentFromMove(v[l]));
                        printf("\nCost\t\t");
                        for (l = 0; l < k; ++l) {
                            getObjectiveLBIncrement(&cost, v[l], s1, ADD);
                            printf("%.1lf\t", cost);
                        }
                        printf("\nPheromone Value\t");
                        for (l = 0; l < k; ++l)
                            printf("%.4lf\t", ph[getComponentFromMove(v[l])]);
                        printf("\nDesirability\t");
                        for (l = 0; l < k; ++l)
                            printf("%.4lf\t", enumDes[l]);
                        printf("\nSum of Desirability for Enumerated Components = %.12lf\n\n", sum);
#endif
                        l = rouletteWheelSelection(enumDes, k, sum);
#if __DEBUG__
                        printMove(v[l]);
#endif
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
                if (cost < mincost) {
                    copySolution(best, s0);
                    mincost = cost;
                    printSolution(best);
                }

                /* Ant ranking */
                if (cnt < w - 1) { /* avoid copying */
                    tmp = ant[cnt];
                    ant[cnt] = s0;
                    s0 = tmp;
                    antCost[cnt++] = cost;

                    /* Build full heap */
                    if (cnt == w - 1)
                        buildHeap(h, antCost, w - 1);
                }
                else if (findHeapMax(h) > cost) { /* avoid copying */
                    k = deleteHeapMax(h);
                    tmp = ant[k];
                    ant[k] = s0;
                    s0 = tmp;
                    antCost[k] = cost;
                    insertHeap(h, k, cost);
                }
            }

            /* Pheromone update */
            for (j = 0; j < m; ++j)
                ph[j] = (1.0 - evp_rate) * ph[j];
#if __DEBUG__
            if (w - 1 > 0)
                printf("Ant Ranking (Top %d)\n\n", w - 1);
            auxcost = -DBL_MAX;
#endif
            for (j = 0; j < w - 1; ++j) {
                k = deleteHeapMin(h);
#if __DEBUG__
                printf("Position: %d\n", j + 1);
                printSolution(ant[k]);
                if (auxcost > antCost[k])
                    exit(EXIT_FAILURE);
                auxcost = antCost[k];
#endif
                while ((l = enumSolutionComponents(ant[k], PRESENT)) > -1)
                    ph[l] += (w - j - 1.0) / antCost[k];
            }
            while ((j = enumSolutionComponents(best, PRESENT)) > -1)
                ph[j] += w / mincost;
        }

        /* Report result */
        printSolution(best);

        /* Clean up */
        freeSolution(s0);
        freeSolution(s1);
        for (i = 0; i < w - 1; ++i)
            freeSolution(ant[i]);
        free(ant);
        freeSolution(best);
        for (i = 0; i <= n; ++i)
            freeMove(v[i]);
        free(v);
        freeHeap(h);
        free(ph);
        free(enumDes);
        free(antCost);
        freeProblem(p);
    }
    gsl_rng_free(rng);
    return 0;
}