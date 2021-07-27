/* aco_acs.c
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

/* Ant Colony Optimisation - Ant Colony System */

/* Inspired on the pseudocode (Algorithm 112) in Essentials of Metaheuristics by Sean Luke
 * https://cs.gmu.edu/~sean/book/metaheuristics/Essentials.pdf
 */

#include "problem.h"
#include <float.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>
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
    int pop_size, max_iter, n, m, i, j, k, l, r;
    double ph_val, cost_wt, q, local_evp_rate, global_evp_rate, cost, auxcost, mincost = DBL_MAX, *phOff, *phOn, *enumDes, sum, aux;

    if (argc < 9) {
        fprintf(stderr, "Usage: %s <file name> <pop size> <ph val> <cost wt> <q> <local evp rate> <global evp rate> <max iter>\n", argv[0]);
        return 0;
    }

    /* Input arguments */
    pop_size = atoi(argv[2]);
    ph_val = atof(argv[3]);
    cost_wt = atof(argv[4]);
    q = atof(argv[5]);
    local_evp_rate = atof(argv[6]);
    global_evp_rate = atof(argv[7]);
    max_iter = atoi(argv[8]);

    /* Input error handling */
    if (pop_size < 1) {
        fprintf(stderr, "Invalid population size. Should be a positive value.\n"
                        "Usage: %s <file name> <pop size> <ph val> <cost wt> <q> <local evp rate> <global evp rate> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (ph_val < 0.0) {
        fprintf(stderr, "Invalid initial pheromone value. Should be nonnegative.\n"
                        "Usage: %s <file name> <pop size> <ph val> <cost wt> <q> <local evp rate> <global evp rate> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (q < 0.0 || q > 1.0) {
        fprintf(stderr, "Invalid probability for the Pseudorandom Proportional Action Choice Rule. Should be included in the unit interval [0,1].\n"
                        "Usage: %s <file name> <pop size> <ph val> <cost wt> <q> <local evp rate> <global evp rate> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (local_evp_rate <= 0.0 || local_evp_rate >= 1.0) {
        fprintf(stderr, "Invalid local evaporation rate value. Should be included in the unit interval (0,1).\n"
                        "Usage: %s <file name> <pop size> <ph val> <cost wt> <q> <local evp rate> <global evp rate> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (global_evp_rate <= 0.0 || global_evp_rate > 1.0) {
        fprintf(stderr, "Invalid global evaporation rate value. Should be included in the unit interval (0,1].\n"
                        "Usage: %s <file name> <pop size> <ph val> <cost wt> <q> <local evp rate> <global evp rate> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (max_iter < 1) {
        fprintf(stderr, "Invalid number of iterations. Should be a positive value.\n"
                        "Usage: %s <file name> <pop size> <ph val> <cost wt> <q> <evp rate> <max iter>\n",
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
        best = allocSolution(p);
        v = malloc((n + 1) * sizeof(struct move *));
        for (i = 0; i <= n; ++i)
            v[i] = allocMove(p);
        phOff = malloc(m * sizeof(double));
        phOn = malloc(m * sizeof(double));
        enumDes = malloc(n * sizeof(double));

        /* Run */

        /* Global pheromone trail initialisation */
        for (i = 0; i < m; ++i)
            phOff[i] = ph_val;

        for (i = 0; i < max_iter; ++i) {
#if __DEBUG__
            printf("Pheromone Array:\nComponent\t");
            for (j = 0; j < m; ++j)
                printf("%d\t", j);
            printf("\nPheromone Value\t");
            for (j = 0; j < m; ++j)
                printf("%.4lf\t", phOff[j]);
            printf("\n\n");
#endif

            /* Local pheromone trail initialisation */
            memcpy(phOn, phOff, m * sizeof(double));

            for (j = 0; j < pop_size; ++j) {
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
                    k = l = 0;
                    sum = aux = 0.0;
                    while (enumMove(v[k], s1, ADD) != NULL) {
                        getObjectiveLBIncrement(&cost, v[k], s1, ADD);
                        enumDes[k] = phOn[getComponentFromMove(v[k])] * pow((1.0 / (cost + 1.0)), cost_wt);
                        sum += enumDes[k];
                        if (enumDes[k] > aux) {
                            l = k;
                            aux = enumDes[k];
                        }
                        ++k;
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
                        for (r = 0; r < k; ++r)
                            printf("%ld\t", getComponentFromMove(v[r]));
                        printf("\nCost\t\t");
                        for (r = 0; r < k; ++r) {
                            getObjectiveLBIncrement(&cost, v[r], s1, ADD);
                            printf("%.1lf\t", cost);
                        }
                        printf("\nPheromone Value\t");
                        for (r = 0; r < k; ++r)
                            printf("%.4lf\t", phOn[getComponentFromMove(v[r])]);
                        printf("\nDesirability\t");
                        for (r = 0; r < k; ++r)
                            printf("%.4lf\t", enumDes[r]);
                        printf("\nSum of Desirability for Enumerated Components = %.12lf\n\n", sum);
                        getObjectiveLBIncrement(&cost, v[l], s1, ADD);
                        printf("Highest Component Desirability:\n"
                               "Component\t%ld\n"
                               "Cost\t\t%.1lf\n"
                               "Pheromone Value\t%.4lf\n"
                               "Desirability\t%.4lf\n\n",
                               getComponentFromMove(v[l]), cost, phOn[l], aux);
#endif
                        r = gsl_rng_uniform(rng) < q ? l : rouletteWheelSelection(enumDes, k, sum);
#if __DEBUG__
                        printMove(v[r]);
#endif
                        if (applyMove(s1, v[r], ADD) == NULL)
                            exit(EXIT_FAILURE);
#if __DEBUG__
                        printSolution(s1);
#endif
                        if (getObjectiveVector(&cost, s1) != NULL && cost < auxcost) {
                            copySolution(s0, s1);
                            auxcost = cost;
                        }

                        /* Local pheromone trail update */
                        phOn[getComponentFromMove(v[r])] = (1.0 - local_evp_rate) * phOn[getComponentFromMove(v[r])] + local_evp_rate * ph_val;
#if __DEBUG__
                        printf("Pheromone Array:\nComponent\t");
                        for (j = 0; j < m; ++j)
                            printf("%d\t", j);
                        printf("\nPheromone Value\t");
                        for (j = 0; j < m; ++j)
                            printf("%.4lf\t", phOn[j]);
                        printf("\n\n");
#endif
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

            /* Global pheromone trail update */
            resetEnumSolutionComponents(best, PRESENT);
            while ((j = enumSolutionComponents(best, PRESENT)) > -1)
                phOff[j] = (1.0 - global_evp_rate) * phOff[j] + global_evp_rate * (1.0 / mincost);
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
        free(phOff);
        free(phOn);
        free(enumDes);
        freeProblem(p);
    }
    gsl_rng_free(rng);
    return 0;
}