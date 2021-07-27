/* aco_mmas.c
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

/* Ant Colony Optimisation - Max-Min Ant System */

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
    struct solution *s0, *s1, *iter_best, *global_best, *best, *tmp;
    struct move **v;
    int pop_size, gb_iter, stag, max_iter, n, m, ch, cnt, i, j, k, l;
    double ph_val, ph_wt, cost_wt, p_best, evp_rate, delta, cost, auxcost, itercost, mincost = DBL_MAX, *ph, *enumDes, ph_min = 0.0, ph_max = DBL_MAX, p_dec, sum, avg;

    if (argc < 12) {
        fprintf(stderr, "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <p best> <evp rate> <gb iter> <stag> <delta> <max iter>\n", argv[0]);
        return 0;
    }

    /* Input arguments */
    pop_size = atoi(argv[2]);
    ph_val = atof(argv[3]);
    ph_wt = atof(argv[4]);
    cost_wt = atof(argv[5]);
    p_best = atof(argv[6]);
    evp_rate = atof(argv[7]);
    gb_iter = atoi(argv[8]);
    stag = atoi(argv[9]);
    delta = atof(argv[10]);
    max_iter = atoi(argv[11]);

    /* Input error handling */
    if (pop_size < 1) {
        fprintf(stderr, "Invalid population size. Should be a positive value.\n"
                        "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <p best> <evp rate> <gb iter> <stag> <delta> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (ph_val < 0.0) {
        fprintf(stderr, "Invalid initial pheromone value. Should be nonnegative.\n"
                        "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <evp rate> <Q> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (p_best < 0.0 || p_best > 1.0) {
        fprintf(stderr, "Invalid probability of constructing best solution. Should be included in the unit interval [0,1].\n"
                        "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <p best> <evp rate> <gb iter> <stag> <delta> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (evp_rate <= 0.0 || evp_rate > 1.0) {
        fprintf(stderr, "Invalid evaporation rate value. Should be included in the unit interval (0,1].\n"
                        "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <p best> <evp rate> <gb iter> <stag> <delta> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (gb_iter < 1) {
        fprintf(stderr, "Invalid frequency at which the global best solution updates the pheromones. Should be a positive value.\n"
                        "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <p best> <evp rate> <gb iter> <stag> <delta> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (stag < 0) {
        fprintf(stderr, "Invalid number of iterations until stagnation. Should be a nonnegative value.\n"
                        "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <p best> <evp rate> <gb iter> <stag> <delta> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (delta < 0.0 || delta > 1.0) {
        fprintf(stderr, "Invalid smoothing rate value. Should be included in the unit interval [0,1].\n"
                        "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <p best> <evp rate> <gb iter> <stag> <delta> <max iter>\n",
                argv[0]);
        return 0;
    }
    if (max_iter < 1) {
        fprintf(stderr, "Invalid number of iterations. Should be a positive value.\n"
                        "Usage: %s <file name> <pop size> <ph val> <ph wt> <cost wt> <p best> <evp rate> <gb iter> <stag> <delta> <max iter>\n",
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
        iter_best = allocSolution(p);
        global_best = allocSolution(p);
        v = malloc((n + 1) * sizeof(struct move *));
        for (i = 0; i <= n; ++i)
            v[i] = allocMove(p);
        ph = malloc(m * sizeof(double));
        enumDes = malloc(n * sizeof(double));

        /* Run */

        /* Probability of "right" decision */
        p_dec = pow(p_best, 1.0 / getMaxSolutionSize(p));

        /* Pheromone initialisation */
        for (i = 0; i < m; ++i)
            ph[i] = ph_val;

        for (i = 0, avg = 0.0, ch = 0, cnt = 0; i < max_iter; ++i, ++cnt) {
#if __DEBUG__
            printf("Pheromone Array:\nComponent\t");
            for (j = 0; j < m; ++j)
                printf("%d\t", j);
            printf("\nPheromone Value\t");
            for (j = 0; j < m; ++j)
                printf("%.4lf\t", ph[j]);
            printf("\n\n");
#endif
            for (j = 0, itercost = DBL_MAX; j < pop_size; ++j) {
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
                    if (!k) {
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
                        avg += (getNeighbourhoodSize(s1, ADD) - avg) / ++ch;
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
                if (cost < itercost) { /* avoid copying */
                    tmp = iter_best;
                    iter_best = s0;
                    s0 = tmp;
                    itercost = cost;
                    if (itercost < mincost) {
                        copySolution(global_best, iter_best);
                        mincost = itercost;
                        printSolution(global_best);

                        /* Pheromone limits update */
                        ph_max = 1.0 / (evp_rate * mincost);
                        ph_min = ph_max * (1.0 - p_dec) / ((avg - 1.0) * p_dec);
#if __DEBUG__
                        printf("Pheromone Limits Were Updated.\n"
                               "Lower Pheromone Limit = %.4lf\n"
                               "Upper Pheromone Limit = %.4lf\n\n",
                               ph_min, ph_max);
#endif
                        cnt = -1;
                    }
                }
            }

            /* Pheromone evaporation */
            for (j = 0; j < m; ++j) {
                ph[j] *= 1.0 - evp_rate;
                if (ph[j] < ph_min)
                    ph[j] = ph_min;
            }

            /* Amount of pheromone laid */
            best = (i + 1) % gb_iter ? iter_best : global_best;
            cost = (i + 1) % gb_iter ? itercost : mincost;
            resetEnumSolutionComponents(best, PRESENT);
            while ((j = enumSolutionComponents(best, PRESENT)) > -1) {
                ph[j] += 1.0 / cost;
                if (ph[j] > ph_max)
                    ph[j] = ph_max;
            }
#if __DEBUG__
            if ((i + 1) % gb_iter)
                printf("Pheromones Updated With Iteration-Best Solution.\n\n");
            else
                printf("Pheromones Updated With Global-Best Solution.\n\n");
            printSolution(best);
#endif

            /* Pheromone reinitialisation */
            if (cnt >= stag) {
#if __DEBUG__
                printf("Stagnation Was Reached. Smoothing Pheromone Trails...\n\n");
#endif
                for (j = 0; j < m; ++j)
                    ph[j] += delta * (ph_max - ph[j]);
            }
        }

        /* Report result */
        printSolution(global_best);

        /* Clean up */
        freeSolution(s0);
        freeSolution(s1);
        freeSolution(iter_best);
        freeSolution(global_best);
        for (i = 0; i <= n; ++i)
            freeMove(v[i]);
        free(v);
        free(ph);
        free(enumDes);
        freeProblem(p);
    }
    gsl_rng_free(rng);
    return 0;
}