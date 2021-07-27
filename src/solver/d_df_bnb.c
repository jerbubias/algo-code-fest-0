/* d_df_bnb.c
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

/* Dichotomic Depth-First Branch-and-Bound */

#include "problem.h"
#include <float.h>
#include <gsl/gsl_rng.h>
#include <time.h>

gsl_rng *rng; /* The single rng instance used by the whole code */

struct solution *expandSolution(double *objUB, struct solution **best, struct solution *s0, struct solution **s1, struct move **v, int i)
{
    struct solution *tmp;
    double obj;

#if __DEBUG__
    printf("---------- Depth Level %d  (Bound) ----------\n\n", i);
#endif
    getObjectiveLB(&obj, s0);
#if __DEBUG__
    printSolution(s0);
#endif

    /* Bounding */
#if __DEBUG__
    printf("Lower Bound Value = %.1lf, Upper Bound Value = %.1lf\n\n", obj, *objUB);
#endif
    if (obj > *objUB)
        return s0;

    /* Heuristic search */
    copySolution(*s1, s0);
    if (heuristicSolution(*s1) != NULL) {
        getObjectiveVector(&obj, *s1);
        if (obj < *objUB) { /* avoid copying */
            tmp = *best;
            *best = *s1;
            *s1 = tmp;
            *objUB = obj;
            printSolution(*best);
        }
    }

    if (heuristicMove(v[i], s0, ADD) != NULL) { /* non-terminal node */
#if __DEBUG__
        printf("---------- Depth Level %d    (Add) ----------\n\n", i);
#endif
        getObjectiveLBIncrement(&obj, v[i], s0, ADD);
#if __DEBUG__
        printMove(v[i]);
#endif
        if (applyMove(s0, v[i], ADD) == NULL)
            exit(EXIT_FAILURE);
#if __DEBUG__
        printSolution(s0);
#endif

        expandSolution(objUB, best, s0, s1, v, i + 1);

#if __DEBUG__
        printf("---------- Depth Level %d (Remove) ----------\n\n", i);
#endif
        getObjectiveLBIncrement(&obj, v[i], s0, REMOVE);
#if __DEBUG__
        printMove(v[i]);
#endif
        if (applyMove(s0, v[i], REMOVE) == NULL)
            exit(EXIT_FAILURE);
#if __DEBUG__
        printSolution(s0);
#endif

#if __DEBUG__
        printf("---------- Depth Level %d (Forbid) ----------\n\n", i);
#endif
        getObjectiveLBIncrement(&obj, v[i], s0, FORBID);
#if __DEBUG__
        printMove(v[i]);
#endif
        if (applyMove(s0, v[i], FORBID) == NULL)
            exit(EXIT_FAILURE);
#if __DEBUG__
        printSolution(s0);
#endif

        expandSolution(objUB, best, s0, s1, v, i + 1);

#if __DEBUG__
        printf("---------- Depth Level %d (Permit) ----------\n\n", i);
#endif
        getObjectiveLBIncrement(&obj, v[i], s0, PERMIT);
#if __DEBUG__
        printMove(v[i]);
#endif
        if (applyMove(s0, v[i], PERMIT) == NULL)
            exit(EXIT_FAILURE);
#if __DEBUG__
        printSolution(s0);
#endif
    }
    else { /* terminal node */
        getObjectiveVector(&obj, s0);
        if (obj < *objUB) {
            copySolution(*best, s0);
            *objUB = obj;
            printSolution(*best);
        }
    }
    return s0;
}

int main(int argc, char **argv)
{
    struct problem *p;
    struct solution *s0, *s1, *best;
    struct move **v;
    int m, i;
    double objUB = DBL_MAX;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <file name>\n", argv[0]);
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
        m = getNumComponents(p);
        s0 = allocSolution(p);
        s1 = allocSolution(p);
        best = allocSolution(p);
        v = malloc((m + 1) * sizeof(struct move *));
        for (i = 0; i <= m; ++i)
            v[i] = allocMove(p);

        /* Run */
        emptySolution(s0);
        expandSolution(&objUB, &best, s0, &s1, v, 0);

        /* Report result */
        printSolution(best);

        /* Clean up */
        freeSolution(s0);
        freeSolution(s1);
        freeSolution(best);
        for (i = 0; i <= m; ++i)
            freeMove(v[i]);
        free(v);
        freeProblem(p);
    }
    gsl_rng_free(rng);
    return 0;
}
