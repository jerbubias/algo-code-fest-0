/* hsc.c
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

/* Heuristic Solution Construction */

#include "problem.h"
#include <gsl/gsl_rng.h>
#include <time.h>

gsl_rng *rng; /* The single rng instance used by the whole code */

int main(int argc, char **argv)
{
    struct problem *p;
    struct solution *s;
    double cost;

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
        s = allocSolution(p);

        /* Run */
        emptySolution(s);
        if (heuristicSolution(s) != NULL) {
            getObjectiveVector(&cost, s);

            /* Report result */
            printSolution(s);
        }
        else
            fprintf(stderr, "No feasible solution found for %s.\n", argv[0]);

        /* Clean up */
        freeSolution(s);
        freeProblem(p);
    }
    gsl_rng_free(rng);
    return 0;
}
