/* test.c
 *
 * (C) 2021 Andreia P. Guerreiro <andreia.guerreiro@tecnico.ulisboa.pt>
 *          Carlos M. Fonseca <cmfonsec@dei.uc.pt>
 *          Samuel B. Outeiro <souteiro@student.dei.uc.pt>
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

#include "init/init.h"
#include "problem.h"
#include "property.h"
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_REP 10 /* number of times each test is repeated */

gsl_rng *rng; /* The single rng instance used by the whole code */

int main(int argc, char **argv)
{
    struct problem *p;

    /* Set up random number generation */
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(0));

    /* Problem instantiation */
    p = initProblem(argc, argv);

    if (p != NULL) {
        /* Run tests */
        runPropertyBased(p, MAX_REP);

        /* Clean up */
        freeProblem(p);
    }
    gsl_rng_free(rng);
    return 0;
}
