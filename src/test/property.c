/* property.c
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

#include "problem.h"
#include <gsl/gsl_rng.h>
#include <stdio.h>

/* Text format */
#define OK "\x1b[1;34m"
#define NOTOK "\x1b[1;31m"
#define INFO "\x1b[0;90m"
#define WARN "\x1b[1;92m"
#define BOLD "\033[1m"
#define RESET "\033[0m"

#define MAX_MOVE 10000

extern gsl_rng *rng; /* The single rng instance used by the whole code */

/**********************************/
/* -------- Enumerations -------- */
/**********************************/

enum ProblemFunction { _getNumComponents_,
                       _getMaxSolutionSize_,
                       _getMaxAddNeighbourhoodSize_,
                       _getMaxRemoveNeighbourhoodSize_,
                       _getMaxForbidNeighbourhoodSize_,
                       _getMaxPermitNeighbourhoodSize_,
                       _allocSolution_,
                       _allocMove_,
                       _freeSolution_,
                       _freeMove_,
                       _printProblem_,
                       _printSolution_,
                       _printMove_,
                       _emptySolution_,
                       _heuristicSolution_,
                       _getObjectiveVector_,
                       _getObjectiveLB_,
                       _getNumPresentSolutionComponents_,
                       _getNumUndefinedSolutionComponents_,
                       _getNumForbiddenSolutionComponents_,
                       _isFeasible_,
                       _isComplete_,
                       _isDominant_,
                       _checkDominance_,
                       _getAddNeighbourhoodSize_,
                       _getRemoveNeighbourhoodSize_,
                       _getForbidNeighbourhoodSize_,
                       _getPermitNeighbourhoodSize_,
                       _enumPresentSolutionComponents_,
                       _enumUndefinedSolutionComponents_,
                       _enumForbiddenSolutionComponents_,
                       _resetEnumPresentSolutionComponents_,
                       _resetEnumUndefinedSolutionComponents_,
                       _resetEnumForbiddenSolutionComponents_,
                       _enumAddMove_,
                       _enumRemoveMove_,
                       _enumForbidMove_,
                       _enumPermitMove_,
                       _resetEnumAddMove_,
                       _resetEnumRemoveMove_,
                       _resetEnumForbidMove_,
                       _resetEnumPermitMove_,
                       _heuristicAddMove_,
                       _heuristicRemoveMove_,
                       _heuristicForbidMove_,
                       _heuristicPermitMove_,
                       _heuristicAddMoveWOR_,
                       _heuristicRemoveMoveWOR_,
                       _heuristicForbidMoveWOR_,
                       _heuristicPermitMoveWOR_,
                       _resetHeuristicAddMoveWOR_,
                       _resetHeuristicRemoveMoveWOR_,
                       _resetHeuristicForbidMoveWOR_,
                       _resetHeuristicPermitMoveWOR_,
                       _randomAddMove_,
                       _randomRemoveMove_,
                       _randomForbidMove_,
                       _randomPermitMove_,
                       _randomAddMoveWOR_,
                       _randomRemoveMoveWOR_,
                       _randomForbidMoveWOR_,
                       _randomPermitMoveWOR_,
                       _resetRandomAddMoveWOR_,
                       _resetRandomRemoveMoveWOR_,
                       _resetRandomForbidMoveWOR_,
                       _resetRandomPermitMoveWOR_,
                       _copySolution_,
                       _applyAddMove_,
                       _applyRemoveMove_,
                       _applyForbidMove_,
                       _applyPermitMove_,
                       _copyMove_,
                       _getComponentFromMove_,
                       _getAddHeuristicValue_,
                       _getRemoveHeuristicValue_,
                       _getForbidHeuristicValue_,
                       _getPermitHeuristicValue_,
                       _getAddObjectiveLBIncrement_,
                       _getRemoveObjectiveLBIncrement_,
                       _getForbidObjectiveLBIncrement_,
                       _getPermitObjectiveLBIncrement_ };

/*********************************/
/* ------ Data structures ------ */
/*********************************/

struct property_test {
    char *name;
    char *description;
    int (*function)(struct problem *);
    enum ProblemFunction *pf;
    int pfn;
};

/**********************************/
/* ------ Global variables ------ */
/**********************************/

static int pfmacro[] = {GETNUMCOMPONENTS, GETMAXSOLUTIONSIZE, GETMAXADDNEIGHBOURHOODSIZE, GETMAXREMOVENEIGHBOURHOODSIZE, GETMAXFORBIDNEIGHBOURHOODSIZE,
                        GETMAXPERMITNEIGHBOURHOODSIZE,
                        ALLOCSOLUTION, ALLOCMOVE, FREESOLUTION, FREEMOVE,
                        PRINTPROBLEM, PRINTSOLUTION, PRINTMOVE,
                        EMPTYSOLUTION, HEURISTICSOLUTION,
                        GETOBJECTIVEVECTOR, GETOBJECTIVELB, GETNUMPRESENTSOLUTIONCOMPONENTS, GETNUMUNDEFINEDSOLUTIONCOMPONENTS, GETNUMFORBIDDENSOLUTIONCOMPONENTS,
                        ISFEASIBLE, ISCOMPLETE, ISDOMINANT, CHECKDOMINANCE, GETADDNEIGHBOURHOODSIZE, GETREMOVENEIGHBOURHOODSIZE, GETFORBIDNEIGHBOURHOODSIZE,
                        GETPERMITNEIGHBOURHOODSIZE, ENUMPRESENTSOLUTIONCOMPONENTS, ENUMUNDEFINEDSOLUTIONCOMPONENTS, ENUMFORBIDDENSOLUTIONCOMPONENTS,
                        RESETENUMPRESENTSOLUTIONCOMPONENTS, RESETENUMUNDEFINEDSOLUTIONCOMPONENTS, RESETENUMFORBIDDENSOLUTIONCOMPONENTS,
                        ENUMADDMOVE, ENUMREMOVEMOVE, ENUMFORBIDMOVE, ENUMPERMITMOVE, RESETENUMADDMOVE, RESETENUMREMOVEMOVE, RESETENUMFORBIDMOVE,
                        RESETENUMPERMITMOVE, HEURISTICADDMOVE, HEURISTICREMOVEMOVE, HEURISTICFORBIDMOVE, HEURISTICPERMITMOVE, HEURISTICADDMOVEWOR,
                        HEURISTICREMOVEMOVEWOR, HEURISTICFORBIDMOVEWOR, HEURISTICPERMITMOVEWOR, RESETHEURISTICADDMOVEWOR, RESETHEURISTICREMOVEMOVEWOR,
                        RESETHEURISTICFORBIDMOVEWOR, RESETHEURISTICPERMITMOVEWOR, RANDOMADDMOVE, RANDOMREMOVEMOVE, RANDOMFORBIDMOVE, RANDOMPERMITMOVE,
                        RANDOMADDMOVEWOR, RANDOMREMOVEMOVEWOR, RANDOMFORBIDMOVEWOR, RANDOMPERMITMOVEWOR, RESETRANDOMADDMOVEWOR, RESETRANDOMREMOVEMOVEWOR,
                        RESETRANDOMFORBIDMOVEWOR, RESETRANDOMPERMITMOVEWOR,
                        COPYSOLUTION, APPLYADDMOVE, APPLYREMOVEMOVE, APPLYFORBIDMOVE, APPLYPERMITMOVE,
                        COPYMOVE,
                        GETCOMPONENTFROMMOVE, GETADDHEURISTICVALUE, GETREMOVEHEURISTICVALUE, GETFORBIDHEURISTICVALUE, GETPERMITHEURISTICVALUE, GETADDOBJECTIVELBINCREMENT,
                        GETREMOVEOBJECTIVELBINCREMENT, GETFORBIDOBJECTIVELBINCREMENT, GETPERMITOBJECTIVELBINCREMENT};

static char *pfname[] = {"getNumComponents()", "getMaxSolutionSize()", "getMaxAddNeighbourhoodSize()", "getMaxRemoveNeighbourhoodSize()", "getMaxForbidNeighbourhoodSize()",
                         "getMaxPermitNeighbourhoodSize()",
                         "allocSolution()", "allocMove()", "freeSolution()", "freeMove()",
                         "printProblem()", "printSolution()", "printMove()",
                         "emptySolution()", "heuristicSolution()",
                         "getObjectiveVector()", "getObjectiveLB()", "getNumPresentSolutionComponents()", "getNumUndefinedSolutionComponents()", "getNumForbiddenSolutionComponents()",
                         "isFeasible()", "isComplete()", "isDominant()", "checkDominance()", "getAddNeighbourhoodSize()", "getRemoveNeighbourhoodSize()", "getForbidNeighbourhoodSize()",
                         "getPermitNeighbourhoodSize()", "enumPresentSolutionComponents()", "enumUndefinedSolutionComponents()", "enumForbiddenSolutionComponents()",
                         "resetEnumPresentSolutionComponents()", "resetEnumUndefinedSolutionComponents()", "resetEnumForbiddenSolutionComponents()",
                         "enumAddMove()", "enumRemoveMove()", "enumForbidMove()", "enumPermitMove()", "resetEnumAddMove()", "resetEnumRemoveMove()", "resetEnumForbidMove()",
                         "resetEnumPermitMove()", "heuristicAddMove()", "heuristicRemoveMove()", "heuristicForbidMove()", "heuristicPermitMove()", "heuristicAddMoveWOR()",
                         "heuristicRemoveMoveWOR()", "heuristicForbidMoveWOR()", "heuristicPermitMoveWOR()", "resetHeuristicAddMoveWOR()", "resetHeuristicRemoveMoveWOR()",
                         "resetHeuristicForbidMoveWOR()", "resetHeuristicPermitMoveWOR()", "randomAddMove()", "randomRemoveMove()", "randomForbidMove()", "randomPermitMove()",
                         "randomAddMoveWOR()", "randomRemoveMoveWOR()", "randomForbidMoveWOR()", "randomPermitMoveWOR()", "resetRandomAddMoveWOR()", "resetRandomRemoveMoveWOR()",
                         "resetRandomForbidMoveWOR()", "resetRandomPermitMoveWOR()",
                         "copySolution()", "applyAddMove()", "applyRemoveMove()", "applyForbidMove()", "applyPermitMove()",
                         "copyMove()",
                         "getComponentFromMove()", "getAddHeuristicValue()", "getRemoveHeuristicValue()", "getForbidHeuristicValue()", "getPermitHeuristicValue()", "getAddObjectiveLBIncrement()",
                         "getRemoveObjectiveLBIncrement()", "getForbidObjectiveLBIncrement()", "getPermitObjectiveLBIncrement()"};

/*********************************/
/* ----- Utility functions ----- */
/*********************************/

/*
 * Return true if all problem functions of a given property-based test are
 * implemented or false if otherwise
 */
static int is_implemented(const enum ProblemFunction *pf, const int pfn)
{
    int implemented = 1, i;
    for (i = 0; i < pfn; ++i)
        if (!pfmacro[pf[i]]) {
            printf(INFO "\t-> Skipping tests using %s function! <-\n", pfname[pf[i]]);
            implemented = 0;
        }
    return implemented;
}

static int expand(struct solution *s1, struct solution *s2, struct move **v, int i)
{
    int passed = 1;
    double obj1, obj2, obji;

    if (randomMove(v[i], s1, ADD) != NULL) {
        getObjectiveLB(&obj1, s1);
        applyMove(s1, v[i], ADD);
        getObjectiveLBIncrement(&obji, v[i], s2, ADD);
        applyMove(s2, v[i], ADD);
        getObjectiveLB(&obj2, s2);
        if (obj1 + obji != obj2) {
            printf("f(s1) = %.1f, f(s2) = %.1f\n\n", obj1 + obji, obj2);
            printSolution(s1);
            printSolution(s2);
            printMove(v[i]);
            passed = 0;
            goto RETURN;
        }

        if (!(passed = expand(s1, s2, v, i + 1)))
            goto RETURN;

        getObjectiveLB(&obj1, s1);
        applyMove(s1, v[i], REMOVE);
        getObjectiveLBIncrement(&obji, v[i], s2, REMOVE);
        applyMove(s2, v[i], REMOVE);
        getObjectiveLB(&obj2, s2);
        if (obj1 + obji != obj2) {
            printf("f(s1) = %.1f, f(s2) = %.1f\n\n", obj1 + obji, obj2);
            printSolution(s1);
            printSolution(s2);
            printMove(v[i]);
            passed = 0;
            goto RETURN;
        }

        getObjectiveLB(&obj1, s1);
        applyMove(s1, v[i], FORBID);
        getObjectiveLBIncrement(&obji, v[i], s2, FORBID);
        applyMove(s2, v[i], FORBID);
        getObjectiveLB(&obj2, s2);
        if (obj1 + obji != obj2) {
            printf("f(s1) = %.1f, f(s2) = %.1f\n\n", obj1 + obji, obj2);
            printSolution(s1);
            printSolution(s2);
            printMove(v[i]);
            passed = 0;
            goto RETURN;
        }

        if (!(passed = expand(s1, s2, v, i + 1)))
            goto RETURN;

        getObjectiveLB(&obj1, s1);
        applyMove(s1, v[i], PERMIT);
        getObjectiveLBIncrement(&obji, v[i], s2, PERMIT);
        applyMove(s2, v[i], PERMIT);
        getObjectiveLB(&obj2, s2);
        if (obj1 + obji != obj2) {
            printf("f(s1) = %.1f, f(s2) = %.1f\n\n", obj1 + obji, obj2);
            printSolution(s1);
            printSolution(s2);
            printMove(v[i]);
            passed = 0;
            goto RETURN;
        }
    }
RETURN:
    return passed;
}

/**********************************/
/* ------- Property tests ------- */
/**********************************/

static int t_incremental_evaluation(struct problem *p)
{
#if GETNUMCOMPONENTS && ALLOCSOLUTION && ALLOCMOVE && FREESOLUTION && FREEMOVE && PRINTSOLUTION && PRINTMOVE && EMPTYSOLUTION && GETOBJECTIVELB && RANDOMADDMOVE && APPLYADDMOVE && APPLYREMOVEMOVE && APPLYFORBIDMOVE && APPLYPERMITMOVE && GETADDOBJECTIVELBINCREMENT && GETREMOVEOBJECTIVELBINCREMENT && GETFORBIDOBJECTIVELBINCREMENT && GETPERMITOBJECTIVELBINCREMENT
    struct solution *s1, *s2;
    struct move **v;
    int m, passed, i;

    /* Memory allocation */
    m = getNumComponents(p);
    s1 = allocSolution(p);
    s2 = allocSolution(p);
    v = malloc((m + 1) * sizeof(struct move *));
    for (i = 0; i <= m; ++i)
        v[i] = allocMove(p);

    /* Run */
    emptySolution(s1);
    emptySolution(s2);
    passed = expand(s1, s2, v, 0);

    /* Clean up */
    freeSolution(s1);
    freeSolution(s2);
    for (i = 0; i <= m; ++i)
        freeMove(v[i]);
    free(v);

    return passed;
#else
    return 0;
#endif
}

static int t_incremental_evaluation_add_remove(struct problem *p)
{
#if ALLOCSOLUTION && ALLOCMOVE && FREESOLUTION && FREEMOVE && PRINTSOLUTION && PRINTMOVE && EMPTYSOLUTION && GETOBJECTIVELB && RANDOMADDMOVE && RANDOMREMOVEMOVE && APPLYADDMOVE && APPLYREMOVEMOVE && GETADDOBJECTIVELBINCREMENT && GETREMOVEOBJECTIVELBINCREMENT
    struct solution *s1, *s2;
    struct move *v;
    int passed = 1, i;
    double obj1, obj2, obji;

    /* Memory allocation */
    s1 = allocSolution(p);
    s2 = allocSolution(p);
    v = allocMove(p);

    /* Run */
    emptySolution(s1);
    emptySolution(s2);
    for (i = 0; i < MAX_MOVE; ++i) {
        if (gsl_rng_uniform_int(rng, 2)) {
            if (randomMove(v, s1, ADD) != NULL) {
                getObjectiveLB(&obj1, s1);
                applyMove(s1, v, ADD);
                getObjectiveLBIncrement(&obji, v, s2, ADD);
                applyMove(s2, v, ADD);
                getObjectiveLB(&obj2, s2);
                if (obj1 + obji != obj2) {
                    printf("f(s1) = %.1f, f(s2) = %.1f\n\n", obj1 + obji, obj2);
                    printSolution(s1);
                    printSolution(s2);
                    printMove(v);
                    passed = 0;
                    break;
                }
            }
        }
        else {
            if (randomMove(v, s1, REMOVE) != NULL) {
                getObjectiveLB(&obj1, s1);
                applyMove(s1, v, REMOVE);
                getObjectiveLBIncrement(&obji, v, s2, REMOVE);
                applyMove(s2, v, REMOVE);
                getObjectiveLB(&obj2, s2);
                if (obj1 + obji != obj2) {
                    printf("f(s1) = %.1f, f(s2) = %.1f\n\n", obj1 + obji, obj2);
                    printSolution(s1);
                    printSolution(s2);
                    printMove(v);
                    passed = 0;
                    break;
                }
            }
        }
    }

    /* Clean up */
    freeSolution(s1);
    freeSolution(s2);
    freeMove(v);

    return passed;
#else
    return 0;
#endif
}

static int t_incremental_evaluation_add(struct problem *p)
{
#if ALLOCSOLUTION && ALLOCMOVE && FREESOLUTION && FREEMOVE && PRINTSOLUTION && PRINTMOVE && EMPTYSOLUTION && GETOBJECTIVELB && RANDOMADDMOVE && APPLYADDMOVE && GETADDOBJECTIVELBINCREMENT
    struct solution *s1, *s2;
    struct move *v;
    int passed = 1;
    double obj1, obj2, obji;

    /* Memory allocation */
    s1 = allocSolution(p);
    s2 = allocSolution(p);
    v = allocMove(p);

    /* Run */
    emptySolution(s1);
    emptySolution(s2);
    while (randomMove(v, s1, ADD) != NULL) {
        getObjectiveLB(&obj1, s1);
        applyMove(s1, v, ADD);
        getObjectiveLBIncrement(&obji, v, s2, ADD);
        applyMove(s2, v, ADD);
        getObjectiveLB(&obj2, s2);
        if (obj1 + obji != obj2) {
            printf("f(s1) = %.1f, f(s2) = %.1f\n\n", obj1 + obji, obj2);
            printSolution(s1);
            printSolution(s2);
            printMove(v);
            passed = 0;
            break;
        }
    }

    /* Clean up */
    freeSolution(s1);
    freeSolution(s2);
    freeMove(v);

    return passed;
#else
    return 0;
#endif
}

static struct property_test test[] = {
    {.name = "Test incremental lower bound evaluation (ADD and REMOVE only)",
     .description = "",
     .function = &t_incremental_evaluation_add_remove,
     .pf = (enum ProblemFunction[14]){_allocSolution_, _allocMove_, _freeSolution_, _freeMove_, _printSolution_, _printMove_, _emptySolution_, _getObjectiveLB_, _randomAddMove_, _randomRemoveMove_, _applyAddMove_, _applyRemoveMove_, _getAddObjectiveLBIncrement_, _getRemoveObjectiveLBIncrement_},
     .pfn = 14},
    {.name = "Test incremental lower bound evaluation (ADD only)",
     .description = "",
     .function = &t_incremental_evaluation_add,
     .pf = (enum ProblemFunction[11]){_allocSolution_, _allocMove_, _freeSolution_, _freeMove_, _printSolution_, _printMove_, _emptySolution_, _getObjectiveLB_, _randomAddMove_, _applyAddMove_, _getAddObjectiveLBIncrement_},
     .pfn = 11}};

/**********************************/
/* -------- Main program -------- */
/**********************************/

/*
 * 
 * Status: CHECK
 */
int runPropertyBased(struct problem *p, int rep)
{
    struct property_test t;
    int ntest = sizeof(test) / sizeof(struct property_test), skipped = 0, passedbytype = 0, result, passed, i, j, k;

    printf("----\nWelcome to property-based testing!\n----\n");
    printf("Number of tests: %d\n", ntest);
    for (i = 0; i < ntest; ++i) {
        t = test[i]; /* property-based test */

        /* Print property-based test information */
        printf("-----------------------------\n");
        printf(BOLD "Test %d: %s\n" RESET, i + 1, t.name);
        printf("\t%s\n", t.description);
        printf("----\n");
        printf("Requires:\n");
        for (j = 0; j < t.pfn; ++j)
            printf("\t%s\n", pfname[t.pf[j]]);

        /* Run property-based test if it is implemented */
        if (is_implemented(t.pf, t.pfn)) {
            for (k = 0, passed = 0; k < rep; ++k) {
                printf("----\n" INFO);
                printf("Iteration %d\n", k + 1);
                result = (*t.function)(p);
                if (result)
                    printf(OK "\tPassed!\n" RESET);
                else
                    printf(NOTOK "\tFailed!\n" RESET);
                passed += result;
            }
            if (passed == rep)
                passedbytype++;
        }
        else {
            printf("----\n");
            printf(WARN "SKIPPED!\n " RESET);
            skipped++;
        }
        printf("-----------------------------\n");
    }

    /* Print summary */
    printf("----------------------------- \n");
    printf("Total number of tests available: %d\n", ntest);
    if (passedbytype == ntest - skipped)
        printf(OK "PASSED ALL TESTS (%d)!\n" RESET, ntest - skipped);
    else
        printf(NOTOK "FAILED %d/%d TESTS!\n" RESET, ntest - passedbytype, ntest);
    if (skipped > 0)
        printf(WARN "%d TEST(S) SKIPPED!\n" RESET, skipped);
    printf("-----------------------------\n");

    return 0;
}
