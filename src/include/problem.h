/* problem.h
 *
 * (C) 2021 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
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

/* This header file contains all problem independent definitions */

/****************/
/* Enumerations */
/****************/

enum ComponentState { PRESENT = 1, FORBIDDEN = 2, UNDEFINED = 4 };
enum SubNeighbourhood { ADD = 1, REMOVE = 2, FORBID = 4, PERMIT = 8 };

/*******************/
/* Data structures */
/*******************/

struct problem;
struct solution;
struct move;

/*************************/
/* Problem instantiation */
/*************************/

/* newProblem() allocates a problem structure and initialises it. Being
 * problem-specific, the function arguments are deliberately left unspecified.
 * The function returns a pointer to the problem if instantiation succeeds or
 * NULL if it fails.
 */
struct problem *newProblem();

/**********************/
/* Problem inspection */
/**********************/

/* getNumComponents() returns the size of the ground set of a problem instance.
 */
long getNumComponents(const struct problem *p);

/* getMaxSolutionSize() returns the largest number of components that a solution
 * can potentially have.
 */
long getMaxSolutionSize(const struct problem *p);

/* getMaxNeighbourhoodSize() returns the largest possible number of direct
 * neighbours in a given subneighbourhood (or a larger number).
 */
long getMaxNeighbourhoodSize(const struct problem *p, const enum SubNeighbourhood nh);

/*********************/
/* Memory management */
/*********************/

/* allocSolution() allocates memory for a solution structure.
 * The function returns a pointer to the solution if allocation succeeds or NULL
 * if it fails.
 */
struct solution *allocSolution(struct problem *p);

/* allocMove() allocates memory for a move structure.
 * The function returns a pointer to the move if allocation succeeds or NULL if
 * it fails.
 */
struct move *allocMove(struct problem *p);

/* freeProblem() deallocates the memory used by a problem structure.
 */
void freeProblem(struct problem *p);

/* freeSolution() deallocates the memory used by a solution structure.
 */
void freeSolution(struct solution *s);

/* freeMove() deallocates the memory used by a move structure.
 */
void freeMove(struct move *v);

/*************/
/* Reporting */
/*************/

/* printProblem() prints a user-formatted representation of a problem instance.
 */
void printProblem(const struct problem *p);

/* printSolution() prints a user-formatted representation of a solution.
 */
void printSolution(const struct solution *s);

/* printMove() prints a user-formatted representation of a move.
 */
void printMove(const struct move *v);

/***********************/
/* Solution generation */
/***********************/

/* emptySolution() initialises a solution structure as an empty solution.
 * The input argument must be a pointer to a solution previously allocated with
 * allocSolution(), which is modified in place.
 * If any of the following functions:
 *     - enumSolutionComponents()
 *     - enumMove()
 *     - heuristicMoveWOR()
 *     - randomMoveWOR()
 * is implemented, emptySolution() must also reset the corresponding state by
 * performing the equivalent to:
 *     - resetEnumSolutionComponents()
 *     - resetEnumMove()
 *     - resetHeuristicMoveWOR()
 *     - resetRandomMoveWOR()
 * respectively.
 * The function returns its input argument.
 */
struct solution *emptySolution(struct solution *s);

/* heuristicSolution() heuristically constructs a feasible solution, preserving
 * all present and forbidden components in a given solution, which is modified
 * in place.
 * If any of the following functions:
 *     - enumSolutionComponents()
 *     - enumMove()
 *     - heuristicMoveWOR()
 *     - randomMoveWOR()
 * is implemented, heuristicSolution() must also reset the corresponding state
 * by performing the equivalent to:
 *     - resetEnumSolutionComponents()
 *     - resetEnumMove()
 *     - resetHeuristicMoveWOR()
 *     - resetRandomMoveWOR()
 * respectively.
 * The function returns its input argument if a new feasible solution is
 * generated or NULL if no feasible solution is found, in which case the
 * original input argument is lost.
 */
struct solution *heuristicSolution(struct solution *s);

/***********************/
/* Solution inspection */
/***********************/

/* getObjectiveVector() supports single or multiple objective full and/or
 * incremental solution evaluation. This function modifies the array passed as
 * first input argument to contain objective values. Once a solution is
 * evaluated, results may be cached in the solution itself so that a subsequent
 * call to this function simply returns the precomputed value(s) and/or the
 * solution can be re-evaluated more efficiently after it is modified by one or
 * more calls to applyMove(). Therefore, the formal argument is not const.
 * The function returns its first input argument if a given solution is feasible
 * or NULL if it is unfeasible, in which case the first input argument is left
 * unspecified (in particular, it may have been modified).
 */
double *getObjectiveVector(double *objv, struct solution *s);

/* getObjectiveLB() supports single or multiple objective full and/or
 * incremental lower bound evaluation. This function modifies the array passed
 * as first input argument to contain lower bound values. The lower bound of a
 * solution must be less than or equal to the lower bound of another solution if
 * the sets of present and forbidden components of the former are both subsets
 * of the corresponding sets of components of the latter. Once the lower bound
 * of a solution is evaluated, results may be cached in the solution itself so
 * that a subsequent call to this function simply returns the precomputed
 * value(s) and/or the solution can be re-evaluated more efficiently after it is
 * modified by one or more calls to applyMove(). Therefore, the formal argument
 * is not const.
 * The function returns its first input argument.
 */
double *getObjectiveLB(double *objLB, struct solution *s);

/* getNumSolutionComponents() returns the number of components of a solution
 * that are in a given state.
 */
long getNumSolutionComponents(const struct solution *s, const enum ComponentState st);

/* isFeasible() returns true if a given solution is feasible or false if it is
 * unfeasible.
 * The result may be cached in the solution in order to speed up future
 * operations. Therefore, the input argument is not const.
 */
int isFeasible(struct solution *s);

/* isComplete() returns true if no components can be added to a given feasible
 * solution or false if there are components that can be added to it.
 * The result may be cached in the solution in order to speed up future
 * operations. Therefore, the input argument is not const.
 */
int isComplete(struct solution *s);

/* getNeighbourhoodSize() returns the number of direct neighbours in a given
 * subneighbourhood of a solution.
 * The result may be cached in the solution in order to speed up future
 * operations. Therefore, the first input argument is not const.
 */
long getNeighbourhoodSize(struct solution *s, const enum SubNeighbourhood nh);

/* enumSolutionComponents() implements the enumeration of the components of a
 * solution that are in a given state, in unspecified order.
 * The function returns a unique component identifier in the range 0..|G|-1 if a
 * new component is enumerated or a negative integer if there are no components
 * left.
 */
long enumSolutionComponents(struct solution *s, const enum ComponentState st);

/* resetEnumSolutionComponents() resets the enumeration of the components of a
 * solution that are in a given state, so that the next call to
 * enumSolutionComponents() will start the enumeration of the components in that
 * state from the beginning. The function returns its first input argument.
 */
struct solution *resetEnumSolutionComponents(struct solution *s, const enum ComponentState st);

/*******************/
/* Move generation */
/*******************/

/* enumMove() implements enumeration of a given subneighbourhood of a solution,
 * in an unspecified order. This is intended to support fast neighbourhood
 * exploration and evaluation with getObjectiveLBIncrement(), particularly when
 * a large part of the neighbourhood is to be visited.
 * The first input argument must be a pointer to a move previously allocated
 * with allocMove(), which is modified in place. The function returns this
 * pointer if a new move is generated or NULL if there are no moves left.
 */
struct move *enumMove(struct move *v, struct solution *s, const enum SubNeighbourhood nh);

/* resetEnumMove() resets the enumeration of a given subneighbourhood of a
 * solution, so that the next call to enumMove() will start the enumeration of
 * that subneighbourhood from the beginning. The function returns its input
 * argument.
 */
struct solution *resetEnumMove(struct solution *s, const enum SubNeighbourhood nh);

/* heuristicMove() generates a heuristic move corresponding to a given
 * subneighbourhood of a solution. Calling heuristicMove() with the same
 * arguments multiple times may generate different moves. Heuristic moves may or
 * may not be greedy with respect to how they affect the lower bound.
 * The first input argument must be a pointer to a move previously allocated
 * with allocMove(), which is modified in place. The function returns this
 * pointer if a new move is generated or NULL if the subneighbourhood is empty.
 */
struct move *heuristicMove(struct move *v, struct solution *s, const enum SubNeighbourhood nh);

/* heuristicMoveWOR() implements heuristic sampling of a given subneighbourhood
 * of a solution, without replacement. Heuristic moves may or may not be
 * generated in ascending order of the corresponding lower bound increment.
 * The first input argument must be a pointer to a move previously allocated
 * with allocMove(), which is modified in place. The function returns this
 * pointer if a new move is generated or NULL if there are no moves left.
 */
struct move *heuristicMoveWOR(struct move *v, struct solution *s, const enum SubNeighbourhood nh);

/* resetHeuristicMoveWOR() resets the heuristic sampling without replacement of
 * a given subneighbourhood of a solution, so that the next call to
 * heuristicMoveWOR() will start the heuristic sampling from the beginning. The
 * order of enumeration may not be preserved by such a reset. The function
 * returns its input argument.
 */
struct solution *resetHeuristicMoveWOR(struct solution *s, const enum SubNeighbourhood nh);

/* randomMove() implements uniform random sampling with replacement of a given
 * subneighbourhood of a solution. This may be an empty set for some
 * subneighbourhoods, in which case the function returns NULL.
 * The first input argument must be a pointer to a move previously allocated
 * with allocMove(), which is modified in place.
 */
struct move *randomMove(struct move *v, struct solution *s, const enum SubNeighbourhood nh);

/* randomMoveWOR() implements uniform random sampling of a given
 * subneighbourhood of a solution, without replacement.
 * The first input argument must be a pointer to a move previously allocated
 * with allocMove(), which is modified in place. The function returns this
 * pointer if a new move is generated or NULL if there are no moves left.
 */
struct move *randomMoveWOR(struct move *v, struct solution *s, const enum SubNeighbourhood nh);

/* resetRandomMoveWOR() resets the uniform random sampling without replacement
 * of a given subneighbourhood of a solution, so that any move corresponding to
 * that subneighbourhood can be generated by the next call to randomMoveWOR().
 * The function returns its input argument.
 */
struct solution *resetRandomMoveWOR(struct solution *s, const enum SubNeighbourhood nh);

/***************************/
/* Operations on solutions */
/***************************/

/* copySolution() copies the contents of the second argument to the first
 * argument, which must have been previously allocated with allocSolution(). The
 * copied solution is functionally indistinguishable from the original solution.
 * The function returns its first input argument.
 */
struct solution *copySolution(struct solution *dest, const struct solution *src);

/* applyMove() modifies a solution in place by applying a move to it. It is
 * assumed that the move was generated for, and possibly evaluated with respect
 * to, that particular solution and the given subneighbourhood. In addition,
 * once a move is applied to a solution, it can be reverted by applying it again
 * with the opposite subneighbourhood. For example, after an ADD move generated
 * for a given solution is applied to that solution, it may be applied again as
 * a REMOVE move to the resulting solution in order to recover the original
 * solution. The result of applying a move to a solution in any other
 * circumstances is undefined.
 * If any of the following functions:
 *     - enumSolutionComponents()
 *     - enumMove()
 *     - heuristicMoveWOR()
 *     - randomMoveWOR()
 * is implemented, applyMove() must also reset the corresponding state
 * by performing the equivalent to:
 *     - resetEnumSolutionComponents()
 *     - resetEnumMove()
 *     - resetHeuristicMoveWOR()
 *     - resetRandomMoveWOR()
 * respectively, for all subneighbourhoods.
 * The function returns its first input argument.
 */
struct solution *applyMove(struct solution *s, const struct move *v, const enum SubNeighbourhood nh);

/***********************/
/* Operations on moves */
/***********************/

/* copyMove() copies the contents of the second argument to the first argument,
 * which must have been previously allocated with allocMove(). The copied move
 * is functionally indistinguishable from the original move.
 * The function returns its first input argument.
 */
struct move *copyMove(struct move *dest, const struct move *src);

/*******************/
/* Move inspection */
/*******************/

/* getComponentFromMove() returns the unique component identifier in the range
 * 0..|G|-1 with respect to a given move.
 */
long getComponentFromMove(const struct move *v);

/* getObjectiveLBIncrement() supports single or multiple objective move
 * evaluation with respect to the solution for which it was generated, before it
 * is actually applied to that solution (if it ever is). This function modifies
 * the array passed as first input argument to contain lower bound increment
 * values. The result of evaluating a move with respect to a solution other than
 * that for which it was generated (or to a pristine copy of it) is undefined.
 * Once a move is evaluated, results may be cached in the move itself, so that
 * they can be used by applyMove() to update the evaluation state of the
 * solution more efficiently. In addition, results may also be cached in the
 * solution in order to speed up evaluation of future moves. Consequently,
 * neither formal argument is const. The function returns its first input
 * argument.
 */
double *getObjectiveLBIncrement(double *obji, struct move *v, struct solution *s, const enum SubNeighbourhood nh);
