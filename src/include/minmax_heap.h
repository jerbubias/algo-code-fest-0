/* minmax_heap.h
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

/* Adapted from the min-max heap implementation by ilbanshee in
 * https://github.com/ilbanshee/min-max-heap */

/*******************/
/* Data structures */
/*******************/

struct heap;

/**********************/
/* Heap instantiation */
/**********************/

/* newHeap() allocates a heap structure and initialises it as an empty heap.
 */
struct heap *newHeap(int size);

/*******************/
/* Heap inspection */
/*******************/

/* isHeapFull() returns true if a given heap structure is full or false if it is
 * not full.
 */
int isHeapFull(struct heap *h);

/*********************/
/* Memory management */
/*********************/

/* freeHeap() deallocates a heap structure.
 */
void freeHeap(struct heap *h);

/*******/
/* I/O */
/*******/

/* printHeap() outputs the heap structure data.
 */
void printHeap(struct heap *h);

/*******************/
/* Heap generation */
/*******************/

/* buildHeap() constructs a heap structure from a given array of values. The
 * first input argument must be a pointer to a heap previously allocated with
 * newHeap(), which is modified in place. The function returns its first input
 * argument if the heap structure is successfully built or NULL if the array
 * size is incompatible with the heap size.
 */
struct heap *buildHeap(struct heap *h, double *val, int size);

/***********************/
/* Operations on heaps */
/***********************/

/* insertHeap() inserts a given element to a heap structure. The function
 * returns its first input argument if the element is successfully inserted or
 * NULL if the element key is not unique or in the range 0..N-1, where N denotes
 * the size of the heap structure.
 */
struct heap *insertHeap(struct heap *h, int key, double val);

/* deleteHeapMin() returns the key of the minimum element of a given heap
 * structure, and removes it from the heap.
 */
int deleteHeapMin(struct heap *h);

/* deleteHeapMin() returns the key of the maximum element of a given heap
 * structure, and removes it from the heap.
 */
int deleteHeapMax(struct heap *h);

/* findHeapMin() returns the value of the minimum element of a given heap
 * structure.
 */
double findHeapMin(struct heap *h);

/* findHeapMax() returns the value of the maximum element of a given heap
 * structure.
 */
double findHeapMax(struct heap *h);
