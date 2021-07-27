/* minmax_heap.c
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

#include "minmax_heap.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define is_min(n) (((int)log2(n) & 1) == 0)
#define parent(n) (n / 2)
#define left_child(n) (n * 2)
#define right_child(n) (n * 2 + 1)

struct heap {
    int *key;
    double *val;
    int *bit;
    int count;
    int size;
};

/*********************************/
/* ----- Utility functions ----- */
/*********************************/

static void swap(struct heap *h, int i, int j)
{
    if (i == j)
        return;
    /* swap keys */
    int key = h->key[i];
    h->key[i] = h->key[j];
    h->key[j] = key;
    /* swap values */
    double val = h->val[i];
    h->val[i] = h->val[j];
    h->val[j] = val;
}

static void bubbleUpMin(struct heap *h, int i)
{
    int pp_idx = parent(parent(i));
    if (pp_idx <= 0) /* no grandparent */
        return;
    if (h->val[i] < h->val[pp_idx]) { /* move up */
        swap(h, i, pp_idx);
        bubbleUpMin(h, pp_idx);
    }
}

static void bubbleUpMax(struct heap *h, int i)
{
    int pp_idx = parent(parent(i));
    if (pp_idx <= 0) /* no grandparent */
        return;
    if (h->val[i] > h->val[pp_idx]) { /* move down */
        swap(h, i, pp_idx);
        bubbleUpMax(h, pp_idx);
    }
}

static void bubbleUp(struct heap *h, int i)
{
    int p_idx = parent(i);
    if (p_idx <= 0) /* no parent */
        return;
    if (is_min(i)) { /* in min-level */
        if (h->val[i] > h->val[p_idx]) { /* move down on max-levels */
            swap(h, i, p_idx);
            bubbleUpMax(h, p_idx);
        }
        else /* move up on min-levels */
            bubbleUpMin(h, i);
    }
    else { /* in max-level */
        if (h->val[i] < h->val[p_idx]) { /* move up on min-levels */
            swap(h, i, p_idx);
            bubbleUpMin(h, p_idx);
        }
        else /* move down on max-levels */
            bubbleUpMax(h, i);
    }
}

static int minChildOrGrandchild(struct heap *h, int i)
{
    int m = 0;
    int l = left_child(i);
    int r = right_child(i);
    int ll = left_child(l);
    int lr = right_child(l);
    int rl = left_child(r);
    int rr = right_child(r);
    if (l <= h->count)
        m = l;
    if (r <= h->count && h->val[r] < h->val[m])
        m = r;
    if (ll <= h->count && h->val[ll] < h->val[m])
        m = ll;
    if (lr <= h->count && h->val[lr] < h->val[m])
        m = lr;
    if (rl <= h->count && h->val[rl] < h->val[m])
        m = rl;
    if (rr <= h->count && h->val[rr] < h->val[m])
        m = rr;
    return m;
}

static int maxChildOrGrandchild(struct heap *h, int i)
{
    int m = 0;
    int l = left_child(i);
    int r = right_child(i);
    int ll = left_child(l);
    int lr = right_child(l);
    int rl = left_child(r);
    int rr = right_child(r);
    if (l <= h->count)
        m = l;
    if (r <= h->count && h->val[r] > h->val[m])
        m = r;
    if (ll <= h->count && h->val[ll] > h->val[m])
        m = ll;
    if (lr <= h->count && h->val[lr] > h->val[m])
        m = lr;
    if (rl <= h->count && h->val[rl] > h->val[m])
        m = rl;
    if (rr <= h->count && h->val[rr] > h->val[m])
        m = rr;
    return m;
}

static void trickleDownMin(struct heap *h, int i)
{
    int m = minChildOrGrandchild(h, i);
    if (m == 0)
        return;
    int p_idx = parent(m);
    if (m > right_child(i)) { /* m is grandchild */
        if (h->val[m] < h->val[i]) {
            swap(h, i, m);
            if (h->val[m] > h->val[p_idx])
                swap(h, m, p_idx);
            trickleDownMin(h, m);
        }
    }
    else { /* m is child */
        if (h->val[m] < h->val[i])
            swap(h, i, m);
    }
}

static void trickleDownMax(struct heap *h, int i)
{
    int m = maxChildOrGrandchild(h, i);
    if (m == 0)
        return;
    int p_idx = parent(m);
    if (m > right_child(i)) { /* m is grandchild */
        if (h->val[m] > h->val[i]) {
            swap(h, i, m);
            if (h->val[m] < h->val[p_idx])
                swap(h, m, p_idx);
            trickleDownMax(h, m);
        }
    }
    else { /* m is child */
        if (h->val[m] > h->val[i])
            swap(h, i, m);
    }
}

static void trickleDown(struct heap *h, int i)
{
    if (is_min(i)) /* in min-level */
        trickleDownMin(h, i);
    else /* in max-level */
        trickleDownMax(h, i);
}

/**********************/
/* Heap instantiation */
/**********************/

/*
 * Heap structure instantiation
 */
struct heap *newHeap(int size)
{
    int i;
    struct heap *h = malloc(sizeof(struct heap));
    h->key = malloc((size + 1) * sizeof(int));
    h->val = malloc((size + 1) * sizeof(double));
    h->bit = malloc(size * sizeof(int));
    for (i = 0; i < size; ++i)
        h->bit[i] = 0;
    h->count = 0;
    h->size = size;
    return h;
}

/*******************/
/* Heap inspection */
/*******************/

/*
 * Return true if a given heap structure is full or false if it is not full
 */
int isHeapFull(struct heap *h)
{
    return h->count == h->size;
}

/*********************/
/* Memory management */
/*********************/

/*
 * Free the memory used by a heap structure
 */
void freeHeap(struct heap *h)
{
    free(h->key);
    free(h->val);
    free(h->bit);
    free(h);
}

/*******/
/* I/O */
/*******/

/*
 * Prints the heap structure data.
 */
void printHeap(struct heap *h)
{
    int i, j;
    double pow_2, log_2 = log2(h->count);
    printf("Number of elements is %d in heap structure of size %d\n", h->count, h->size);
    for (i = 0; i < log_2; ++i) {
        pow_2 = pow(2, i);
        for (j = 0; j < pow_2 && j + pow_2 <= h->count; ++j)
            printf("%d(%.1lf) ", h->key[j + (int)pow_2], h->val[j + (int)pow_2]);
        printf("\n");
    }
    printf("\n");
}

/*******************/
/* Heap generation */
/*******************/

/*
 * Constructs a heap structure from a given array of values
 */
struct heap *buildHeap(struct heap *h, double *val, int size)
{
    if (size != h->size) {
        fprintf(stderr, "buildHeap() error: Invalid size, incompatible with heap size.\n");
        return NULL;
    }
    int i;
    for (i = 0; i < size; ++i) {
        h->key[i + 1] = i;
        h->bit[i] = 1;
    }
    memcpy(h->val + 1, val, size * sizeof(double));
    h->count = size;
    /* Floyd's linear-time heap construction */
    for (i = size / 2; i > 0; --i)
        trickleDown(h, i);
    return h;
}

/***********************/
/* Operations on heaps */
/***********************/

/*
 * Inserts a given element to a heap structure
 */
struct heap *insertHeap(struct heap *h, int key, double val)
{
    if (key < 0 || key >= h->size || h->bit[key]) { /* invalid key */
        fprintf(stderr, "insertHeap() error: Invalid key.\n");
        return NULL;
    }
    h->key[++h->count] = key;
    h->val[h->count] = val;
    h->bit[key] = 1;
    bubbleUp(h, h->count);
    return h;
}

/*
 * Returns the key of the minimum element of a given heap structure, and removes
 * it from the heap
 */
int deleteHeapMin(struct heap *h)
{
    int key;
    if (h->count > 2) {
        key = h->key[1];
        h->key[1] = h->key[h->count];
        h->val[1] = h->val[h->count--];
        h->bit[key] = 0;
        trickleDown(h, 1);
        return key;
    }
    else if (h->count == 2) {
        key = h->key[1];
        h->key[1] = h->key[h->count];
        h->val[1] = h->val[h->count--];
        h->bit[key] = 0;
        return key;
    }
    else if (h->count == 1) {
        h->bit[h->key[1]] = 0;
        return h->key[h->count--];
    }
    return -1;
}

/*
 * Returns the key of the maximum element of a given heap structure, and removes
 * it from the heap
 */
int deleteHeapMax(struct heap *h)
{
    int i, key;
    if (h->count > 3) {
        i = h->val[2] > h->val[3] ? 2 : 3;
        key = h->key[i];
        h->key[i] = h->key[h->count];
        h->val[i] = h->val[h->count--];
        h->bit[key] = 0;
        trickleDown(h, i);
        return key;
    }
    else if (h->count == 3) {
        i = h->val[2] > h->val[3] ? 2 : 3;
        key = h->key[i];
        h->key[i] = h->key[h->count];
        h->val[i] = h->val[h->count--];
        h->bit[key] = 0;
        return key;
    }
    else if (h->count > 0) {
        h->bit[h->key[h->count]] = 0;
        return h->key[h->count--];
    }
    return -1;
}

/*
 * Returns the value of the minimum element of a given heap structure
 */
double findHeapMin(struct heap *h)
{
    if (h->count > 0)
        return h->val[1];
    return DBL_MAX;
}

/*
 * Returns the value of the maximum element of a given heap structure
 */
double findHeapMax(struct heap *h)
{
    if (h->count > 2)
        return h->val[2] > h->val[3] ? h->val[2] : h->val[3];
    else if (h->count == 2)
        return h->val[2];
    else if (h->count == 1)
        return h->val[1];
    return -DBL_MAX;
}
