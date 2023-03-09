/*
 * File: observerc.c
 * Description: C-routines to monitor and save data.
 *
 *
 *
 */

#include <stdio.h>
#include <stdint.h>

void print_matf(double *mat, int32_t ni, int32_t nj)
 {
   int32_t i, j;
   printf("ni: %d, nj: %d\n", ni, nj);
   for (i=0; i<nj; i++)
   {
     printf("f=====");
   }
   printf("\n");

   for (i=0; i<ni; i++)
   {
     for (j=0; j<nj; j++)
     {
       printf("%6.2f", *(mat + j*ni + i));
     }
     printf("\n");
   }
 }

void print_mat(double **mat, int32_t ni, int32_t nj)
 {
   int i, j;
   for (i=0; i<nj; i++)
   {
     printf("======");
   }
   printf("\n");

   for (i=0; i<ni; i++)
   {
     for (j=0; j<nj; j++)
     {
       printf("%6.2f", mat[i][j]);
     }
     printf("\n");
   }
 }

