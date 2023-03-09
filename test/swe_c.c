/*
 * File: swe_c.c 
 * Description: Program to test libswe functions in C99
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "equationsc.h"
#include "observerc.h"

/*
void print_mat(double **mat, int ni, int nj)
{
  int i, j;
  for (i=0; i<ni; i++)
  {
    for (j=0; j<nj; j++)
    {
      //printf("%6.2f", *(*(mat+j*ni)+i));
      printf("%6.2f", mat[i][j]);
    }
    printf("\n");
  }
}
*/

int main(int argc, char* argv[])
{
  int status = 0;
  int i, j, count;
  const int ni=10, nj=5;
  //double arr2d[ni][nj];
  double **arr2d;

  /* allocate array */
  
  
  count = 0;
  // allocate memory
  arr2d = malloc(ni*sizeof(double));
  for (i=0; i<ni; i++) {
      //arr2d[i][j] = (double)count;
      arr2d[i] = malloc(nj*sizeof(double));
  }

  // initialize matrix
  for (i=0; i<ni; i++)
  {
    for (j=0; j<nj; j++)
    {
      //*(*(arr2d+j)+i) = (double)count;
      arr2d[i][j] = (double)count;
      count++;
    }        
  }

 print_mat(arr2d, ni, nj);
 status = twice_arr_2d(arr2d, ni, nj);
 //status = twice_arr_2d(ni, nj);
 print_mat(arr2d, ni, nj);

  // free memory
  for (i=0; i<ni; i++) {
    free(arr2d[i]);
  }
  free(arr2d);

  printf("No problem detected.\n");


  return status;
}
