/*
 * File: equations.c
 * Description: Collection of C routines that implement the functionality of 
 *              Shallow Water Equations (SWE)
 *
 */
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "equationsc.h"


/*
 * Function: say_hi
 * Description: Test function that returns int. To be deleted.
 * TODO: delete when finished testing
 */
int say_hi(void){
  printf("Hello from C.\n");
  return 0;
}

/*
 * Function: twice_arr
 * Description: test function to double the value of a 1d array 
 * TODO: delete when finished testing
 */
int twice_arr(double *arr, int32_t arr_size){
  int32_t i;

  for (i=0; i<arr_size; i++)
  {
    *arr = 2.0* (*arr);
    arr++;
  }

  return 0; // no error

}

/*
 * Function: avg_arr
 * Description: test function to average array with its nearest neighbors
 * TODO: delete when finished testing
 */
int avg_arr(double *arr, int32_t arr_size){
  int32_t i;
  double *arr_tmp = malloc(sizeof(double) * arr_size);  // temporary space for averaging
  double *arr_head, *arr_tmp_head;

  // backup array pointers
  arr_head = arr;
  arr_tmp_head = arr_tmp;
  
  for (i=0; i<arr_size; i++)
  {
    *arr_tmp = *arr;
    arr++;
    arr_tmp++;
  }

  // reset arr pointers
  arr = arr_head;
  arr_tmp = arr_tmp_head;

  // left-most value 
  *arr = (*(arr_tmp+1) + *(arr_tmp) + *(arr_tmp + arr_size-1))/3.0;

  // right-most value
  *(arr + arr_size -1)  = (*(arr_tmp + arr_size-2) + *(arr_tmp + arr_size - 1) + *(arr_tmp))/3.0;

  for (i=1; i<(arr_size-1); i++)
  {
    *(arr+i)= (*arr_tmp + *(arr_tmp+i) + *(arr_tmp-i))/3.0;
  }

  // cleanup 
  free(arr_tmp);

  return 0; // no error

}

/*
 * Function: twice_arr_2d
 * Description: double the values of a 2d array
 * 
 */
int twice_arr_2d(double *arr, int32_t ni, int32_t nj)
{
  int32_t i, j;
  for (i=0; i<ni; i++)
  {
    for (j=0; j<nj; j++)
    {
      *arr = 2.0 * (*arr);
      arr++;
    }      
  }
  return 0;
}

/*
 * Function: diff_center_x
 * Description: Returns the central difference in x-direction of a 2d array
 *              with periodic boundary conditions.
 *
 * Input: r -- 2d array of values
 *        dr -- 2d array, x-directed central difference 
 *        dx -- x resolution
 *        nx -- x dimension of r, dr
 *        ny -- y dimension of r, dr
 * Returns: none 
 *
 */

//void diff_center_x(double *r, double *dr, double dx, int32_t nx, int32_t ny ) {
//  int32_t i, j;
//
//   /* difference of left-most indices */
//  for (j=0; j<(ny); j++)
//  {
//    *(dr*j) = (0.5*dx) * (*(r*(j+1)) - *(r*(j+nx-1)));
//  }
//
//  /* difference of right-most indices */
//  for (j=0; j<(ny); j++)
//  {
//    *(dr*(j+nx-1)) = (0.5*dx) * (*(r*(j+nx-2)) - *(r*j));
//  }
//  
//  /* difference of central region */
//  for (i=1; i<(nx-1); i++)
//  {
//    for (j=0; j<ny; j++)
//    {
//      *((dr*j)+i) = (0.5*dx) * (*(r*((j+1)+i)) - *(r*((j-1)+i))); 
//    }
//  }
//}
