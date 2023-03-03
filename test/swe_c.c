/*
 * File: swe_c.c
 * Description: Test_routine to test swe equations.
 */

#include <stdint.h>
#include <stdio.h>
#include "equationsc.h"

int main(int argc, char* argv[])
{
  int32_t i, status;
  const int32_t narray = 10;
  double a[narray];

  for (i=0; i<narray; i++)
  {
    a[i] = (double)(i)*(double)(i);
  }

  for (i=0; i<narray; i++)
  {
    printf("%.2f ", a[i]);
  }
  printf("\n");

  status = twice_arr(a, narray);

  for (i=0; i<narray; i++)
  {
    printf("%.2f ", a[i]);
  }
  printf("\n");

  return status;
}
