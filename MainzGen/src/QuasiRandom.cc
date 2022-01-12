#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "QuasiRandom.h"

SobolSequence::SobolSequence(const int d)
{
  const short initmdeg[] = {0, 
			1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6};
  const int initip[] = {0, 
			 0, 1, 1, 2, 1, 4, 2, 4, 7,11,13,14, 1,13,16,19,22,25};
  const int initiv[] = {0, 
			 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
			 3, 1, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
			 5, 7, 7, 3, 3, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
			 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
			 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
			 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1};
  
  dim = (d > MAXDIM ? MAXDIM : d);
  xSet = new double[dim];
  in = 0;

  for (int k=0; k<=MAXDIM; k++) mdeg[k] = initmdeg[k]; // just copy
  for (int k=0; k<=MAXDIM; k++) ip[k]   = lint(initip[k],0);
  for (int k=1; k<=MAXDIM; k++) ix[k]   = lint(0,0);
  for (int k=0; k<=MAXDIM*mdeg[MAXDIM]; k++) 
    iv[k] = lint(initiv[k],0);

  for (int j=1,k=0; j <= MAXBIT; j++,k += MAXDIM) iu[j]=&iv[k]; //for 2d access

  for (int k=1; k<=MAXDIM; k++) {
    for (int j=1; j<=mdeg[k]; j++)   // normalize first elements
      iu[j][k] <<= (MAXBIT-j);
    for (int j=mdeg[k]+1; j<=MAXBIT; j++) {   // recurrency relation for vi
      lint ipp = ip[k];
      lint i   = iu[j-mdeg[k]][k];
      i ^= i >> mdeg[k];
      for(int l = mdeg[k]-1; l>=1; l--) {
	if (ipp & lint(1,0)) i ^= iu[j-l][k];
	ipp >>= 1;
      }
      iu[j][k]=i;
    }
  }
}  

double 
SobolSequence::operator()()
{
  if (access >= dim) {
    std::cerr << "ERROR: Sequence contains " << dim
	      << " values!" << std::endl << std::flush;
    exit(-1);
  }
  return xSet[access++];
}

void 
SobolSequence::operator()(double x[])
{
  const double fac = 1.0 / (lint(1,0) << MAXBIT);
  lint im = in++;
  int j = 0 ; 
  for(j = 1; j <= MAXBIT; j++) if (im & lint(1,0)) im >>= 1; else break;
  if (j > MAXBIT) {
    std::cerr << "\aERROR: Max length of Sobol' Sequence reached: "
	      << in << std::endl;
    exit(-1);
  }
  for (int k=1; k<=dim; k++) {
    ix[k]  ^= iv[(j-1) * MAXDIM + k];
    x[k-1] =  ix[k] * fac;
  }
}

void
SobolSequence::init(double start)
{
  in = start;
  for (int k=1; k<=dim; k++) {
    lint graycode = lint(in) ^ (lint(in) >> 1);
    ix[k] = lint(0,0);
    for (int i = 0; i<MAXBIT; i++, graycode >>= 1) 
      if (graycode & lint(1,0)) ix[k] ^= iv[i * MAXDIM + k];
  }
}
