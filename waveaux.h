/*
Copyright (c) 2014, Rafat Hussain
Copyright (c) 2016, Holger Nahrstaedt
*/
#ifndef WAVEAUX_H_
#define WAVEAUX_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#ifdef __cplusplus
extern "C" {
#endif


int filtlength(const char* name);

int filtcoef(const char* name, double *lp1, double *hp1, double *lp2, double *hp2);

void copy_reverse(const double *in, int N, double *out);

void qmf_even(const double *in, int N, double *out);

void qmf_wrev(const double *in, int N, double *out);

void copy(const double *in, int N, double *out);

int wmaxiter(int sig_len, int filt_len);

void conv_direct(double *inp1, int N, double *inp2, int L, double *oup);

#ifdef __cplusplus
}
#endif


#endif /* WAVEAUX_H_ */