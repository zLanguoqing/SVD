#ifndef _SVD_H
#define _SVD_H
#define MINNum (1e-30)

#include "stdlib.h"
#include "math.h"
#include<vector>
#include<iostream>

int SVD(double a[],int m,int n,double u[],double v[],double eps,int ka);
void QRCompute(double a[],double e[],double s[],double v[],int m,int n);
void Givens(double fg[2],double cs[2]);
void TriMul(double a[], double b[], int m, int n, int k, double c[]);
#endif
