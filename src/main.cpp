 	#include "SVD.h"
	#include "stdio.h"
	#include <cstring>
	
	#define M 4
	#define N 3
	
	int main()
	{
	int i,j;
	double a[M*N]={1.0,1.0,-1.0,2.0,1.0,0.0,1.0,-1.0,0.0,-1.0,2.0,1.0};
	double b[M*N]={1.0,1.0,-1.0,-1.0,2.0,1.0,0.0,2.0,1.0,-1.0,0.0,1.0};
	double u[M*M],v[N*N],c[M*N],d[M*N];
	// double aa[4] = {3.0,1.0,4.0,3.0};
	int M_N = M*N;
	std::memset(u,0,sizeof(double)*M_N);
	double eps=0.000001;
	i=SVD(a,M,N,u,v,eps,M+1);
	printf("\n");
	printf("i=%d\n",i);
	printf("\nMAT U Is:\n");
	for(i=0;i<=N;i++)
	{
	for(j=0;j<=N;j++)
	printf("%e ",u[i*M+j]);
	printf("\n");
	}
	printf("\n");
	printf("MAT V IS:\n");
	for(i=0;i<N;i++)
	{
	for(j=0;j<N;j++)
	printf("%e ",v[i*N+j]);
	printf("\n");
	}
	printf("\n");
	printf("MAT A Is:\n");
	for(i=0;i<M;i++)
	{
	for(j=0;j<N;j++)
	printf("%e ",a[i*N+j]);
	printf("\n");
	}
	printf("\n----the nrom is :%f\n",a[0]);
	TriMul(u,a,M,M,N,c);
	TriMul(c,v,M,N,N,a);
	printf("\nMAT UAV IS:\n");
	for(i=0;i<M;i++)
	{
	for(j=0;j<N;j++)
	printf("%e ",a[i*N+j]);
	printf("\n");
	}
	printf("\n\n");
	return 0;
	
	}
