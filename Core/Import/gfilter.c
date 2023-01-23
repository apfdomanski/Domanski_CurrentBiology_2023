/* (c) 2004 Daniel Durstewitz */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stddef.h>
#include "mex.h"

#define pi 3.14159265358979


void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
  /*    double *v,*vf,*ff,ffsum;
	long int N1,N2,N,L,s,i,j,m,k,k1; */
  double *v,*vf,*ff,ffsum,s;
  long int N1,N2,N,L,i,j,m,k,k1;

    if (nrhs!=3) mexErrMsgTxt("3 inputs required: v,s,L");
    v=mxGetPr(prhs[0]);
    N1=mxGetM(prhs[0]);
    N2=mxGetN(prhs[0]);
    N=fmax(N1,N2);
    s=*(mxGetPr(prhs[1]));
    L=*(mxGetPr(prhs[2]));
    plhs[0]=mxCreateDoubleMatrix(N,1,mxREAL);
    vf=mxGetPr(plhs[0]);
    m=(int)(L/2);
    L=2*m+1;
    ff=(double *) mxMalloc(L*sizeof(double));

    /*    for (i=0,ffsum=0.0;i<L;i++) {
	ff[i]=exp(-(i-m)*(i-m)/(2.0*s*s))/(sqrt(2*pi)*s);
	ffsum+=ff[i]; }
    for (i=0;i<L;i++) ff[i]/=ffsum;
    for (i=0;i<N;i++) {
	vf[i]=0.0;
	k=fmax(0,i-m);
	k1=fmax(0,m-i)-k;
	for (j=k;j<=fmin(N-1,i+m);j++)
	vf[i]+=v[j]*ff[j+k1]; } */

    for (i=0;i<L;i++) ff[i]=exp(-(i-m)*(i-m)/(2.0*s*s))/(sqrt(2*pi)*s);
    for (i=0;i<N;i++) {
	vf[i]=0.0;
	k=fmax(0,i-m);
	k1=fmax(0,m-i)-k;
	ffsum=0.0;
	for (j=k;j<=fmin(N-1,i+m);j++) ffsum+=ff[j+k1];
	for (j=k;j<=fmin(N-1,i+m);j++) vf[i]+=v[j]*ff[j+k1]/ffsum;
    }

    return;
}
