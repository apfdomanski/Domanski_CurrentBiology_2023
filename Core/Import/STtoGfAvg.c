/* (c) 2006 Daniel Durstewitz, CTCN, Univ. Of Plymouth
   converts spike train matrix into Gaussian-filtered activity matrix
   NOTE: ST=0 is interpreted as an end-of-spike-train signal!!!
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stddef.h>
#include "mex.h"

#define miss 0
#define cutoff 1.0e-7
#define pi 3.14159265358979
#define resol 100000

/* ----- NR components --------------------------------------------------------------- */

#define NR_END 1
#define FREE_ARG char*
static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

long *ivector(long nl, long nh)
/* allocate an long vector with subscript range v[nl..nh] */
{
	long *v;

	v=(long *)mxMalloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) mexErrMsgTxt("allocation failure in ivector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)mxMalloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) mexErrMsgTxt("allocation failure in dvector()");
	return v-nl+NR_END;
}

void free_ivector(long *v, long nl, long nh)
/* free an long vector allocated with ivector() */
{
	mxFree((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	mxFree((FREE_ARG) (v+nl-NR_END));
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) mxMalloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) mexErrMsgTxt("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) mxMalloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) mexErrMsgTxt("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

long **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a long matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	long **m;

	/* allocate pointers to rows */
	m=(long **) mxMalloc((size_t)((nrow+NR_END)*sizeof(long*)));
	if (!m) mexErrMsgTxt("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(long *) mxMalloc((size_t)((nrow*ncol+NR_END)*sizeof(long)));
	if (!m[nrl]) mexErrMsgTxt("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	mxFree((FREE_ARG) (m[nrl]+ncl-NR_END));
	mxFree((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(long **m, long nrl, long nrh, long ncl, long nch)
/* free an long matrix allocated by imatrix() */
{
	mxFree((FREE_ARG) (m[nrl]+ncl-NR_END));
	mxFree((FREE_ARG) (m+nrl-NR_END));
}

/* ----------------------------------------------------------------------------------- */


void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
  double *STMtx,*Tinv,sGf,*GfST,dt,*CumNDF,x1,x2,binw,T2,k1,k2,s1,s2;
    long i,j,m,k,Ncells,Nst,Lpdf,Lmax;

    if (nrhs<3) mexErrMsgTxt("at least 3 inputs required: STMtx,Tinv,sGf,(binw)");
    STMtx=mxGetPr(prhs[0]);
    Ncells=(long)mxGetM(prhs[0]);
    Nst=(long)mxGetN(prhs[0]);
    Tinv=mxGetPr(prhs[1]);
    i=(long)mxGetM(prhs[1]);
    j=(long)mxGetN(prhs[1]);
    m=(long)FMAX(i,j)-1;
    sGf=*(mxGetPr(prhs[2]));
    if (nrhs>3) binw=*(mxGetPr(prhs[3]));
    dt=sGf/resol;
    plhs[0]=mxCreateDoubleMatrix(Ncells,m,mxREAL);
    GfST=mxGetPr(plhs[0]);

    Lmax=(long)(sqrt(-log(cutoff)*2)*resol-0.5+1.0);
    Lpdf=(long)((Tinv[m]-Tinv[0])/dt+1);
    Lpdf=FMIN(Lpdf,Lmax);
    CumNDF=dvector(0,Lpdf-1);
    CumNDF[Lpdf-1]=dt*exp(-(Lpdf-0.5)*dt*(Lpdf-0.5)*dt/(2.0*sGf*sGf))/(sqrt(2.0*pi)*sGf);
    for (i=Lpdf-2;i>=0;i--) CumNDF[i]=CumNDF[i+1]+dt*
	exp(-(i+0.5)*dt*(i+0.5)*dt/(2.0*sGf*sGf))/(sqrt(2.0*pi)*sGf);

    for (i=0;i<Ncells;i++)
	for (j=0;j<m;j++) {
	    GfST[i+j*Ncells]=0.0;
	    if (nrhs>3) T2=Tinv[j]+binw;
	    else T2=Tinv[j+1];
	    for (k=0;k<Nst;k++)
		if ((STMtx[i+k*Ncells]>Tinv[j]-Lpdf*dt) &&
		    (STMtx[i+k*Ncells]<T2+Lpdf*dt)) {
		    if (STMtx[i+k*Ncells]<=miss) break;
		    if (Tinv[j]<STMtx[i+k*Ncells]) s1=-0.5; else s1=0.5;
		    if (T2<STMtx[i+k*Ncells]) s2=-0.5; else s2=0.5;
		    k1=(double)((Tinv[j]-STMtx[i+k*Ncells])/dt+s1);
		    k2=(double)((T2-STMtx[i+k*Ncells])/dt+s2);

		    /* if (i+j*Ncells==200+1) mexPrintf("%lf %lf\n",k1,k2);
		       if (i+j*Ncells==200+1) mexPrintf("%ld %ld\n",(long)k1,(long)k2); */

		    if (fabs(k1)<Lpdf) x1=CumNDF[(long)fabs(k1)]; else x1=0.0;
		    if (fabs(k2)<Lpdf) x2=CumNDF[(long)fabs(k2)]; else x2=0.0;
		    if (k1*k2>=0.0) GfST[i+j*Ncells]+=fabs(x1-x2);
		    else GfST[i+j*Ncells]+=(1.0-x1-x2);
		    /* else GfST[i+j*Ncells]+=(2*CumNDF[0]-x1-x2); */
		}
	    GfST[i+j*Ncells]/=T2-Tinv[j];
	}

    free_dvector(CumNDF,0,Lpdf-1);
    return;
}
