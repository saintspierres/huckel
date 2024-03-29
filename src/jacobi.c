#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
        a[k][l]=h+s*(g-h*tau);

void jacobi(double **a, int n, double d[], double **v, double *nrot)
/*
Computes all eigenvalues and eigendvectors of a real symmetric ddmatrix a[1..n][1..n]. 
On output, elements of a above the diagonal are destroyed. d[1..n] returns the 
eigenvalues of a. v[1..n][1..n] is a ddmatrix whose columns contain, on output, 
the normalized eigendvectors of a. nrot returns the number of Jacobi rotations 
that were required.
*/
{
        int j,iq,ip,i;
        double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

        b=dvector(1,n);
        z=dvector(1,n);
        for (ip=1;ip<=n;ip++) {
                for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
                v[ip][ip]=1.0;
        }
        for (ip=1;ip<=n;ip++) {
                b[ip]=d[ip]=a[ip][ip];
                z[ip]=0.0;
        }
        *nrot=0;
        for (i=1;i<=50;i++) {
                sm=0.0;
                for (ip=1;ip<=n-1;ip++) {
                        for (iq=ip+1;iq<=n;iq++)
                                sm += fabs(a[ip][iq]);
                }
                if (sm == 0.0) {
                        free_dvector(z,1,n);
                        free_dvector(b,1,n);
                        return;
                }
                if (i < 4)
                        tresh=0.2*sm/(n*n);
                else
                        tresh=0.0;
                for (ip=1;ip<=n-1;ip++) {
                        for (iq=ip+1;iq<=n;iq++) {
                                g=100.0*fabs(a[ip][iq]);
                                if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
                                        && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
                                        a[ip][iq]=0.0;
                                else if (fabs(a[ip][iq]) > tresh) {
                                        h=d[iq]-d[ip];
                                        if ((double)(fabs(h)+g) == (double)fabs(h))
                                                t=(a[ip][iq])/h;
                                        else {
                                                theta=0.5*h/(a[ip][iq]);
                                                t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                                                if (theta < 0.0) t = -t;
                                        }
                                        c=1.0/sqrt(1+t*t);
                                        s=t*c;
                                        tau=s/(1.0+c);
                                        h=t*a[ip][iq];
                                        z[ip] -= h;
                                        z[iq] += h;
                                        d[ip] -= h;
                                        d[iq] += h;
                                        a[ip][iq]=0.0;
                                        for (j=1;j<=ip-1;j++) {
                                                ROTATE(a,j,ip,j,iq)
                                        }
                                        for (j=ip+1;j<=iq-1;j++) {
                                                ROTATE(a,ip,j,j,iq)
                                        }
                                        for (j=iq+1;j<=n;j++) {
                                                ROTATE(a,ip,j,iq,j)
                                        }
                                        for (j=1;j<=n;j++) {
                                                ROTATE(v,j,ip,j,iq)
                                        }
                                        ++(*nrot);
                                }
                        }
                }
                for (ip=1;ip<=n;ip++) {
                        b[ip] += z[ip];
                        d[ip]=b[ip];
                        z[ip]=0.0;
                }
        }
        nrerror("Too many iterations in routine jacobi");
}
#undef ROTATE
#undef NRANSI

void eigsrt(double d[], double **v, int n)
/*
Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] 
as output from jacobi (§11.1) or tqli (§11.3), this routine 
sorts the eigenvalues into descending order, and rearranges 
the columns of v correspondingly. The method is straight 
insertion.
*/
{
        int k,j,i;
        double p;

        for (i=1;i<n;i++) {
                p=d[k=i];
                for (j=i+1;j<=n;j++)
                        if (d[j] >= p) p=d[k=j];
                if (k != i) {
                        d[k]=d[i];
                        d[i]=p;
                        for (j=1;j<=n;j++) {
                                p=v[j][i];
                                v[j][i]=v[j][k];
                                v[j][k]=p;
                        }
                }
        }
}
