/*
 *  huckel.c
 *  
 *
 *  Created by Eric Bremond on 12/10/11.
 *  Copyright 2011 huckel.c. All rights reserved.
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<nrutil.h>

#include"huckel.h"

int main(int argc, char *argv[])
{
    int i, j;

    int natom  = 0;
    int nhelts = 0;
    int elts[TAB_SIZE]  = {0};
    int helts[TAB_SIZE] = {0};
    
    double xyz[TAB_SIZE][CRD_SIZE]  = {{0.}};
    double hmat[TAB_SIZE][TAB_SIZE] = {{0.}};

    int nrot = 0;
    double eigval[TAB_SIZE] = {0};
    double eigvct[TAB_SIZE][TAB_SIZE] = {{0}};

    rxyz(argv[1], &natom, elts, xyz);
/*
    printf(" %s: %d atoms\n", argv[1], natom);

    for (i=0; i<natom; i++)
    {
       printf("%12.6lf%12.6lf%12.6lf\n", xyz[i][0], xyz[i][1], xyz[i][2]);
    }
*/
    satom(natom, elts, &nhelts, helts);

    mhuckel(nhelts, helts, elts, xyz, hmat);
/*
    printf("Input matrix\n");

    for (i=0; i<nhelts; i++)
    {
        for (j=0; j<nhelts; j++)
        {
            printf("%12.6lf", hmat[i][j]);
        }
        printf("\n");
    }
*/
    diagonalize(hmat, nhelts, eigval, eigvct);
/*
    printf("Eigen vectors\n");

    for (i=0; i<nhelts; i++)
    {
        for (j=0; j<nhelts; j++)
        {
            printf("%12.6lf", eigvct[i][j]);
        }
        printf("\n");
    }

    printf("Eigen values\n");

    for (i=0; i<nhelts; i++)
    {
        printf("%12.6lf", eigval[i]);
    }

    printf("\n");

    printf("Rotation of the molecule\n");
*/
    rotate(natom, xyz);
/*
    for (i=0; i<natom; i++)
    {
       printf("%12.6lf%12.6lf%12.6lf\n", xyz[i][0], xyz[i][1], xyz[i][2]);
    }
*/
    plot3d(nhelts, helts, natom, elts, xyz, eigvct, 0);
}

void rxyz(char *filename, int *natom, int *elts, double xyz[][CRD_SIZE])
{
    int i;
    char lign[LGN_SIZE] = "";

    FILE *pfile = NULL;
    pfile = fopen(filename, "r");

    if (pfile==NULL)
    {
        fprintf(stderr, "huckel: %s: No such file or directory\n", filename);
        exit(0);
    }

    fscanf(pfile, "%d", natom);
    fgets(lign, LGN_SIZE, pfile);

    for (i=0; i<*natom; i++)
    {
        fscanf(pfile, "%d%lf%lf%lf",  &elts[i], &xyz[i][0], &xyz[i][1], &xyz[i][2]);
    }

    fclose(pfile);
}

void satom(int natom, int *elts, int *nhelts, int *helts)
{

    int i, n = 0;

    for (i=0; i<natom; i++)
    {
        if (elts[i]!=1)
        {
            helts[i] = i;
            n+=1;
        }
    }

    *nhelts=n;
}

void mhuckel(int nhelts, int *helts, int *elts, double xyz[][CRD_SIZE], double hmat[][TAB_SIZE])
{
    int i, j, k;

    double hhuk[11] =      {0.,0.,0.,0.,0.,0.,0.00,1.37,0.97,0.,0.};
    double khuk[11][11] = {{0.,0.,0.,0.,0.,0.,0.00,0.00,0.00,0.,0.},
                           {0.,0.,0.,0.,0.,0.,0.00,0.00,0.00,0.,0.},
                           {0.,0.,0.,0.,0.,0.,0.00,0.00,0.00,0.,0.},
                           {0.,0.,0.,0.,0.,0.,0.00,0.00,0.00,0.,0.},
                           {0.,0.,0.,0.,0.,0.,0.00,0.00,0.00,0.,0.},
                           {0.,0.,0.,0.,0.,0.,0.00,0.00,0.00,0.,0.},
                           {0.,0.,0.,0.,0.,0.,1.00,0.89,1.06,0.,0.},
                           {0.,0.,0.,0.,0.,0.,0.89,0.98,1.13,0.,0.},
                           {0.,0.,0.,0.,0.,0.,1.06,1.13,1.26,0.,0.},
                           {0.,0.,0.,0.,0.,0.,0.00,0.00,0.00,0.,0.},
                           {0.,0.,0.,0.,0.,0.,0.00,0.00,0.00,0.,0.}};
    double rvdw[11] = {0.00,0.31,0.28,1.28,0.96,0.84,0.76,0.00,0.66,0.57,0.58};

    double rclc = 0.;
    double rthr = 0.;

    for (i=0; i<nhelts; i++)
    {
        for (j=0; j<nhelts; j++)
        {
            rclc=0.;
            for (k=0; k<CRD_SIZE; k++)
            {
                rclc+=pow(xyz[helts[i]][k]-xyz[helts[j]][k], 2.);
            }
            rclc=pow(rclc, .5);

            rthr=rvdw[elts[helts[i]]] + rvdw[elts[helts[j]]];

            if (rclc<=1.1*rthr)
            {
                if (i==j)
                {
                    hmat[i][j]=hhuk[elts[helts[i]]];
                }
                else
                {
                    hmat[i][j]=khuk[elts[helts[i]]][elts[helts[j]]];
                }
            }
        }
    }
}

void diagonalize(double hmat[][TAB_SIZE], int nhelts, double *eigval, double eigvct[][TAB_SIZE])
{
    int i, j;

    int nrot=0;

    double **a, **v, *d;

    a=dmatrix(1,nhelts,1,nhelts);
    v=dmatrix(1,nhelts,1,nhelts);
    d=dvector(1,nhelts);

    for (i=1; i<=nhelts; i++)
    {
        for (j=1; j<=nhelts; j++)
        {
            a[i][j]=hmat[i-1][j-1];
        }
    }

    jacobi(a, nhelts, d, v, &nrot);
    eigsrt(d, v, nhelts);

    for (i=1; i<=nhelts; i++)
    {
        eigval[i-1]=d[i];

        for (j=1; j<=nhelts; j++)
        {
            eigvct[i-1][j-1]=v[i][j];
        }
    }

    free_dmatrix(a,1,nhelts,1,nhelts);
    free_dmatrix(v,1,nhelts,1,nhelts);
    free_dvector(d,1,nhelts);
}

void vctp(double *vec1, double *vec2, double *vec12)
{
    int i, j;

    double norm=0.;

    vec12[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    vec12[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    vec12[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];

    dotp(vec12, vec12, &norm);
    norm = pow(norm, .5);

    for(i=0; i<CRD_SIZE; i++)
    {
        vec12[i]=vec12[i]/norm;
    }
}

void dotp(double *vec1, double *vec2, double *vec12)
{
    int i;

    *vec12=0.;
    for(i=0; i<CRD_SIZE; i++)
    {
        *vec12+=vec1[i]*vec2[i];
    }
}

void rmat(double *axe, double *cth, double R[][CRD_SIZE])
{

    int i, j;

    double sth = 0.;

    double P[CRD_SIZE][CRD_SIZE] = {{0.}};
    double I[CRD_SIZE][CRD_SIZE] = {{0.}};
    double Q[CRD_SIZE][CRD_SIZE] = {{0.,-axe[2],axe[1]},
                                    {axe[2],0.,-axe[0]},
                                    {-axe[1],axe[0],0.}};

    sth=pow(1. - pow(*cth, 2.), .5);

    for(i=0; i<CRD_SIZE; i++)
    {
        for(j=0; j<CRD_SIZE; j++)
        {
            P[i][j] = axe[i] * axe[j];

            if(i==j)
            {
                I[i][j]=1.;
            }
        }
    }

    for(i=0; i<CRD_SIZE; i++)
    {
        for(j=0; j<CRD_SIZE; j++)
        {
            R[i][j] = P[i][j]+(I[i][j]-P[i][j])*(*cth)+Q[i][j]*sth;
        }
    }

}

void rotate(int natom, double xyz[][CRD_SIZE])
{
    int i, j, k;

    double a0 = .52917720859;

    double ref[CRD_SIZE]={0.,0.,1.}, axe[CRD_SIZE]={0.};
    double R[CRD_SIZE][CRD_SIZE] = {{0.}};
    double u[CRD_SIZE]={0.}, v[CRD_SIZE]={0.}, n[CRD_SIZE]={0.};
    double cth=0.;
    double rxyz[TAB_SIZE][CRD_SIZE] = {{0.}};

    double G[CRD_SIZE] = {0.};

    // compute the normal vector n of the plan
    for(i=0; i<CRD_SIZE; i++)
    {
        u[i] = xyz[0][i]-xyz[1][i];
        v[i] = xyz[0][i]-xyz[2][i];
    }

    vctp(u, v, n);

    // compute the cosine angle cth between n and ref
    dotp(n, ref, &cth);

    // compute the axe of the rotation
    vctp(n, ref, axe);

    // compute the rotation matrix R
    rmat(axe, &cth, R);

    // rotate the molecule
    for(i=0; i<natom; i++)
    {
         for(j=0; j<CRD_SIZE; j++)
         {
             rxyz[i][j]=0.;
             for(k=0; k<CRD_SIZE; k++)
             {
                 rxyz[i][j]+=R[j][k]*xyz[i][k];
             }
         }
    }

    // compute the mass center and translate
    for(i=0; i<natom; i++)
    {
        for(j=0; j<CRD_SIZE; j++)
        {
            G[j]+=rxyz[i][j];
        }
    }
    for(i=0; i<CRD_SIZE; i++)
    {
        G[i]=G[i]/((double)natom);
    }

    for(i=0; i<natom; i++)
    {
         for(j=0; j<CRD_SIZE; j++)
         {
             xyz[i][j]=(rxyz[i][j]-G[j])/a0;
         }
    }

}

void comp2p(double Z, double x, double y, double z, double *psi)
{
    double r=0.;

    r = sqrt( pow(x, 2.) + pow(y, 2.) + pow(z, 2.) );

    *psi= z * Z * exp(-r*Z/2.);
}

void compom(int nhelts, int *helts, int *elts, double xyz[][CRD_SIZE], 
    double x, double y, double z, double eigvct[][TAB_SIZE], int numorb, 
    double *om)
{
    int i;

    double psi=0.;

    *om=0.;
    for(i=0; i<nhelts; i++)
    {
        x=x-xyz[helts[i]][0];
        y=y-xyz[helts[i]][1];
        z=z-xyz[helts[i]][2];
        comp2p((double) elts[helts[i]], x, y, z, &psi);
        *om+=eigvct[i][numorb]*psi;
    }
}

void plot3d(int nhelts, int *helts, int natom, int *elts, double xyz[][CRD_SIZE],
    double eigvct[][TAB_SIZE], int numorb)
{
    int i, j, k, l;
    int Nx=80, Ny=80, Nz=80;

    double x0=-8.4, y0=-8.4, z0=-8.4;

    double Stpx[CRD_SIZE] = {.2,0.,0.};
    double Stpy[CRD_SIZE] = {0.,.2,0.};
    double Stpz[CRD_SIZE] = {0.,0.,.2};

    double x, y, z, om;

    printf("orbitale 1s\nMO coefficients\n");
    printf("%5d%12.6f%12.6f%12.6f\n", -natom, x0, y0, z0);
    printf("%5d%12.6f%12.6f%12.6f\n", Nx, Stpx[0], Stpx[1], Stpx[2]);
    printf("%5d%12.6f%12.6f%12.6f\n", Ny, Stpy[0], Stpy[1], Stpy[2]);
    printf("%5d%12.6f%12.6f%12.6f\n", Nz, Stpz[0], Stpz[1], Stpz[2]);
    for(i=0; i<natom; i++)
    {
        printf("%5d%12.6lf%12.6f%12.6f%12.6f\n", elts[i], (double) elts[i], xyz[i][0], xyz[i][1], xyz[i][2]);
    }
    printf("%5d%5d\n", 1, 1);

    for(i=0; i<Nx; i++)
    {
        for(j=0; j<Ny; j++)
        {
            l=0;
            for(k=0; k<Nz; k++)
            {
                x=x0+i*Stpx[0]+j*Stpx[1]+k*Stpx[2];
                y=y0+i*Stpy[0]+j*Stpy[1]+k*Stpy[2];
                z=z0+i*Stpz[0]+j*Stpz[1]+k*Stpz[2];
                compom(nhelts, helts, elts, xyz, x, y, z, eigvct, numorb, &om);
                //printf("%12.6lf%12.6lf%12.6lf%13.5lE\n", x, y, z, om);
                printf("%13.5lE", om);
                l++;
                if(l==6)
                {
                    printf("\n");
                    l=0;
                }
            }
            printf("\n");
        }
    }
}
