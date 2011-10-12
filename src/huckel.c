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

    printf(" %s: %d atoms\n", argv[1], natom);

    for (i=0; i<natom; i++)
    {
       printf("%12.6lf%12.6lf%12.6lf\n", xyz[i][0], xyz[i][1], xyz[i][2]);
    }

    satom(natom, elts, &nhelts, helts);

    mhuckel(nhelts, helts, elts, xyz, hmat);

    printf("Input matrix\n");

    for (i=0; i<nhelts; i++)
    {
        for (j=0; j<nhelts; j++)
        {
            printf("%12.6lf", hmat[i][j]);
        }
        printf("\n");
    }

    diagonalize(hmat, nhelts, eigval, eigvct);

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
