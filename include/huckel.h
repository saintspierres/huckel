/*
 *  huckel.h
 *  
 *
 *  Created by Eric Bremond on 12/10/11.
 *  Copyright 2011 huckel.h. All rights reserved.
 *
 */

#define TAB_SIZE 100
#define CRD_SIZE 3
#define LGN_SIZE 100

void rxyz(char *filename, int *natom, int *elts, double xyz[][CRD_SIZE]);
void satom(int natom, int *elts, int *nhelts, int *helts);
void mhuckel(int nhelts, int *helts, int *elts, double xyz[][CRD_SIZE], double hmat[][TAB_SIZE]);
void jacobi(double **a, int n, double d[], double **v, int *nrot);
void eigsrt(double d[], double **v, int n);
void diagonalize(double hmat[][TAB_SIZE], int nhelts, double *eigval, double eigvct[][TAB_SIZE]);
void vctp(double *vec1, double *vec2, double *vec12);
void dotp(double *vec1, double *vec2, double *vec12);
void rmat(double *axe, double *cth, double R[][CRD_SIZE]);
void rotate(int natom, double xyz[][CRD_SIZE]);
void comp2p(double Z, double x, double y, double z, double *psi);
void compom(int nhelts, int *helts, int *elts, double xyz[][CRD_SIZE],
    double x, double y, double z, double eigvct[][TAB_SIZE], int numorb,
    double *om);
void plot3d(int nhelts, int *helts, int natom, int *elts, double xyz[][CRD_SIZE],
    double eigvct[][TAB_SIZE], int numorb);
