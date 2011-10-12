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
