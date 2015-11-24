//header file for communication pattern related functions
#include <stdio.h>
#include "mpi.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include "mkl.h"


typedef struct
{
	int local_num_row;
	int local_num_col;
        int global_num_row;
        int global_num_col;
        int start_idx;
        
	double * data;
} denseType;


void binaryAllReduceButterflyPrototype(int, int);

void parseCSV_rowMajor(char* , double** , int , int ); 

void printMatrix (denseType, int, int) ;

void printMatrix_rowMajor (double * , int , int , int , int );

void dataDist_rowMajor(denseType * , int, int, char* , int , int ); 

void localMatrixCopy (double **, double * , int);

void saveData_SendBuffer_RestoreR(double * , double * , double * , int ,int , int );

void initDenseType_rowOrder (denseType* , int , int , int , int ,int);

void matrixChainMM_homo(double ** , int , double * , int , int );

double tsqrResCheck ( denseType , denseType , denseType );

