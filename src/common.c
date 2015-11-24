//source file for communication related functions
//#define GETRHS_DEBUG

#include "common.h"



void binaryAllReduceButterflyPrototype(int myid, int numprocs)
{
	MPI_Status status;
	int sendValue = myid + 1;
	//int sendValue = 1;
	int partnerId;
	int stageCounter;
	printf("id=%d, numproces=%d\n",myid,numprocs);
	int totalLevelNum = (int)log2(numprocs);
	printf ("totalLeveNum= %d\n", totalLevelNum);

	for (stageCounter=0; stageCounter<totalLevelNum; stageCounter++)
	{
		int groupSize = pow (2,stageCounter+1);
		int receiveTemp;
		partnerId = groupSize *  
			(int)floor ( (double)myid / (double) groupSize) +
			(myid%groupSize + groupSize/2) % groupSize;
		if (myid < partnerId)	
		{
			// send
			MPI_Send ((void*)&sendValue, 1, MPI_INT, partnerId, 1, MPI_COMM_WORLD);
			// receive
			MPI_Recv ((void*)&receiveTemp, 1, MPI_INT, partnerId, 1, MPI_COMM_WORLD, &status);
		}
		else
		{
			// receive
			MPI_Recv ((void*)&receiveTemp, 1, MPI_INT, partnerId, 1, MPI_COMM_WORLD, &status);
			// send
			MPI_Send ((void*)&sendValue, 1, MPI_INT, partnerId, 1, MPI_COMM_WORLD);
		}
		sendValue += receiveTemp;
	}
	
	printf ("myId=%d, finalSum=%d \n", myid, sendValue);	
	//MPI_Barrier(MPI_COMM_WORLD);
}

void parseCSV_rowMajor(char* filename, double** output, int numRows, int numCols) {

    int buffersize = 10240;
    char buf[buffersize];


    // make sure we have buf long enough for reading each line
//      int tsize = (sizeof (double) + sizeof (char)) * numCols;
//      printf ("tsize is %d",tsize);
    assert((sizeof (double) + sizeof (char)) * numCols < buffersize);

    FILE * fstream = fopen(filename, "r");
    assert(fstream != 0);
//      printf ("good i am here file opend \n");
    *output = (double*) calloc(numRows*numCols, sizeof (double));

    int rowCounter = 0;
    const char * tok;
    int outputCounter = 0;
    while (rowCounter < numRows && fgets(buf, sizeof (buf), fstream)) {
        //        printf("%s", buf);
        // parse a line of CSV
        int rowEleCounter = 0;
        tok = strtok(buf, ",");
        //        printf("%s\n",tok);
        for (; tok&& *tok; tok = strtok(NULL, ",\n")) {
            // remember *output must be parenthesised
            (*output)[outputCounter++] = atof(tok);
            rowEleCounter++;
        }

        assert(rowEleCounter == numCols);

        rowCounter++;
    }


    fclose(fstream);
        printf ("good i am here file closed \n");
}

void printMatrix (denseType data, int myid, int numprocs)
{
        int rowIdx;
        int colIdx;
        int datacounter=0;
	printf ("Hi, I am rank:%d \n", myid);
        for (rowIdx = 0; rowIdx<data.local_num_row; rowIdx++)
        {
                for (colIdx=0; colIdx<data.local_num_col; colIdx++)
                {
                        printf ("%f,",data.data[datacounter++]);
                }
                printf ("\n");
        }
}

void printMatrix_rowMajor (double * M, int numRows, int numCols, int myid, int numprocs)
{
        int rowIdx;
        int colIdx;
        int datacounter=0;
	printf ("Hi, I am rank:%d \n", myid);
        for (rowIdx = 0; rowIdx<numRows; rowIdx++)
        {
                for (colIdx=0; colIdx<numCols; colIdx++)
                {
                        printf ("%f,",M[datacounter++]);
                }
                printf ("\n");
        }
	
}

void dataDist_rowMajor(denseType * vector, int numRows, int numCols, char* rhsFile, int myid, int numprocs) {
    int idx;
    int local_length, local_length_normal;
    int ierr;
    double * Total_data_buffer;
    int sendCount[numprocs];
    int sendDispls[numprocs];
    int procCounter;
    double normel_ele_num;

    local_length_normal = numRows / numprocs;
    //vector->start_idx = myid * vector->global_num_col * local_length_normal;
    vector->start_idx = myid * local_length_normal;

    if (myid == numprocs - 1)
    {
        local_length = numRows - (numprocs - 1) * local_length_normal;
    }
    else
    {
        local_length = local_length_normal;
    }

    normel_ele_num = local_length_normal * numCols;

    for (procCounter = 0; procCounter < numprocs; procCounter++) {
        sendCount[procCounter] = normel_ele_num;
        sendDispls[procCounter] = procCounter * normel_ele_num;
    }
    sendCount[numprocs - 1] = (numRows - (numprocs - 1) * local_length_normal) * numCols;

#ifdef GETRHS_DEBUG
    printf("in GetRHS.c, myid = %d, local_length=%d\n", myid, local_length);
#endif
    vector->local_num_row = local_length;
    vector->local_num_col = numCols; // only consider
    vector->global_num_row = numRows;
    vector->global_num_col = numCols;

    vector->data = (double *) calloc(vector->local_num_row * vector->local_num_col, sizeof (double));

    int local_num_element = vector->local_num_row * vector->local_num_col;

    // rank 0 read CSV
    if (myid == 0) {
        printf("Reading input matrix  data from %s ... ...\n", rhsFile);
        parseCSV_rowMajor(rhsFile, &Total_data_buffer, numRows, numCols);
	
        printf("Reading input matrix  data from %s done.\n", rhsFile);
#ifdef GenVector_ReadCSV_DB
        //check_csv_array_print(Total_data_buffer, length, num_cols, myid);
        //        exit(0);
#endif

    }
    //    // Scatter data
    ierr = MPI_Scatterv((void*) Total_data_buffer, sendCount, sendDispls,
            MPI_DOUBLE, vector->data, local_num_element,
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //
    // based on the assumption of equal division of rows among processes 
    //vector->start_idx = myid * vector->global_num_col * local_length_normal;

#ifdef GenVector_ReadCSV_DB
    //    if (myid == 0){
    //        check_csv_array_print(vector->data, vector->local_num_row, vector->global_num_col, myid);
    //        printf ("local rows:%d, local cols %d, local nnz:%d\n", vector->local_num_row, vector->local_num_col,local_num_element);
    //    }
    if (myid == numprocs-1) {
        local_dense_mat_print(*vector, myid);
    }
    exit(0);
#endif    
    if (myid == 0) {
        free(Total_data_buffer);
    }

}

void localMatrixCopy(double ** target, double * src, int size)
{
	int eleIdx = 0;
	(*target) = calloc (size, sizeof(double)); 	
	for (; eleIdx<size;eleIdx++)	
	{
		(*target)[eleIdx] = src[eleIdx];	
	} 
}

void saveData_SendBuffer_RestoreR(double * sendBuffer, double * assemblizedR, double *  src, int numCols,int myid, int partnerid)
{
	int commonRowIdx, commonColIdx;
	int sentSizeRIdx=0;	
	for (commonRowIdx=0; commonRowIdx<numCols;commonRowIdx++)
	{
		for (commonColIdx=commonRowIdx;commonColIdx<numCols;commonColIdx++)
		{
			int tempRIdx, tempIdx;
			tempIdx=commonRowIdx*numCols+commonColIdx;

			if (myid>partnerid) //lower part of assembled R
			{
				tempRIdx=numCols * numCols + tempIdx;
			}
			else
			{
				tempRIdx = tempIdx;
			}
			double tempRval = src[tempIdx];
			sendBuffer[sentSizeRIdx++] = tempRval;		
			assemblizedR[tempRIdx] = tempRval; 
		}
	}

}

void saveFinalR_rowOrder (denseType finalR, double * src)
{
	int size = finalR.local_num_col;
	int commonRowIdx, commonColIdx;
	for (commonRowIdx=0; commonRowIdx<size;commonRowIdx++)
	{
		for (commonColIdx=commonRowIdx;commonColIdx<size;commonColIdx++)
		{
			int tempIdx;
			tempIdx = commonRowIdx * size + commonColIdx;
			finalR.data[tempIdx] = src[tempIdx]; 
		}
	}
}

void initDenseType_rowOrder (denseType* M, int localNumRows, int localNumCols, int globalNumRows, int globalNumCols,int startIdx)
{
	M->data = calloc ( localNumRows*localNumCols, sizeof(double));
	M->local_num_row = localNumRows;
	M->local_num_col = localNumCols;
	M->global_num_row = globalNumRows;
	M->global_num_col = globalNumCols;
	M->start_idx = startIdx;
}

void matrixChainMM_homo(double ** matrixChain, int chainSize, double * result, int numRows, int numCols)
{
	int totalElements = numRows*numCols;
	double * buffer_0 = calloc ( totalElements,sizeof(double));
	double * buffer_1 = calloc ( totalElements,sizeof(double));

	if ( chainSize < 2)
	{
		return;
	}
	cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, numRows, numCols,numCols
			, 1.0, matrixChain[0], numCols, matrixChain[1], numCols
			, 0.0, buffer_0, numCols);	

	int matrixIdx;
	for (matrixIdx = 2; matrixIdx<chainSize; matrixIdx++)
	{
		if (matrixIdx%2 == 0)
		{
			cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, numRows, numCols,numCols
					, 1.0, buffer_0, numCols, matrixChain[matrixIdx], numCols
					, 0.0, buffer_1, numCols);	
		}
		else 
		{
			cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, numRows, numCols,numCols
					, 1.0, buffer_1, numCols, matrixChain[matrixIdx], numCols
					, 0.0, buffer_0, numCols);	
		}
	}

	//copy final result 
	int eleIdx;
	if (chainSize % 2 == 0)
	{
		for (eleIdx=0;eleIdx < totalElements; eleIdx++)
			result[eleIdx] = buffer_0[eleIdx];	
	}
	else
	{
		for (eleIdx=0;eleIdx < totalElements; eleIdx++)
			result[eleIdx] = buffer_1[eleIdx];	
	}

	// free buffers
	free (buffer_0);
	free (buffer_1);
}

void restoreDataFromReceiveBufferToIntermediateQRMatrix(double * M, double * receiveBuffer, int numCols,int myid, int partnerid)
{
	int commonRowIdx, commonColIdx;
	int receiveBufferIdx=0;	
	for (commonRowIdx=0; commonRowIdx<numCols;commonRowIdx++)
	{
		for (commonColIdx=commonRowIdx;commonColIdx<numCols;commonColIdx++)
		{
			int tempIdx;
			tempIdx=commonRowIdx*numCols+commonColIdx;

			if (myid<partnerid) //lower part of assembled R
			{
				tempIdx += numCols * numCols;
			}
			double tempVal = receiveBuffer[receiveBufferIdx++];
			M[tempIdx] = tempVal;		
		}
	}

}

double tsqrResCheck ( denseType A, denseType Q, denseType R)
{
	int localMatrixSize = A.local_num_row * A.local_num_col;
	double * resBuffer = calloc(localMatrixSize, sizeof(double));

	
	cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, Q.local_num_row, R.local_num_col, Q.local_num_col
			, 1.0, Q.data, Q.local_num_col, R.data, R.local_num_col 
			, 0.0, resBuffer, Q.local_num_col);	

	int idx;
	double resultDist = 0.0;
	for (idx=0; idx<localMatrixSize; idx++)	
	{
		double gap = resBuffer[idx] - A.data[idx];
		resultDist += gap * gap;
	}

	free (resBuffer);

	return sqrt(resultDist);	

}

