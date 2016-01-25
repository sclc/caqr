// for tsqr all reduction routines

#include "tsqrAllReduction.h"

//phrase 1, 1st QR,the result Q isn't necessary to send, only R comm is necessary
//pharse 2, 2nd-log(numprocs)th QR, processes need to exchange both Q and R
//phrase 3, final phrease, all processes have the Q from the 1st phrase and the same copy of Qlist, and final R
 

void tsqrAllReduction_rowOrder_butterFly_raw (denseType localA, denseType finalQ, denseType finalR, int myid, int numprocs, double * elapsedComm)
{
	MPI_Status status;

	int localNumRows = localA.local_num_row;
	int localNumCols = localA.local_num_col;	
	int commonRowIdx, commonColIdx;

	//create buffer for finalR and finalQ
	//double * finalRData = calloc (localNumCols * localNumCols, sizeof(double));
	//double * finalQData = calloc (localNumRows * localNumCols, sizeof(double));

	//in order not to change the original input matrix, we create a buffer right now
	double * inputMatrixBuffer;
	localMatrixCopy (&inputMatrixBuffer, localA.data, localNumRows*localNumCols);
	double * firstPhraseQ= inputMatrixBuffer;

	//create buffer to restore R
	int sizeIntermediateQ = 2 * localNumCols * localNumCols;
	double * sendBuffer;
	double * receiveBuffer;

	//householder buffer tau
	double * tau = calloc (localNumCols, sizeof(double));
	//double * tau = calloc (localNumRows, sizeof(double));

	//communication the first R, or say, first phrase of butterfly
	int commStageCounter;
	int numTotalStage = (int)log2(numprocs);

	// create send receive buffer 

	int sizeSendReceiveBuffer;
	sizeSendReceiveBuffer = localNumCols * (localNumCols+1) / 2;
	sendBuffer = calloc (sizeSendReceiveBuffer, sizeof(double)) ;
	receiveBuffer = calloc (sizeSendReceiveBuffer, sizeof(double)) ;

	// allocate space for intermediate Q list
	double * intermediateQRList[numTotalStage];	
	double * intermediateQListForRetrievingQ[numTotalStage];

	int sizeIntermediateQR = 2 * localNumCols * localNumCols;
	int sizeIntermediateQForRetrievingQ = localNumCols * localNumCols;

	for (commStageCounter = 0; commStageCounter<numTotalStage; commStageCounter++)
	{
		intermediateQRList[commStageCounter] = calloc (sizeIntermediateQ, sizeof(double));
	//	intermediateQListForRetrievingQ[commStageCounter] = calloc (sizeIntermediateQ, sizeof(double));
	}

	double * QChainMMResult = calloc (sizeIntermediateQForRetrievingQ, sizeof(double));
	
	int idFactor;

	double allreduceCommTimeAccumulator = 0.0;
	struct timespec allreduceStart, allreduceFinish;

	//printf ("I am going to coommStageCounter\n");
	for (commStageCounter = 0; commStageCounter<=numTotalStage; commStageCounter++)
	{
		double allreduceCommTime;
		// calculate parterner id from myid
		int groupSize = pow (2, commStageCounter+1);
		int partnerid = groupSize *  
			(int)floor ( (double)myid / (double) groupSize) +
			(myid%groupSize + groupSize/2) % groupSize;


		//int checkIdx=1;
		//if (commStageCounter == checkIdx &&  myid==1)
		//{
		//	printMatrix_rowMajor (intermediateQRList[checkIdx-1], 2*localNumCols, localNumCols, myid, numprocs);
		//}

		// sequential QR
		if (commStageCounter == 0 )
		{	
			//if (myid == numprocs -1)
			//{
			//	printMatrix_rowMajor (inputMatrixBuffer, localNumRows, localNumCols,myid,numprocs);
			//}

			LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, localNumRows,localNumCols, inputMatrixBuffer, localNumCols, tau);
			if (myid == numprocs -1)
			{
				//printMatrix_rowMajor (inputMatrixBuffer, localNumRows, localNumCols,myid,numprocs);
			}

			if (commStageCounter != numTotalStage)
			{
				saveData_SendBuffer_RestoreR(sendBuffer, intermediateQRList[commStageCounter], inputMatrixBuffer, 
									localNumCols, myid, partnerid);
				//if (myid == 0)
				//	printMatrix_rowMajor (intermediateQRList[0], 2*localNumCols, localNumCols, myid, numprocs);
				//exit(1);
			}
			else
			{
				// restore final R
				saveFinalR_rowOrder (finalR, inputMatrixBuffer);
			}

			// retrive Q from reslut of degeqrf
			LAPACKE_dorgqr(LAPACK_ROW_MAJOR, localNumRows,localNumCols, localNumCols, inputMatrixBuffer, localNumCols, tau);
			
			idFactor=myid < partnerid?0:1;
			//if (myid == numprocs -1)
			//{
			//	//printMatrix_rowMajor (inputMatrixBuffer, localNumRows, localNumCols,myid,numprocs);
			//}
			//printf ("first QR is done\n");

		}
		else //  
		{
			// this part will be replaced by special designed QR for tall R later on
			LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, 2*localNumCols, localNumCols, intermediateQRList[commStageCounter-1], localNumCols, tau);

			//if (myid == numprocs-1)
			//{
			//	//printMatrix_rowMajor (intermediateQRList[commStageCounter-1], 2*localNumCols, localNumCols,myid,numprocs);
			//}

			if (commStageCounter != numTotalStage)
			{
				saveData_SendBuffer_RestoreR(sendBuffer, intermediateQRList[commStageCounter], intermediateQRList[commStageCounter-1], 
									localNumCols, myid, partnerid);
			}
			else
			{
				// restore final R
				saveFinalR_rowOrder (finalR, intermediateQRList[commStageCounter-1]);
			}

			// retrive Q from reslut of degeqrf
			LAPACKE_dorgqr(LAPACK_ROW_MAJOR, 2*localNumCols,localNumCols, localNumCols, intermediateQRList[commStageCounter-1], localNumCols, tau);

			intermediateQListForRetrievingQ[commStageCounter-1]=(intermediateQRList[commStageCounter-1]
										+ idFactor * sizeIntermediateQForRetrievingQ);	

			// if (myid == 0)
			// {
			// //	printMatrix_rowMajor (intermediateQListForRetrievingQ[commStageCounter-1], localNumCols, localNumCols,myid,numprocs);
			// }

			//give intermediateQRMatrix to Qlist after QR
			//if (myid < partnerid)	
			//{
			//	intermediateQListForRetrievingQ[commStageCounter-1]=intermediateQRList[commStageCounter-1];	
			//	//localMatrixCopy(&(intermediateQListForRetrievingQ[commStageCounter-1]) , intermediateQRList[commStageCounter-1], sizeIntermediateQ/2);
			//	if (myid == numprocs-1)
			//	{
			//		//printMatrix_rowMajor (firstPhraseQ, localNumRows, localNumCols,myid,numprocs);
			//		//printMatrix_rowMajor (intermediateQRList[commStageCounter-1], 2*localNumCols, localNumCols,myid,numprocs);
			//		printf("i am fucking in the wrong section, myid:%d,partnerid:%d\n",myid,partnerid);
			//	}
			//}
			//else
			//{
			//	intermediateQListForRetrievingQ[commStageCounter-1]=(intermediateQRList[commStageCounter-1]
			//								+ sizeIntermediateQForRetrievingQ);	
			//	//localMatrixCopy(&(intermediateQListForRetrievingQ[commStageCounter-1]) , intermediateQRList[commStageCounter-1]+(sizeIntermediateQ/2), sizeIntermediateQ/2);
			//}

			idFactor=myid < partnerid ? 0 : 1;
			
		}
		//free(tau);
		//tau = calloc (localNumCols, sizeof(double));

		//clock_gettime(CLOCK_MONOTONIC, &allreduceStart);
		if (commStageCounter != numTotalStage)
		{
			// MPI send and receive communication
			if (myid < partnerid)	
			{
				// send
				MPI_Send ((void*)sendBuffer, sizeSendReceiveBuffer, MPI_DOUBLE, partnerid, 1, MPI_COMM_WORLD);
				// receive
				MPI_Recv ((void*)receiveBuffer, sizeSendReceiveBuffer, MPI_DOUBLE, partnerid, 1, MPI_COMM_WORLD, &status);
			}
			else
			{
				// receive
				MPI_Recv ((void*)receiveBuffer, sizeSendReceiveBuffer, MPI_DOUBLE, partnerid, 1, MPI_COMM_WORLD, &status);
				// send
				MPI_Send ((void*)sendBuffer, sizeSendReceiveBuffer, MPI_DOUBLE, partnerid, 1, MPI_COMM_WORLD);
			}

			// assemble the received intermediate R  intermediateQRMatrix
			restoreDataFromReceiveBufferToIntermediateQRMatrix(intermediateQRList[commStageCounter], receiveBuffer, localNumCols, myid, partnerid);
		}
		//clock_gettime(CLOCK_MONOTONIC, &allreduceFinish);
    	allreduceCommTime = (allreduceFinish.tv_sec - allreduceStart.tv_sec);
    	allreduceCommTime += (allreduceFinish.tv_nsec - allreduceStart.tv_nsec) / 1000000000.0;
    	//if (myid == 0) printf ("%20.18f sec eplapsed\n", elapsed);
    	allreduceCommTimeAccumulator += allreduceCommTime;
	}

	(*elapsedComm)=allreduceCommTimeAccumulator;
	

	// retrive Q from 1st Q and intermediate Q list
	matrixChainMM_homo(intermediateQListForRetrievingQ, numTotalStage, QChainMMResult, localNumCols, localNumCols);

	if (numTotalStage == 0)
	{
		//copy firstPhraseQ to finalQ
		int firstPhraseQidx;
		for (firstPhraseQidx =0; firstPhraseQidx < localNumRows*localNumCols; firstPhraseQidx++)
		{
			finalQ.data[firstPhraseQidx] = inputMatrixBuffer[firstPhraseQidx]; 
		}

	}
	else if (numTotalStage == 1)
	{
		if (myid == numprocs-1)
		{
			//printMatrix_rowMajor (firstPhraseQ, localNumRows, localNumCols,myid,numprocs);
			//printMatrix_rowMajor (intermediateQListForRetrievingQ[0], localNumCols, localNumCols,myid,numprocs);
		}
		cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, localNumRows, localNumCols,localNumCols
				, 1.0, firstPhraseQ, localNumCols, intermediateQListForRetrievingQ[0], localNumCols
				, 0.0, finalQ.data, localNumCols);	
		if (myid == numprocs-1)
		{
			//printMatrix_rowMajor (finalQ.data, localNumRows, localNumCols,myid,numprocs);
		}
	}
	else
	{
		cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, localNumRows, localNumCols,localNumCols
				, 1.0, firstPhraseQ, localNumCols, QChainMMResult, localNumCols
				, 0.0, finalQ.data, localNumCols);	
	}

	

	// release all calloc buffer
	int freeIdx;
	for (freeIdx=0;freeIdx<numTotalStage;freeIdx++)
	{
		free (intermediateQRList[freeIdx]);
	}
	free(tau);
	free(inputMatrixBuffer);
	free(sendBuffer);
	free(receiveBuffer);
	free(QChainMMResult);
}
