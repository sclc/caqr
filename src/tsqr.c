// tsqr main function

#include "tsqr.h"

int main (int argc, char** argv) 
{
	int myid, numprocs;
	int ierr;

        if (argc<4) {printf("Note, we need 3 parameters, 1st: filename, 2nd: numRows, 3rd:numCols\n"); exit(1);}

	ierr = MPI_Init(&argc, &argv);
	//printf ("ierr=%d\n",ierr);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	if (myid == 0) printf ("%d process in total",numprocs);
	//printf ("myid=%d \n",myid);
	//printf ("in mainf numprocs=%d\n",numprocs);
	//binaryAllReduceButterflyPrototype(myid, numprocs);

        int numRows, numCols;

        int filenameLen= strlen(argv[1]);
        int buffersize=100;
        char filenameContainer[buffersize];
        if (filenameLen>buffersize) {printf("too long file name\n");exit(1);}

        numRows = atoi(argv[2]);
        numCols = atoi(argv[3]);

        strcpy (filenameContainer, argv[1]);
	denseType localMatrix;
	dataDist_rowMajor(&localMatrix, numRows, numCols, filenameContainer, myid, numprocs);
        //if(myid == 2) {printMatrix (localMatrix, myid, numprocs);}
        //printMatrix (localMatrix, myid, numprocs);
	//printf ("myid:%d, start from Rows:%d\n", myid, localMatrix.start_idx);

	//call tsqr allreduction routines
	int expCounter;
	int totalExp = 10;
	double elapsedTimeAccmulator = 0.0;
	double elapsedCommTimeAccmulator=0.0;

	struct timespec start, finish;
	for (expCounter=0;expCounter<totalExp;expCounter++)
	{
		double elapsed;
		double elapsedComm;
		denseType finalQ;
		denseType finalR;
	
		initDenseType_rowOrder (&finalQ, localMatrix.local_num_row, localMatrix.local_num_col
				, localMatrix.global_num_row, localMatrix.global_num_col
				, localMatrix.start_idx);
	
		initDenseType_rowOrder (&finalR, localMatrix.local_num_col, localMatrix.local_num_col
				, localMatrix.local_num_col, localMatrix.local_num_col
				, 0);
		
		clock_gettime(CLOCK_MONOTONIC, &start);
	
		tsqrAllReduction_rowOrder_butterFly_raw (localMatrix, finalQ, finalR, myid, numprocs, &elapsedComm);
		MPI_Barrier(MPI_COMM_WORLD);
	
		clock_gettime(CLOCK_MONOTONIC, &finish);
    	elapsed = (finish.tv_sec - start.tv_sec);
    	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    	//printf ("%f sec eplapsed\n", elapsed);
    	elapsedTimeAccmulator += elapsed;
    	elapsedCommTimeAccmulator += elapsedComm;

	
		//double resDistance = tsqrResCheck ( localMatrix, finalQ, finalR);
		//printf ("myid: %d, resDistance=%f\n",myid,resDistance );
	
		// if (myid == 0)
		// {
		// 	//printMatrix(finalR, myid, numprocs);
		// 	//printMatrix(finalQ, myid, numprocs);
		// }

		free (finalQ.data);
		free (finalR.data);
	}

	elapsedTimeAccmulator/=totalExp;
	elapsedCommTimeAccmulator /= totalExp;

	if (myid == 0) 
	{
		printf ("averaged time cost for tsqrAllReduction_rowOrder_butterFly_raw: %15.12f\n", elapsedTimeAccmulator);
		printf ("averaged time cost for communicaiton in tsqrAllReduction_rowOrder_butterFly_raw: %15.12f\n", elapsedCommTimeAccmulator);
	}

	ierr = MPI_Finalize();
}

