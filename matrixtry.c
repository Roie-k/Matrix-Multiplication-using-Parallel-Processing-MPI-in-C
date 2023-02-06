#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define NRA 62 /// Number of rows in matrix A
#define NCA 15/// number of columns in matrix A
#define NCB 7/// Number of columns in matrix B
#define MASTER 0/// Task id of first task
#define FROM_MASTER 1 ///setting a message type
#define FROM_WORKER 2 ///setting a message type
#define NRB NRA ///Number of rows in B = Number of columns

void printMatrixA(double matrix[][NCA],int nrows,int ncolumns, char* title) {
    /* Print results */
    printf("******************************************************\n");
    printf("%s:\n",title);
    for (int i=0; i<nrows; i++) {
        printf("\n");
        for (int j=0; j<ncolumns; j++) {
			printf("%6.2f   ", matrix[i][j]);
        }
    }
    printf("\n******************************************************\n");
}

void printMatrixB(double matrix[][NCB],int nrows,int ncolumns, char* title) {
    /* Print results */
    printf("******************************************************\n");
    printf("%s:\n",title);
    for (int i=0; i<nrows; i++) {
        printf("\n");
        for (int j=0; j<ncolumns; j++) {
			printf("%6.2f   ", matrix[i][j]);
        }
    }
    printf("\n******************************************************\n");
}

void printMatrixC(double matrix[][NCB],int nrows,int ncolumns, char* title) {
    /* Print results */
    printf("******************************************************\n");
    printf("%s:\n",title);
    for (int i=0; i<nrows; i++) {
        printf("\n");
        for (int j=0; j<ncolumns; j++) {
			printf("%6.2f   ", matrix[i][j]);
        }
    }
    printf("\n******************************************************\n");
}



int main (int argc, char *argv[]) {
	int nproc; /// number of tasks in partition
	int id; /// a task identifier
	int numslaves; /// numver of worker tasks
	int source; /// task id of message source
	int dest; /// task id of message destination
	int mtype; /// message type
	int rows; /// rows of matrix A sent to each worker
	int averow; /// used to determine rows sent to each worker
	int extra; /// used to determine rows sent to each worker 
	int offset; /// used to determine rows sent to each worker
	int i; ///misc
	int j; ///misc
	int k; ///misc 
	int rc; ///misc
	double a[NRA][NCA]; /// matrix A to be multiplied
	double b[NRB][NCB]; /// matrix B to be multiplied
	double c[NRA][NCB]; /// result matrix C
   
    MPI_Status status;


    //Seed for rand function
    srand(time(NULL));

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&id);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);

	if (nproc < 2) {
		printf("Need at least two MPI tasks. Quitting ...\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		exit(1);
		}

    /* Number of slave = number of processes - 1 */
    numslaves = nproc-1;

    /// Master Initialization
    if (id == MASTER) {
        printf("MPI_MM has started with %d tasks, & %d slaves.\n",nproc,numslaves);
        printf("Initializing arrays...\n");
        for (i=0; i<NRA; i++) {
            for (j=0; j<NCA; j++) {
                a[i][j]= i+j;
            }
        }
        for (i=0; i<NRB; i++) {
            for (j=0; j<NCB; j++) {
                b[i][j]= i*j;
            }
        }

        printMatrixA(a,NRA,NCA,"Matrix A");
        printMatrixB(b,NRB,NCB,"Matrix B");

        /* Comute start time */
        double start = MPI_Wtime();

        /// Send matrix data to the worker tasks
        averow = NRA/numslaves;
        extra = NRA%numslaves;/* Marker to identify which processes will receive extra 					lines */
        offset = 0; /* serves to mark what has already been sent */
        mtype = FROM_MASTER;
        for (dest=1; dest<=numslaves; dest++) {
            rows = (dest <= extra) ? averow+1 : averow;
            printf("Sending %d rows to task %d offset=%d\n",rows,dest,offset);
            /* Send the offset value */
            MPI_Send(&offset, 1, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
            /* Send the numbers of lines that will be sent */
            MPI_Send(&rows, 1, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
            /* Send the lines of matrix A */
            MPI_Send(&a[offset][0], rows*NCA, MPI_DOUBLE, dest, FROM_MASTER,
                     MPI_COMM_WORLD);
            /* Send matrix B */
            MPI_Send(&b, NRA*NCB, MPI_DOUBLE, dest, FROM_MASTER, MPI_COMM_WORLD);
            offset = offset + rows;
        }

        /// Receive results from worker tasks
        mtype = FROM_WORKER;
        for (i=1; i<=numslaves; i++) {
            source = i;
            /* Receives offset value and number of lines to be able to place the values 		in matrix C */
            MPI_Recv(&offset, 1, MPI_INT, i, FROM_WORKER, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, i, FROM_WORKER, MPI_COMM_WORLD, &status);
            MPI_Recv(&c[offset][0], rows*NCB, MPI_DOUBLE, i, FROM_WORKER,
                     MPI_COMM_WORLD, &status);
            printf("Getting lines from slave %d\n",source);
        }

        printMatrixC(c,NRA,NCB,"Final Result");

        /* Compute finish time and display total process time */
        double finish = MPI_Wtime();
        printf("Completed in %f seconds.\n", finish - start);
    }


    /* Saves Process - ID > 0 */
    if (id > MASTER) {
        mtype = FROM_MASTER;
        MPI_Recv(&offset, 1, MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);
        MPI_Recv(&a, rows*NCA, MPI_DOUBLE, MASTER, FROM_MASTER, MPI_COMM_WORLD, 		&status);
        MPI_Recv(&b, NRB*NCB, MPI_DOUBLE, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);

        for (k=0; k<NCB; k++) {
            for (i=0; i<rows; i++) {
                c[i][k] = 0.0;
                for (j=0; j<NCA; j++) {
                    c[i][k] = c[i][k] + a[i][j] * b[j][k];
                }
            }
        }
        mtype = FROM_WORKER;
        MPI_Send(&offset, 1, MPI_INT, MASTER, FROM_WORKER, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, FROM_WORKER, MPI_COMM_WORLD);
        MPI_Send(&c, rows*NCB, MPI_DOUBLE, MASTER, FROM_WORKER, MPI_COMM_WORLD);
    }
    MPI_Finalize();
}
