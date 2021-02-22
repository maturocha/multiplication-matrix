/*

Algoritm: Matrix multiplication
Type: MPI
Compile: mpicc mpi_multiplication.c -std=c99 -lm -o mpi
Execute: mpirun -np [nro_process] ./mpi
Autor: Matur

*/

#include <mpi.h> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//Const 
#define SIZE 512 
#define NRO_TESTS 10
#define DEBUG TRUE
#define MAESTRO 0 //Master ID
#define FROM_MASTER 0 //From Master to slave
#define FROM_SLAVE 1 //From Slave to Master
#define STEPS = 5

//Functions
void init_vector();
void print_vector(double vec[],int n);
void init_matrix(double mat[SIZE][SIZE],int n);
void print_matrix(double mat[SIZE][SIZE],int n);
void master_task(int slaves);
void slave_task();

// Declare Matrix and Vector
double	A[SIZE][SIZE],  // Matrix A
	      B[SIZE][SIZE],  // Matrix B
	      C[SIZE][SIZE],  // Matrix C
        W[SIZE];     // Vector W

int idProceso; // ID Proceso
int cant_procesos; // Cantidad de Procesos


int main (int argc, char *argv[]) 
{
    int error, slaves;

    MPI_Init(&argc, &argv); //Init MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &idProceso); //Get process_id
    MPI_Comm_size(MPI_COMM_WORLD, &cant_procesos); //Get #process
    if (cant_procesos < 2 ) {
       printf("Error. Se necesita al menos dos procesos...\n");
       MPI_Abort(MPI_COMM_WORLD, error);
       exit(1);
    }
    slaves = cant_procesos-1;

    if (idProceso == 0) 
      //master id = 0
          master_task(slaves); //Init master task
    else
    //slaves id != 0
          slave_task(); //Init slave task
    
    MPI_Finalize(); // Finalizar operaciones MPI
   
   return 0;
}

/*--------------------------------------
	VECTOR Functions
*/

void init_vector(){

	for (int i=0; i<SIZE; i++) {
      W[i] = 1;
  } 
  
}

void print_vector(double vec[],int n){
	for (int i = 0; i < n; ++i)
	{
		printf("%6.2f  ",vec[i] );
	}
}

/*--------------------------------------
	MATRIX Functions
*/


void init_matrix(double mat[SIZE][SIZE],int n){
	
  srand(time(NULL));
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			mat[i][j]= rand() % 10 + 1;
		}
	}
}

void print_matrix(double mat[SIZE][SIZE],int n){

	for (int i = 0; i < n; ++i)
	{
    printf("\n\t\t");
		for (int j = 0; j < n; ++j)
		{
			printf("%6.2f  ",mat[i][j]);
		}
	}
  printf("\n\n");

}

/*--------------------------------------
	Master TASK
*/

void master_task(int slaves) {

  //Vars
  double time_init, time_final; 
  int step; 
  int id_slave, rows, offset, rest, i, j;

  //Vars MPI
  MPI_Status status; // MPI_Recv Status
  MPI_Datatype dt_aux, dt_column; //Special type of data to rcv columns


  offset = 0; //Offset var in rows and columns
  rows = (SIZE / slaves); // Part of Matrix to process
  rest = SIZE % slaves; // Last part of Matrix
  if (rest == 0)
    rest = rows;
  else
    rest = rest + rows;

  for (int test=1;test<=NRO_TESTS;test++) {
    printf("\nInicio de la test= %i\n", test);
    //Init structures
    init_vector();
    print_vector(W,SIZE);
    init_matrix(A,SIZE);
    print_matrix(A,SIZE);
    init_matrix(B,SIZE);
    print_matrix(B,SIZE);


   step=1; //Init step
   time_init = MPI_Wtime();
   for (step= 1; step <= STEPS; step++) 
   {
     printf("\nSteP: %i\n", step);

     //Send task to slaves
     for (id_slave= 1; id_slave <= slaves; id_slave++) {           
       //Last slave
       if (id_slave==slaves)
         rows = rest;
      //Send MPI
       MPI_Send(&step, 1, MPI_INT, id_slave, FROM_MASTER,MPI_COMM_WORLD); //Send step num.
       MPI_Send(&rows, 1, MPI_INT, id_slave, FROM_MASTER,MPI_COMM_WORLD); //Send rows num.
       MPI_Send(&offset, 1, MPI_INT, id_slave, FROM_MASTER,MPI_COMM_WORLD); //Send offset
       //Depends on the step -> the action to do
       switch (step) {
          case 1:
          case 3:
          case 5:
               MPI_Send(&A[offset][0], rows*SIZE, MPI_DOUBLE, id_slave, FROM_MASTER,MPI_COMM_WORLD); //Send part of rows of A
               MPI_Send(&B, SIZE*SIZE, MPI_DOUBLE, id_slave, FROM_MASTER,MPI_COMM_WORLD); //Send B
               MPI_Send(&W, SIZE, MPI_DOUBLE, id_slave, FROM_MASTER,MPI_COMM_WORLD); // Send W
          break; 
          case 2: 
               MPI_Send(&C[offset][0], rows*SIZE, MPI_DOUBLE, id_slave, FROM_MASTER,MPI_COMM_WORLD); //Send part of rows of vector C
          break;
          case 4: 
               MPI_Type_vector(SIZE, 1, SIZE, MPI_DOUBLE, &dt_aux); 
               MPI_Type_create_resized(dt_aux, 0, sizeof(double), &dt_column);
               MPI_Type_commit(&dt_column);
               MPI_Send(&C[0][offset], rows, dt_column, id_slave, 0, MPI_COMM_WORLD);  //Send part of columns of vector C
          break;
          }  
       offset = offset+ rows;          
     }

    //Recive task to slaves
     for (id_slave= 1; id_slave <= slaves; id_slave++)
     {   
        MPI_Recv(&offset  , 1, MPI_INT, id_slave, FROM_SLAVE, MPI_COMM_WORLD, &status); //Rcv offset
        MPI_Recv(&rows, 1, MPI_INT, id_slave, FROM_SLAVE, MPI_COMM_WORLD, &status); //Rcv rows num
        switch (step) {
          case 1:
          case 3:
          case 5:
               MPI_Recv(&C[offset][0], rows*SIZE, MPI_DOUBLE, id_slave, FROM_SLAVE, MPI_COMM_WORLD, &status); 
          break; 
          case 2: 
               MPI_Recv(&W[offset], rows, MPI_DOUBLE, id_slave, FROM_SLAVE, MPI_COMM_WORLD, &status); 
          break;
          case 4:
               MPI_Recv(&W[offset], rows, MPI_DOUBLE, id_slave, FROM_SLAVE, MPI_COMM_WORLD, &status); 
          break; 
        }  
      }
     //print_matrix(C,SIZE);
     //print_vector(W,SIZE);
     
      // Finalize slaves
      if ((test==NRO_TESTS) && (step = 5)) {

        for (id_slave= 1; id_slave<=slaves;id_slave++)
          MPI_Send(&step, 1, MPI_INT, id_slave, FROM_MASTER,MPI_COMM_WORLD);      
      }          
      rows = (SIZE / slaves); 
      offset = 0;  
  }
 
  time_final = MPI_Wtime();
  //print_matrix(C,SIZE);
  //print_vector(W,SIZE);
  printf("\nTime: %f\n\n", time_final - time_init);
  }
}

/*--------------------------------------
	Slave TASK
*/

void slave_task(){
  int step=1, rows, offset, i, j, k;
  MPI_Status status; // MPI_Recv Status
  MPI_Datatype dt_column,dt_aux;  //Special type of data to rcv columns
  double ac;

  for (step= 1; step <= STEPS; step++) { 
      MPI_Recv(&step, 1, MPI_INT, MAESTRO, FROM_MASTER, MPI_COMM_WORLD, &status); //Rcv step nro
      switch (step) {
        case 1: 
        case 3: 
        case 5: 
          MPI_Recv(&rows, 1, MPI_INT, MAESTRO, FROM_MASTER, MPI_COMM_WORLD, &status); 
          MPI_Recv(&offset  , 1, MPI_INT, MAESTRO, FROM_MASTER, MPI_COMM_WORLD, &status); 
          MPI_Recv(&A,rows*SIZE, MPI_DOUBLE, MAESTRO, FROM_MASTER, MPI_COMM_WORLD, &status); 
          MPI_Recv(&B,SIZE*SIZE, MPI_DOUBLE, MAESTRO, FROM_MASTER, MPI_COMM_WORLD, &status); 
          MPI_Recv(&W,SIZE, MPI_DOUBLE, MAESTRO, FROM_MASTER, MPI_COMM_WORLD, &status);

          //C(i,j) = ∑ √(A(i,k) − W(k))**2 * (B(k,j) − W(k))**2 
          for (i=0; i<rows; i++)  
            for (j=0; j<SIZE; j++)
              {
                ac = 0;
                for (k=0; k<SIZE; k++)
                  ac += sqrt( pow(A[i][k]-W[k],2) * pow(B[k][j]-W[k],2) );
                C[i][j]= ac;
              }
          MPI_Send(&offset, 1, MPI_INT, MAESTRO, FROM_SLAVE, MPI_COMM_WORLD); //Send offset to master
          MPI_Send(&rows, 1, MPI_INT, MAESTRO, FROM_SLAVE, MPI_COMM_WORLD);  //Send rows number
          MPI_Send(&C, rows*SIZE, MPI_DOUBLE, MAESTRO, FROM_SLAVE, MPI_COMM_WORLD); //Send part of C
        break;
        case 2: 
           MPI_Recv(&rows, 1, MPI_INT, MAESTRO, FROM_MASTER, MPI_COMM_WORLD, &status); //Rcv rows numbers
           MPI_Recv(&offset  , 1, MPI_INT, MAESTRO, FROM_MASTER, MPI_COMM_WORLD, &status); //Rcv offset
           MPI_Recv(&C,rows*SIZE, MPI_DOUBLE, MAESTRO, FROM_MASTER, MPI_COMM_WORLD, &status); //Rcv part of C
           //W[i] = Prom row i of C
           for (i=0; i<rows; i++) {   
               ac = 0;
               for (j=0; j<SIZE; j++) {
                    ac = ac + C[i][j];
                 }
               W[i] = ac / SIZE;
             }     
           MPI_Send(&offset, 1, MPI_INT, MAESTRO, FROM_SLAVE, MPI_COMM_WORLD);
           MPI_Send(&rows, 1, MPI_INT, MAESTRO, FROM_SLAVE, MPI_COMM_WORLD);  
           MPI_Send(&W, rows, MPI_DOUBLE, MAESTRO, FROM_SLAVE, MPI_COMM_WORLD);
        break;
        case 4:
           MPI_Recv(&rows, 1, MPI_INT, MAESTRO, FROM_MASTER, MPI_COMM_WORLD, &status); 
           MPI_Recv(&offset  , 1, MPI_INT, MAESTRO, FROM_MASTER, MPI_COMM_WORLD, &status); 
           MPI_Type_vector(SIZE, 1, SIZE, MPI_DOUBLE, &dt_aux); 
           MPI_Type_create_resized(dt_aux, 0, sizeof(double), &dt_column);
           MPI_Type_commit(&dt_column);
           MPI_Recv(&C,rows, dt_column, MAESTRO, FROM_MASTER, MPI_COMM_WORLD, &status);
           //W[i] = Prom column i of C
           for (j=0; j<rows; j++) {   
             ac = 0;
             for (i=0; i<SIZE; i++){
                ac = ac + C[i][j];
                 }
             W[j] = ac / SIZE;
             } 
           MPI_Send(&offset, 1, MPI_INT, MAESTRO, FROM_SLAVE, MPI_COMM_WORLD); 
           MPI_Send(&rows, 1, MPI_INT, MAESTRO, FROM_SLAVE, MPI_COMM_WORLD);  
           MPI_Send(&W, rows, MPI_DOUBLE, MAESTRO, FROM_SLAVE, MPI_COMM_WORLD);
        break;
      }              
   }
}