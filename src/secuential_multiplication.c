/*

Algoritm: Matrix multiplication
Type: Secuential
Compile: gcc secuential_multiplication.c -std=c99 -lm -o secuential
Execute: ./secuential
Autor: Matur

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//Const 
#define SIZE 512 
#define NRO_TESTS 10
#define DEBUG TRUE

//Functions
void init_vector();
void print_vector(double vec[],int n);
void init_matrix(double mat[SIZE][SIZE],int n);
void print_matrix(double mat[SIZE][SIZE],int n);
void multiplication();

// Declare Matrix and Vector
double	A[SIZE][SIZE],  // Matrix A
	      B[SIZE][SIZE],  // Matrix B
	      C[SIZE][SIZE],  // Matrix C
        W[SIZE];     // Vector W

int main(int argc,  char *argv[]) {

  int i, j, k, fase;
  float ac_prom_w;
  clock_t time_init, time_final;

  //for (int test=1; test<=NRO_TESTS;test++) {
    //printf("\nTest Nro= %i\n", test);
    init_vector();
    print_vector(W,SIZE);
    init_matrix(A,SIZE);
    print_matrix(A,SIZE);
    init_matrix(B,SIZE);
    print_matrix(B,SIZE);

    time_init = clock();

    /*
      FASE 1
    */
    //Step 1
    multiplication();
    print_matrix(C,SIZE);

    /*
      FASE 2
    */
   //Step 2
   //W[i] = Prom row i of C
    for(int i=0;i<SIZE;i++){
        ac_prom_w = 0;
        for(int j=0;j<SIZE;j++){
          ac_prom_w += C[i][j];
        }
        W[i] = ac_prom_w / SIZE;
    }

    print_vector(W,SIZE);
    //Step 3
    multiplication();
    print_matrix(C,SIZE);

    /*
      FASE 3
    */
    //Step 4
     //W[i] = Prom column i of C
    for(int j=0;j<SIZE;j++){
        ac_prom_w = 0;
        for(int i=0;i<SIZE;i++){
          ac_prom_w += C[i][j];
        }
        W[j] = ac_prom_w / SIZE;
    }

    print_vector(W,SIZE);
    //Step 5
    multiplication();

    //FINAL

    time_final = clock();
    printf("Time: %f\n", (double)(time_final - time_init) / CLOCKS_PER_SEC);

  
 //}

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

void multiplication()
{

  float ac = 0;
  float term1 = 0;
  float term2 = 0;

  //C(i,j) = ∑ √(A(i,k) − W(k))**2 * (B(k,j) − W(k))**2 
  
  for(int i=0;i<SIZE;i++){
      for(int j=0;j<SIZE;j++){
        for(int k=0;k<SIZE;k++){

          term1 = pow( (A[i][k] - W[k]) ,2);
          term2 = pow( (B[k][j] - W[k]) ,2);
          ac += sqrt( term1*term2 );

        }

        C[i][j] = ac;
        ac = 0;

      }
  }

}