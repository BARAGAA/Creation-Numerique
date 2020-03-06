#include "matrix.hpp"
#include "matrix.hpp"
#include <cstdint>
#include <iostream>
#define taillBloc n/2,n/2
#define nouvelleMatrice allocateMatrix(taillBloc)
 using namespace std;

/*
    Nous allons stocker des matrices dans des tablaux à 1 dimension:
    Par exemple, matrice de taille n x m est stockué comme:
        A_tab = [A_11, A_12, ... A_1m, A_21, ... A_2m, ... , A_1n, ..., A_nm]

    Par conséquant, un élement A_ij aurait quelle indice dans le tableau A_tab ? 
    answer (i * j)-1

 */

/* Memory allocation for a matrix of size n x m and initilization to 0  */
double *allocateMatrix(uint64_t n,uint64_t m) {
  double *A;
  A = (double *) calloc (n * m, sizeof(double));
  return A;
}


/* Frees the memory allocated to matrix A
*/
void freeMatrix(double *A) {
    free(A);
}

/* Allocates a n sized vector and initializes all entries to 0 
*/
double *allocateVector(uint64_t n) {
  double *v; 
  v = (double *) calloc(n, sizeof(double));
  return v;
}

/* Trees the memory allocated to a vector
*/
void freeVector(double *v) {
  free(v);
}


/* Sets a n * m matrix A to all zeros */
void setMatrixZero(double *A, uint64_t n, uint64_t m) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
        /* Note that for a n x m matrix flattened to a 1D array, 
        element A_ij has index i * m + j
        */
      A[i * m + j] = 0.0; 
    }
  }
}

/* Sets a n * n matrix A to identity */
void setMatrixIdentity(double *A, uint64_t n) {
  uint64_t i, j;

  for (i=1;i<n;i++) {
    for (j=1;j<n;j++) {
     A[i * n + j] = 0.0;
    }
    A[i * n + i] = 1.0;
  }
}



/* Copies a matrix  
*/
void copyMatrix(double *B, double *A, uint64_t n, uint64_t m) {
  uint64_t i,j;

  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      B[i * m + j] = A[i * m + j]; 
    }
  }
}

/*
Writes a matrix to a stream. For example, writing a matrix to standard output is
writeMatrix(stdout, A, n, m);
A sream can also be a file. 
*/
void writeMatrix(FILE *stream, double *A, uint64_t n, uint64_t m)
{
	fprintf(stream, "%d %d \n", (int)n, (int)m);
	int i, j;
	for(i = 0; i < n; ++i)
	{
	      for(j = 0; j < m; ++j)
	      {
		      fprintf(stream, "%f \t", A[i * m + j]);
	      }
	      fprintf(stream, "\n");
	}
}



//The function computes the element-by-element abs of matrix A
void absMatrix(double *Aabs,double *A, uint64_t n, uint64_t m)
{
	uint64_t i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
            Aabs[i*m + j] = fabs(Aabs[i*m + j]);
		}
	}

}


/*
Performs addition of two matrix A (size n x m) and B (size n x m).
The result S = A + B is a n x m matrix.
We consider that S is allocated outside the function.
*/
void matrixAdd(double *S, double *A, double *B, uint64_t n, uint64_t m){
    uint64_t i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
            S[i*m + j] = A[i*m + j] + B[i*m + j];
		}
	}
}

/*
Performs subtraction of two matrix A (size n x m) and B (size n x m).
The result S = A - B is a n x m matrix.
We consider that S is allocated outside the function.
*/
void matrixSub(double *S, double *A, double *B, uint64_t n, uint64_t m){
    uint64_t i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
            S[i*m + j] = A[i*m + j] - B[i*m + j];
		}
	}
}



/* For a double m x n matrix A the function returns its maximum in absolute value
element. */
double getMaxInMatrix(double max, double *A, uint64_t n, uint64_t m)
{
	double maxA = fabs(A[0]);
	double current = fabs(A[0]);
	int i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			current = fabs(A[i * m + j]);
			if(current > maxA)
				maxA = current;
		}
	}
    return maxA;

}


/* Rajouter les prototypes de vos méthodes ici. Par exemple */

/* Performs naive multiplication of matrix A (size p x k) by a matrix B (size k x r).
The result matrix S = A*B  is of size (k x r).
We assume that S has already been allocated outside the function.
*/
void        matrixMultiplyNaive (double *S, double *A, double *B, uint64_t p, uint64_t k, uint64_t r){ 

  setMatrixZero(S,k,r); // si elle n' est pas déja au zéro
  for(int i = 0 ; i < p ; i++ ){
    for(int j = 0 ; j < r ; j++ ){
      for(int l = 0 ; l < k ; l++ ){
          S[(i *r) + j ] += A[ ( i * k ) + l ] * B[ (l * r) + j] ;
        }
      }
    }
    
}
bool ValidN(uint64_t n){
  if(n>0)
    {
        while(n%2 == 0)
        {
            n = n/2;
        }
        return (n==1);
    }
}


/* Performs a multiplication of two sqaure matrices A and B (size n x n) by Strassen algorithm.
    We assume that S has already been allocated outside the function.
*/
void        matrixMultiplyStrassen (double *S, double *A, double *B, uint64_t n){
    setMatrixZero(S,n,n); // mettre à zero la matrice de resultat au cas où
    if(n==1){
      S[0] = A[0] * B[0] ;
    }
    else if(n == 2){         // si la matrice est en forme 2x2 alors l'algorithme strassen de base est fait
      double st1 = (A[0] + A[3]) * (B[0] + B[3]);
      double st2 = (A[2] + A[3]) * B[0];
      double st3 = (B[1] - B[3]) * A[0];
      double st4 = (B[2] - B[0]) * A[3];
      double st5 = (A[0] + A[1]) * B[3];
      double st6 = (A[2] - A[0]) * (B[0]+ B[1]);
      double st7 = (A[1] - A[3]) * (B[2]+ B[3]);
      // sommation de Strassen (7 multiplications et 8 additions)
      S[0] = st1 + st4 - st5 + st7;
      S[1] = st3 + st5;
      S[2] = st2 + st4;
      S[3] = st1 - st2 + st3 + st6;

    }else if (ValidN(n)){     //sinon la matrice est coupé vers  autres matrices plus petits  si le taille est une puissance de 2

      //Simulations des valeures dans la 1er brance de si 
      double* A0 = nouvelleMatrice;
      double* A1 = nouvelleMatrice;
      double* A2 = nouvelleMatrice;
      double* A3 = nouvelleMatrice;
      double* B0 = nouvelleMatrice;
      double* B1 = nouvelleMatrice;
      double* B2 = nouvelleMatrice;
      double* B3 = nouvelleMatrice;

      //couper les matrices vers  des sous matrices plus petits et les remplir
      for(int i = 0 ;i <= n-1; i++){
       if(i < n/2){ 
         for(int j = 0 ; j <= n-1 ; j++){
           if(j< n / 2){
              A0[(i * (n/2)) + j] = A[(i * n) + j ];

              B0[(i * (n/2)) + j] = B[(i * n) + j ];
           }else{
              A1[(i * (n/2)) + (j % (n/2))] = A[(i * n) + j ];

              B1[(i * (n/2)) + (j % (n/2))] = B[(i * n) + j ];
           }
         }
       }else{
         for(int j = 0 ; j <= n-1 ; j++){
           if(j< n / 2){
              A2[((i%2) * (n/2)) + j] = A[(i * n) + j ];
              B2[((i%2) * (n/2)) + j] = B[(i * n) + j ];
           }else{
              A3[((i%2) * (n/2)) + (j % (n/2))] = A[(i * n) + j ];
              B3[((i%2) * (n/2)) + (j % (n/2))] = B[(i * n) + j ];
           }
         }
       }
      }

      // matrices auxilieres pour les opérations  en simulations de ce qui se passe dans le 1er brance de si 
      double* aux1 = nouvelleMatrice;
      double* aux2 = nouvelleMatrice;
       
      double* st1 = nouvelleMatrice;
      double* st2 = nouvelleMatrice;
      double* st3 = nouvelleMatrice;
      double* st4 = nouvelleMatrice;
      double* st5 = nouvelleMatrice;
      double* st6 = nouvelleMatrice;
      double* st7 = nouvelleMatrice;
      //st1
      matrixAdd(aux1 , A0 , A3 ,taillBloc);
      matrixAdd(aux2 , B0 , B3 ,taillBloc);
      matrixMultiplyStrassen(st1, aux1, aux2, n/2);

      //st2
      matrixAdd(aux1 , A2 , A3 ,taillBloc);
      matrixMultiplyStrassen(st2, aux1, B0, n/2);

      //st3
      matrixSub(aux1 , B1 , B3 ,taillBloc);
      matrixMultiplyStrassen(st3, A0, aux1, n/2);

      //st4
      matrixSub(aux1 , B2 , B0 ,taillBloc);
      matrixMultiplyStrassen(st4,A3,aux1,n/2);

      //st5
      matrixAdd(aux1 , A0 , A1 ,taillBloc);
      matrixMultiplyStrassen(st5, aux1, B3, n/2);

      //st6
      matrixSub(aux1 , A2 , A0 ,taillBloc);
      matrixAdd(aux2 , B0 , B1 ,taillBloc);
      matrixMultiplyStrassen(st6, aux1, aux2, n/2);

      //st7
      matrixSub(aux1 , A1 , A3 ,taillBloc);
      matrixAdd(aux2 , B2 , B3 ,taillBloc);
      matrixMultiplyStrassen(st7, aux1, aux2, n/2);
       // allocation des blockes 

      double* bloc1 = nouvelleMatrice;
      double* bloc2 = nouvelleMatrice;
      double* bloc3 = nouvelleMatrice;
      double* bloc4 = nouvelleMatrice;

      // sommation des matrices an tant que des elements dans la branche si

      //bloc 1 
      matrixAdd(aux1,st1,st4, taillBloc);
      matrixSub(aux2,st7,st5, taillBloc);
      matrixAdd(bloc1, aux1, aux2, taillBloc);
      //bloc 2
      matrixAdd(bloc2,st3,st5, taillBloc);

      //bloc 3
      matrixAdd(bloc3,st2, st4, taillBloc);

      //bloc4
      matrixAdd(aux1,st3,st6, taillBloc);
      matrixSub(aux2,st1,st2, taillBloc);
      matrixAdd(bloc4 ,aux1 ,aux2 ,taillBloc);

      //Remplir la matrice S par les blocs 
      for(int i = 0 ;i <= n-1; i++){
       if(i < n/2){
         for(int j = 0 ; j <= n-1 ; j++){
           if(j< n / 2){
             S[(i * n) + j ] =  bloc1[(i * (n/2)) + j] ;
           }else{
             S[(i * n) + j ] = bloc2[(i * (n/2)) + (j % (n/2))] ;
           }
         }
       }else{
         for(int j = 0 ; j <= n-1 ; j++){
           if(j< n / 2){
             S[(i * n) + j ] = bloc3[((i%(2)) * (n/2)) + j];
           }else{
             S[(i * n) + j ] = bloc4[((i%(2))* (n/2)) + (j % (n/2))];
           }
         }
       }
       /*désallocation des tableaux dynamics utilisés dans la processus*/
      free(A0);
      free(A1);
      free(A2);
      free(A3);
      free(B0);
      free(B1);
      free(B2);
      free(B3);
      free(st1);
      free(st2);
      free(st3);
      free(st4);
      free(st5);
      free(st6);
      free(st7);
      free(aux1);
      free(aux2);
      free(bloc1);
      free(bloc2);
      free(bloc3);
      free(bloc4);
      }
    }else
      cout<< "taille non possible"<<endl;  
}

/* 
    Solves a system of linear equations Ax=b for a double-precision matrix A (size n x n).
    Uses iterative ascension algorithm. 
    After the procedure, x contains the solution of Ax=b.
    We assume that x has been allocated outside the function.
*/


void        SolveTriangularSystemUP   (double *x, double *A, double *b, uint64_t n){
  for(int i = n-1 ; i >= 0 ; i--){ //l'indice de colonne en ordre décroissant
      // standardiser  à 1  le pivot et modifier la valeure correspondant dans le vecteur  
        b[i] = b[i] / A[(i * n) + i]; // le pivot
        x[i] = b[i];  
        A[(i * n) + i] = 1 ; // le pivot  =1 
      for(int j = i-1 ; j >= 0; j--){ //l'indice du ligne commencant par l'un haut de pivot 
        b[j] -=(A[j*n+i]*b[i]); //modification de l'element correspondant dans le vecteur de resultat
        A[(j * n) +i] -= (A[(j * n) + i]); // performer l'opération sur l'element haut de pivot dans la matrice

      }
    }
  }
 
    /*Performs Gauss elimination for given a matrix A (size n x n) and a vector b (size n).
    Modifies directly matrix A and vector b.
    In the end of the procedure, A is upper truangular and b is modified accordingly.
    Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is impossible to triangularize. 
*/

bool        Triangularize           (double *A, double *b, uint64_t n){
  int pivotSuiv = n+1; // passer à la pivot suivant dans une structrue 1D  de matrice
  for (int i =0 ;i < n*n ; i += pivotSuiv){ // utiliser la variable pivot suiv pour passer à la pivot suivant
    if(A[i] == 0){            // un pivot mort (=0) 
      int j;
      for( j=i  ; A[j] == 0 && j< n*n ;j +=n  ){/*Rien*/}// incrementation of j by n to pass to a line that have a non-null pivot if there is any 
      
      if(j>n*n){
        return false; // aucun pivot vivant (!=0) existe donc la matrice n'est pas inversible.
      }
      for(int s = j%n ; s < n ; s++){     //echanger le ligne de pivot mort avec l'autre trouvé avec un pivot vivant
      /*element par element*/ 
        double aux = A[i + s] ;
        A[i + s] = A[j + s];
        A[j + s] = aux;
      }
      double aux = b[i/n]; // echangement de la valeur qui corresponde dans le vecteur
      b[i/n] = b[j/n];
      b[j/n] = aux;
    }
    for(int l =i+n ; l <(n*n) ; l+=n){
        double sousPivot =A[l]; //les elems sous le pivot
        int lineEnd = n - (l%n); //dernier dans la ligne
        for (int s = 0 ; s <lineEnd ; s++){
          A[l+s] -=(sousPivot/A[i]*A[i+s]); //Echangement des valeures dessous le pivot  pour fourmer la matrice (upper treiangle).
        }
        b[l/n] -=(sousPivot/A[i])*b[i/n];
    }
  }
  return true;
}
  

/*
    Solves a system of linear equations Ax=b, given a matrix A (size n x n) and vector b(size n).
    Uses Gauss elimination algorithm based on truangularization and the ascension solving.
    After the procedure, vector x contains the solution to Ax=b.
    We assume that x has been allocated outside the function.
        Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is of rank <n .
*/
bool        SolveSystemGauss        (double *x, double *A, double *b, uint64_t n){
    bool inversible = Triangularize(A,b,n); // triangulise la matrice et vérifie si elle est incversible ( triangulisé avec succes) 
    if(inversible)
      SolveTriangularSystemUP(x,A,b,n); // resoudre par l'élimination de gausse

    return inversible;
}



bool LU_Decomposition(double *A, uint64_t n){
  for (int i =0 ;i < n*n ; i +=(n+1)){ //  passer à la pivot suivant
      if(A[i] == 0)
        return false;

      for(int l =i+n ; l <(n*n) ; l+=n){
        double sousPivot =A[l]; //les elems sous le pivot
        int lineEnd = n - (l%n); //dernier dans la ligne
        for (int s = 0 ; s <lineEnd ; s++){
          A[l+s] -=(sousPivot/A[i]*A[i+s]); //Echangement des valeures dessous le pivot  pour fourmer la matrice (upper treiangle).
        }
    }
  }
  return true;
}

double det(double *A,uint64_t n){
  if (LU_Decomposition(A,n))
  {
  double d =1;
  for(int i =0; i<= n ;i+= (n + 1))
    d *= A[i];
  return d;
}
return 0;

