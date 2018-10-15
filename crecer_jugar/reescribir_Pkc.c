////////////////////////////////////////////////////////////////////
//// Programa para reescribir los archivos              ////////
//// Pkct_10-1-b2.10-e0.99_size1000_1000iter.dat         //////
//// eliminando las filas con Pkc=0 (para evitar        ////////
//////  picos en la representacion gnuplot)         ///////////////
///////////////////////////////////////////////////////////



#include <stdio.h>
#include <math.h>
#include <stdlib.h>


//        Pkct_10-1-b2.10-e0.99_size1000_1000iter.dat




#define b 2.22
#define eps  0.99


#define tauT 10 
#define tauD 1



#define N 1000
#define Niter 10


# define filas 37000     //que contiene el archivo de entrada
# define col 3 


int J,I,c1,c2;
double M[filas][col],c3;




int main()
{ 
  
  FILE *fich2;
  FILE *fich3;

  char file2[256], file3[256];
  
  
  


//archivo que leere 
 
sprintf(file3,"Pkct_%d-%d-b%.2lf-e%.2lf_size%d_%diter.dat",tauT,tauD,b,eps,N,Niter);
  
 

//archivo que creare

 sprintf(file2,"Pkct_%d-%d-b%.2lf-e%.2lf_size%d_%diter_.dat",tauT,tauD,b,eps,N,Niter);
  fich2=fopen(file2,"wt");
  fclose(fich2);
  
  
  
  for(J=0;J<filas;J++)
    {
      for(I=0;I<col;I++)
	{
	  M[J][I]=0.0;
	  
	}    
   //printf ("J=%d",J);
   }


 printf ("vector inicializado a cero\n");


 fich3=fopen(file3,"r");     //leo el fichero y lo guardo en un array

   for(J=0;J<filas;J++)
    {
	fscanf(fich3,"%d   %d   %lf\n",&c1,&c2,&c3);
	M[J][0]=c1;
	M[J][1]=c2;
	M[J][2]=c3;
//printf ("J=%d",J);
    }
  
    fclose(fich3);
  



    printf ("archivo almacenado\n");


 



 
    fich2=fopen(file2,"at");     //imprimo el histograma  pero solo las lineas en las cuales Pkc(t)!=0
    for (J=0;J<filas;J++)
    {     
	if(M[J][2]!=0.000000)
	{
	    fprintf(fich2,"%.0lf  %.0lf  %lf\n",M[J][0],M[J][1],M[J][2]);      
	}			
	

//	printf ("J=%d",J);
	if(M[J][0]!=M[J-1][0])   //la separacion de una linea entre bloques de tiempos distintos
	{
	fprintf(fich2,"  \n"); 	
	}
	
    }    
  fclose (fich2);
  
  





}



 




  
 
