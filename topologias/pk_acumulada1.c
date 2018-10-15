///////////////////////////////////////////////////////////
//// Programa para trasnformar una distribucion p(k)  ////
//// (leyendo del archivo) en la distribucion acumulada ////
///////////////////////////////////////////////////////////


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 1000
#define alfa 0.0
#define filas 100
#define b 2.500



int J,I,i,c1;
double   M[filas+1][4], pk_acum[filas+1],c2,norma[filas+1], suma;


int main()
{ 

  FILE *fich2;
  FILE *fich3;

  char file2[256], file3[256];

 
 /////////  LECTURA DEL FICHERO  ////////
 


  
  //sprintf(file3,"pk_aleat.txt");
  //sprintf(file3,"pk_sf.txt");

  sprintf(file3,"TauCvsK0.00-bini%.3lf.dat",b);
  //sprintf(file3,"TCvsK0.00-bini%.3lf.dat",b);
  
  
  //sprintf(file2,"pk_aleat_acum.dat");
  //sprintf(file2,"pk_sf_acum.dat");

  sprintf(file2,"TauCvsK0.00-bini%.3lf_acum.dat",b);
  //  sprintf(file2,"TCvsK0.00-bini%.3lf_acum.dat",b);
  
  fich2=fopen(file2,"wt");
  fclose(fich2);
  
 

  for(J=0;J<=filas;J++)
  {
      for(I=0;I<=3;I++)
      {
	  M[J][I]=0;
      }

      pk_acum[J]=0;
      norma[J]=0;
  }

  fich3=fopen(file3,"r");
  for(J=1;J<=filas;J++)
  {
      fscanf(fich3,"%d %lf\n",&c1,&c2);
      M[c1][2]=c2;
  }

  fclose(fich3);


  for (J=1;J<=filas;J++)
   {
    printf("%d %f \n",J,M[J][2]);
   }




  ////Calculo distribucion acumulada////
  suma=0;
 
  for(i=1;i<=filas;i++)
    {
      suma+=M[i][2];     
    }


      for(J=1;J<=i;J++)
	{
	  suma=suma-M[J][2];
	  pk_acum[J]=suma;	
	}




 
  fich2=fopen(file2,"at");
  for (J=1;J<=filas;J++)
    {
      printf("pk_acum[%d]: %lf \n",J,pk_acum[J]);
      
      if(pk_acum[J]!=0)
	{
	  fprintf(fich2,"%d  %lf\n",J,pk_acum[J]);      
	}
    }
  
  
  fclose (fich2);


 

  
}








  
 
