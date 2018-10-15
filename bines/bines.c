///////////////////////////////////////////////////////////
//// Programa para trasnformar una distribucion p(k)  ////
//// (leyendo del archivo) haciendo un binning       ////
///////////////////////////////////////////////////////////


#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define alfa 0.0
#define filas 200
#define b 1.100
#define Nbines 25



int J,I,i,c1,n;
double M[filas+1][2],c2;
int k_max,k_min,aux;

double k_max_double, k_min_double, r;     //para el binning
double pk[Nbines+1], k[Nbines+1],double1,double2, suma;
double norma[Nbines+1];

int main()
{ 
  
  FILE *fich2;
  FILE *fich3;

  char file2[256], file3[256];
  
  
  /////////  LECTURA DEL FICHERO  ////////
  
  

  
  //sprintf(file3,"pk_sf.txt");
  //sprintf(file3,"pk_aleat.txt");
  //sprintf(file3,"TauCvsK0.00-bini%.3lf.dat",b);
   sprintf(file3,"TCvsK0.00-bini%.3lf.dat",b);
  
  
  //sprintf(file2,"pk_sf_binning.dat");
  //sprintf(file2,"pk_aleat_binning.dat");
  //sprintf(file2,"TauCvsK0.00-bini%.3lf_binning25.dat",b);
  sprintf(file2,"TCvsK0.00-bini%.3lf_binning25.dat",b);

  fich2=fopen(file2,"wt");
  fclose(fich2);
  
  
  
  for(J=0;J<=filas;J++)
    {
      for(I=0;I<=3;I++)
	{
	  M[J][I]=0;
	}       
    }


  for(J=0;J<=Nbines;J++)
    {
      pk[J]=0;
      k[J]=0;
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
  
  
  
  
  //calculo k_max y k_min
  
  k_max=0;
  for(i=1;i<=filas;i++)
    {
      if(M[i][2]!=0)
	{	  
	  if(i > k_max)
	   {
	     k_max=i;
	   //printf("encontrado mejor max %d\n", k_max);
	   }
	}
    }
  
  k_min=10000;
  for(i=1;i<=filas;i++)
    {
      if(M[i][2]!=0)
	{
	  if(i < k_min)
	    k_min= i;
	}
      
    }

  
  printf("k_max: %d  k_min: %d\n",k_max,k_min);
  
  
  
  
  
  //calculo r

  k_max_double=k_max;
  k_min_double=k_min;

  double1=Nbines;
  r=pow((k_max_double/k_min_double),1.0/double1);    //al loro con las potencias!! pow(base,exponente)

  
 printf("r: %lf\n",r);

 for(n=0;n<=Nbines;n++)     //calculo los intervalos para la distribucion nueva
   {
     k[n]=k_min*pow(r,n);
     printf("k[%d]: %lf\n",n,k[n]);
   }
  

 printf("\n");


                           //calculo la nueva distribucion

for(n=0;n<=Nbines-1;n++)
   {    
     for(i=1;i<=filas;i++)
       {
	 if((i > k[n]) && (i <= k[n+1]))
	   {
	     pk[n]=pk[n]+M[i][2];   
	     norma[n]++;

	     //printf(" k[%d]:%lf <  i:%d < k[%d]:%lf\n",n,k[n],i,n,k[n+1]);    
	     //printf("pk[%d]: %lf\n",n,pk[n]);     
	   }
       }       
   }

                      //normalizo 

for(n=0;n<=Nbines-1;n++)
   {   
     if(norma[n] !=0)
       pk[n]=pk[n]/norma[n];
     else       	
       pk[n]=0;
     
     printf("pk[%d]: %lf\n\n",n,pk[n]);
     

   }

 
  printf("\n");


   
  fich2=fopen(file2,"at");
  for (J=0;J<=Nbines;J++)
    {
      printf("k:%lf   pk_bin: %lf \n",k[J],pk[J]);   //ojo con lo que tengo q imprimir!
      

      if(pk[J]!=0)
	{
	  fprintf(fich2,"%lf  %lf\n",k[J],pk[J]);      
	}
    }
  

  
  fclose (fich2);
  
  
  
  
  
}



 




  
 
