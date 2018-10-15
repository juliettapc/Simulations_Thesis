////////////////////////////////////////////////////////////////////
//// Programa para hacer histogramas de las concentraciones  ////////
//// finales (leyendo del archivo) en el programa  especies.c   ///
//// (basado en el programa bines.c                  ///////////////
///////////////////////////////////////////////////////////



#include <stdio.h>
#include <math.h>
#include <stdlib.h>


//        Concentr_vs_c.i._b2.40_bprima1.50_alfa0.0_ro10.15_ro20.70_ro30.15_ro40.00_N4000_1000iter.dat

#define b 2.40
#define bprima 1.50


#define alfa 0.0

#define ro1 0.02
#define ro2 0.9
#define ro3 0.08
#define ro4 0.00

#define N 4000
#define Niter 1000

# define Ninterv  40   // numero de intervalos del histograma (entre c1=0 y c1=1 -->> 1/40=0.025 )   modificar según el rango de c1
# define filas 1000     //que contiene el archivo de entrada

# define S_  1

int J,I,i,j;
float M[filas+1][13];



double p_1[Ninterv+1], p_2[Ninterv+1], p_3[Ninterv+1], p_4[Ninterv+1] ,aux, aux1;
int k_1[Ninterv+1],k_2[Ninterv+1],k_3[Ninterv+1],k_4[Ninterv+1];
double norma_1,norma_2,norma_3,norma_4;

int main()
{ 
  
  FILE *fich2;
  FILE *fich3;

  char file2[256], file3[256];
  




  
  /////////  LECTURA DEL FICHERO  ////////
  
 
  //sprintf(file3,"TCvsK0.00-bini%.3lf.dat",b);
sprintf(file3,"Concentr_vs_c.i._b%.2lf_bprima%.2lf_alfa%.1lf_ro1%.2lf_ro2%.2lf_ro3%.2lf_ro4%.2lf_N%d_%diter.dat",b,bprima,alfa,ro1,ro2,ro3,ro4,N,Niter);
  
 
  //sprintf(file2,"TCvsK0.00-bini%.3lf_binning25.dat",b);
 sprintf(file2,"Histogr_b%.2lf_bprima%.2lf_alfa%.1lf_ro1%.2lf_ro2%.2lf_ro3%.2lf_ro4%.2lf_S_%d_N%d_%diter.dat",b,bprima,alfa,ro1,ro2,ro3,ro4,S_,N,Niter);



  fich2=fopen(file2,"wt");
  fclose(fich2);
  
  
  
  for(J=0;J<=filas;J++)
    {
      for(I=0;I<=12;I++)
	{
	  M[J][I]=0;
	}       
    }


  for(J=0;J<=Ninterv;J++)  
  {
      p_1[J]=0;
      p_2[J]=0;
      p_3[J]=0;
      p_4[J]=0;

      k_1[J]=0;
      k_2[J]=0;     
      k_3[J]=0;
      k_4[J]=0;

     
    }



/*  fich3=fopen(file3,"r");     //leo el fichero y lo guardo en un array

   for(J=1;J<=filas;J++)
    {
      fscanf(fich3,"%d %lf\n",&c1,&c2);
      M[c1][2]=c2;
    }
  
    fclose(fich3);*/
  


//printf("empiezo lectura fichero\n ");
  
  fich3=fopen(file3,"r");     //leo el fichero y lo guardo en un array
  for(J=1;J<=filas;J++)
  {

      fscanf(fich3,"%f   %f    %f    %f   %f   %f    %f    %f   %f   %f   %f   %f\n",&M[J][1],&M[J][2],&M[J][3],&M[J][4],&M[J][5],&M[J][6],&M[J][7],&M[J][8],&M[J][9],&M[J][10],&M[J][11],&M[J][12]);     //esto es asi???????????!!!!!!!!!!!!!!!!    
      
      //    printf("scanf n:%d hecho\n ",J);

  }
  
  fclose(fich3);
  
   
/*  for (J=1;J<=filas;J++)      //imprimo el propio fichero como comprobacion
  { 
      for (I=1;I<=12;I++)
      {
	  
	  printf("%f  ",M[J][I]);
      }
      printf("\n ");
      }*/
  
  
 


 for(i=0;i<=filas;i++)
 {
     for(j=1;j<=Ninterv;j++)       //recuento de eventos en cada intervalo y de cada estrategia
     {
	 aux=Ninterv;
	 aux1=j;


	 if (M[i][5]==0.0)   // para contar por separado los eventos que son cero ESTA BIEN ASÍ (EXCLUYE LOS OTROS IFS??)???????
	 {
	     p_1[0]++;
	     norma_1++;
	     break;
	 }
	 



	 if( (M[i][5]>=(1.0/aux)*(aux1-1.0))  &&   (M[i][5]<(1.0/aux)*aux1) && (M[i][5]!=0)     )
	 {
	     p_1[j]++;
	     norma_1++;
	 }   
	 if (M[i][5]==0.0)   // para contar por separado los eventos que son cero ESTA BIEN ASÍ ???????		 
	 {
	     p_1[0]++;
	     norma_1++;
	 }   
	 
	 
	 if( (M[i][6]>=(1.0/aux)*(aux1-1.0)  &&   M[i][6]<(1.0/aux)*aux1 ) && (M[i][6]!=0)  )
	 {
	     p_2[j]++;
	     norma_2++;
	 } 
	 if (M[i][6]==0.0)   // para contar por separado los eventos que son cero ESTA BIEN ASÍ ???????		 
	 {
	     p_2[0]++;
	     norma_2++;
	 }   
	 
	 
	 if( (M[i][7]>=(1.0/aux)*(aux1-1.0)  &&   M[i][7]<(1.0/aux)*aux1 )  && (M[i][7]!=0) )
	 {
	     p_3[j]++;
	     norma_3++;
	 } 
	 if (M[i][7]==0.0)   // para contar por separado los eventos que son cero ESTA BIEN ASÍ ???????		 
	 {
	     p_3[0]++;
	     norma_3++;
	 }   
	 
	 
	 if( (M[i][8]>=(1.0/aux)*(aux1-1.0)  &&   M[i][8]<(1.0/aux)*aux1 )  && (M[i][8]!=0) )
	 {
	     p_4[j]++;
	     norma_4++;
	 }
	 if (M[i][8]==0.0)   // para contar por separado los eventos que son cero ESTA BIEN ASÍ ???????		 
	 {
	     p_4[0]++;
	     norma_4++;
	 }   
	 
     }
     
 }
 
/* for(i =1; i <= Ninterv; i++)            //recuento histograma
 {
     p_1[k_1[i]]++;
     p_2[k_2[i]]++;
     p_3[k_3[i]]++;
     p_4[k_4[i]]++;
     }*/

 
 for(i =0; i <= Ninterv; i++)            //normalizo
 {
     p_1[i]=p_1[i]/1000.0;
     p_2[i]=p_2[i]/1000.0;
     p_3[i]=p_3[i]/1000.0;
     p_4[i]=p_4[i]/1000.0;
     
 }





 
 fich2=fopen(file2,"at");     //imprimo el histograma
  for (J=0;J<=Ninterv;J++)
    {     
	aux=Ninterv;
	aux1=J;	

	fprintf(fich2,"%lf    %lf    %lf    %lf    %lf    ",(1.0/aux)*aux1,p_1[J],p_2[J],p_3[J],p_4[J]);      
	
	 

/*	if(p_1[J]!=0)      // esto hace falta??????!!!!!!!!	
	    fprintf(fich2,"%lf    ",p_1[J]);      
	else
	    fprintf(fich2,"     ");  

	if(p_2[J]!=0)      // esto hace falta??????!!!!!!!!	
	    fprintf(fich2,"%lf    ",p_2[J]);      
	else
	    fprintf(fich2,"     ");  

	if(p_3[J]!=0)      // esto hace falta??????!!!!!!!!	
	    fprintf(fich2,"%lf    ",p_3[J]);      
	else
	    fprintf(fich2,"     ");  

	if(p_4[J]!=0)      // esto hace falta??????!!!!!!!!	
	    fprintf(fich2,"%lf    ",p_4[J]);      
	else
	fprintf(fich2,"     ");  */
	
	
	fprintf(fich2,"  \n"); 
	
	
    }    
  fclose (fich2);
  
 


}



 




  
 
