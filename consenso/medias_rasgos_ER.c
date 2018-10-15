//16-6-09: rutina para hacer las medias sobre rasgos de distintas magnitudes dentro de una iteracion (basado en media_tpos_y_Smax.c)

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

# define filas 10000     //que contiene el archivo de entrada (estimacion por encima) DEPENDERÁ DEL TAMAÑO DE LA RED!!!!!!!!!!!!!!! 

# define columnas 14    //que contiene el archivo de entrada (12 el de Ncl,    14 el de S_max, Av_p_l y CCoef)
                                                           

# define promediar 5     // el número de la columna del archivo de entrada a partir de la cual empiezan los datos
// de los distintos rasgos sobre los que promediare  (5 para CC, Av, Smax_f, y 3 para N_Cl)


# define Niter 20  //por ahora, no hago media sobre iteraciones

# define N 1000


# define Q    170



# define f 10

int J,I,i,j,iter,final_archivo=0;
float  M[filas+1][columnas+1];
float Media[filas+1];
double aux,aux1,sigma_2,resta_2;





int main()
{ 
  

  FILE *fich2;
  FILE *fich3;
 

  char  file2[256], file3[256];
  
  
  for(iter=1;iter<=Niter;iter++)
  {
 

////////// CREACION DE FICHEROS  ////////////

         sprintf(file2,"medias/Promedio_rasgos_Smax_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter);
  //   sprintf(file2,"medias/Promedio_rasgos_Av_p_l_vs_tpo_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter); 
//sprintf(file2,"medias/Promedio_rasgos_CC_vs_tpo_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter);
      // sprintf(file2,"medias/Promedio_rasgos_N_Cl_vs_tpo_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter);


  fich2=fopen(file2,"wt");
  fclose(fich2);

//en modo "at" escribo sin borrar lo que ya este escrito

 



      /////////  LECTURA DEL FICHERO  ////////    

  
  sprintf(file3,"Consenso_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter);   
 
//    sprintf(file3,"Av_path_length_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter);   
//sprintf(file3,"CC_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter);   

// sprintf(file3,"N_Clusters_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter);   



  for(i=0;i<=filas;i++)
  {
      Media[i]=0.0;
  }

 
     
  for(J=0;J<=filas;J++)
  {
      for(I=0;I<=columnas;I++)
      {
	  M[J][I]=0;
      }       
  }
  
  
  

  
  fich3=fopen(file3,"r");     //leo el fichero y lo guardo en un array
  for(J=1;J<=filas;J++)
  {	  
      for(i=1;i<=columnas;i++)
      {
	  fscanf(fich3,"%f   ",&M[J][i]);
	  
      }   

      /*if(J>2 || M[J][1]==0 )  // 
      {
	  final_archivo=J; //guardo la ultima linea del archivo
	  }	  */
  }  
  fclose(fich3);


  //printf("fichero leido\n ");



  
   /* for (J=1;J<=filas;J++)      //imprimo el propio fichero de entrada como comprobacion
  { 
      for (I=1;I<=columnas;I++)
      {	  
	
	      printf("%.2f  ",M[J][I]);
	  
      }
      
     
	  printf("\n ");
      
  }
  getchar();*/


  
  
  for(I=1;I<=filas;I++)
  {
      for(J=promediar;J<=columnas;J++)
      {
	  Media[I]=Media[I]+M[I][J];	  	
      }     
  }
  


  for(I=1;I<=filas;I++)
  {
      Media[I]=Media[I]/f;
  //printf("Media:%f   (t%d)",Media[I],I); 
  //getchar();     
  }
  
  
  //printf("media hecha\n ");
  
  
  
      //imprimo el fichero de salida
  
  fich2=fopen(file2,"wt");                   
  for(j=1;j<=filas;j++)          //calculo la desviacion estandar de la magnitud
  {
      
      sigma_2=0;
      resta_2=0;   
      
      for(I=promediar;I<=columnas;I++)
      {
	  resta_2=M[j][I]-Media[j];
	  resta_2=resta_2*resta_2;
	  
	  sigma_2=sigma_2+resta_2;	  	  
      }
      
      
      sigma_2=sigma_2/Niter;
      sigma_2=sqrt(sigma_2);
      
      if(Media[j]!=0)
      {
      fprintf(fich2,"%d    %.2f    %.2f\n",j,Media[j],sigma_2); 

     
      
      // printf("%d    %.2f    %.2f\n",j,Media[j],sigma_2);  
      //getchar();
      }
      
  }
   fclose (fich2);   
  
  
  
  }
  
  
}
