//4-5-09: calculo de la media sobre iteraciones de los valores finales de Smax para varios valores de Q y los guardo en un solo archivo
//7-5-09: estudio los tiempos de consenso en funcion de Q: su media y su desviación típica


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

# define filas 5000     //que contiene el archivo de entrada (estimacion por encima) DEPENDERÁ DEL TAMAÑO DE LA RED!!!!!!!!!!!!!!! 

# define columnas 14    //que contiene el archivo de entrada (12 el de Av y el de Ncl, 14 el de S_max)
                                                           

# define promediar 4     // el número de la columna del archivo de entrada para el que haré la media (2 para Av y Ncl, 4 para S_max)


# define Niter 100

# define N 1000


# define Q    500



# define f 10

int J,I,i,j,iter;
float  M[filas+1][columnas+1];
float Media,tpo_medio,tpo[Niter+1],S[Niter+1];
double aux,sigma_2_tpo,resta_2_tpo,sigma_2_S,resta_2_S;





int main()
{ 
  
  FILE *fich;
  FILE *fich2;
  FILE *fich3;
  FILE *fich4;

  char file[256], file2[256], file3[256], file4[256];
  
  

 

////////// CREACION DE FICHEROS  ////////////

//sprintf(file2,"medias/Medias_Consenso_N%d_Q%d_f%d_alfa0.0_Dnoblocked.dat",N,Q,f);
//sprintf(file2,"medias/Medias_Av_path_N%d_Q%d_f%d_alfa0.0_Dnoblocked.dat",N,Q,f);
//sprintf(file2,"medias/Medias_N_Clusters_N%d_Q%d_f%d_alfa0.0_Dnoblocked.dat",N,Q,f);

  
 
sprintf(file2,"medias/Medias_Smax_y_Tpos_vs_Q_tpo_N%d_f%d_SF_%diter.dat",N,f,Niter);
  fich2=fopen(file2,"at");
  fclose(fich2);


  sprintf(file,"medias/Individuales_Smax_vs_Q_N%d_f%d_SF_%diter.dat",N,f,Niter);      //en modo "at" escribo sin borrar lo que ya este escrito
  fich=fopen(file,"at");
  fclose(fich);


  sprintf(file4,"medias/Individuales_tpo_vs_Q_N%d_f%d_SF_%diter.dat",N,f,Niter);      //en modo "at" escribo sin borrar lo que ya este escrito
  fich4=fopen(file4,"at");
  fclose(fich4);





  Media=0.0;
  tpo_medio=0.0;
  sigma_2_tpo=0.0;
  resta_2_tpo=0.0;
  sigma_2_S=0.0;
  resta_2_S=0.0;
 

  for(i=0;i<=Niter;i++)
  {
      tpo[iter]=0.0;
      S[iter]=0.0;
  }
  
  for(iter=1;iter<=Niter;iter++)
  {
      
      /////////  LECTURA DEL FICHERO  ////////      
        sprintf(file3,"Consenso_N%d_Q%d_f%d_alfa0.0_SF_iter%d.dat",N,Q,f,iter);   // para las simus de Q>=220
      // sprintf(file3,"Consenso_N%d_Q%d_f%d_Dnoblocked_iter%d.dat",N,Q,f,iter);
                     
                      
      for(J=0;J<=filas;J++)
      {
	  for(I=0;I<=columnas;I++)
	  {
	      M[J][I]=0;
	  }       
      }
      
      // printf("inicializo a cero\n ");
      

      fich3=fopen(file3,"r");     //leo el fichero y lo guardo en un array
      for(J=1;J<=filas;J++)
      {	
	  for(i=1;i<=columnas;i++)
	  {
	      fscanf(fich3,"%f   ",&M[J][i]);
	  }     	  
      }
      
      fclose(fich3);

      //   printf("fichero leido y guardado\n ");

      // printf("fichero leido\n ");
      
      /*  for (J=1;J<=filas;J++)      //imprimo el propio fichero de entrada como comprobacion
      { 
	  for (I=1;I<=columnas;I++)
	  {	  
	      printf("%.2f  ",M[J][I]);
	  }
	  printf("\n ");
      }
      */
      
  
    
      for(i=filas;i>=1;i--)     //busco desde el final hacia el principio, el último valor no nulo de Smax
      {     
	  if(M[i][promediar]!=0)
	  {	      
	      Media=Media+M[i][promediar];
	     
	      tpo_medio=tpo_medio+M[i][1];

	      tpo[iter]=M[i][1];
	      S[iter]=M[i][promediar];

	      fich=fopen(file,"at");             //guardo los Niter valores finales de Smax
	      fprintf(fich,"%d    %.2f\n",Q,M[i][promediar]);        	 
	      fclose (fich);
	    

	      fich4=fopen(file4,"at");             //guardo los Niter valores finales de tpo
	      fprintf(fich4,"%d    %.2f\n",Q,M[i][1]);        	 
	      fclose (fich4);
	      


	      break;
	  }
      }
      
    
      // printf("busqueda completada\n "); 
      
  }  //fin bucle sobre iteraciones
  
  
  
  aux=Niter;                      //hago la media sobre iteraciones
  
  Media= Media/aux;
  tpo_medio=tpo_medio/aux;

                               
  for(j=1;j<=Niter;j++)          //calculo la desviacion estandar de los tiempos de consenso y de la S_max  
  {
      resta_2_tpo=tpo[j]-tpo_medio;
      resta_2_tpo=resta_2_tpo*resta_2_tpo;

      sigma_2_tpo=sigma_2_tpo+resta_2_tpo;


      resta_2_S=S[j]-Media;
      resta_2_S=resta_2_S*resta_2_S;

      sigma_2_S=sigma_2_S+resta_2_S;

  }
	     
 sigma_2_tpo=sigma_2_tpo/aux;
 sigma_2_tpo=sqrt(sigma_2_tpo);


 sigma_2_S=sigma_2_S/aux;
 sigma_2_S=sqrt(sigma_2_S);
 

  
  fich2=fopen(file2,"at");     //imprimo el fichero de salida (el valor de Q y el de la media correspondiente sobre iteraciones de Smax)  
  fprintf(fich2,"%d    %.2f    %.2f    %.2f    %.2f    %d    %d\n",Q,Media,sigma_2_S,tpo_medio,sigma_2_tpo,Niter,f);        	 
  fclose (fich2);
  
  

 


}
