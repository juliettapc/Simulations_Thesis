//16-6-09: rutina para hacer las medias sobre rasgos de distintas magnitudes dentro de una iteracion (basado en media_tpos_y_Smax.c)
//17-6-09: renormalizar por el tpo de consenso de cada iter, hacer la media y bienes... buffff

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

# define filas 3000     //que contiene el archivo de entrada (estimacion por encima) DEPENDERÁ DEL TAMAÑO DE LA RED y de Q!!!!!!!!!!!!!!! 

# define columnas 3    
                                                           
# define Niter 20  //por ahora, no hago media sobre iteraciones

#define Nbines 5

# define N 1000


# define Q    170


# define f 10

int J,I,i,j,iter,tpo_cons[Niter+1],indice_bines,lim=0;
float  M[filas+1][columnas+1];
float M_renorm[filas+1][columnas+1];
double Media[Nbines],norma[Nbines],delta_tau=0,aux_indice;




int main()
{ 
    
    
    FILE *fich2;
    FILE *fich3;
    FILE *fich4;
    
    
    char  file2[256], file3[256], file4[256];





    
    ////////// CREACION DE FICHEROS  ////////////

    //  sprintf(file4,"medias/Bines%d_Renorm_Promedio_rasgos_Smax_N%d_Q%d_f%d_ER_%diter.dat",Nbines,N,Q,f,Niter);
    //  sprintf(file4,"medias/Bines%d_Renorm_Promedio_rasgos_Av_p_l_vs_tpo_N%d_Q%d_f%d_ER_%diter.dat",Nbines,N,Q,f,Niter); 
    //  sprintf(file4,"medias/Bines%d_Renorm_Promedio_rasgos_CC_vs_tpo_N%d_Q%d_f%d_ER_%diter.dat",Nbines,N,Q,f,Niter);
     sprintf(file4,"medias/Bines%d_Renorm_Promedio_rasgos_N_Cl_vs_tpo_N%d_Q%d_f%d_ER_%diter.dat",Nbines,N,Q,f,Niter);
        
    fich4=fopen(file4,"wt");
    fclose(fich4);
    




    for(i=0;i<Nbines;i++)
    {
	Media[i]=0.;
	norma[i]=0.;
    }
    
    
    
    for(i=0;i<=Niter;i++)
    {
	tpo_cons[i]=0;
    }
    
    
  	
    delta_tau=1.0/Nbines;  //anchura del intervalo de tiempo renormalizado

    printf("delta_tau:%f",delta_tau);
    // getchar();  
   


    for(iter=1;iter<=Niter;iter++)
    {
	


	/////////  LECTURA DEL FICHERO  ////////      
	
// sprintf(file2,"medias/Promedio_rasgos_Smax_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter);
	//   sprintf(file2,"medias/Promedio_rasgos_Av_p_l_vs_tpo_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter); 
//sprintf(file2,"medias/Promedio_rasgos_CC_vs_tpo_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter);	
	sprintf(file2,"medias/Promedio_rasgos_N_Cl_vs_tpo_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter);
	




    
////////// CREACION DE FICHEROS  ////////////

//	  sprintf(file3,"medias/Renorm_Promedio_rasgos_Smax_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter);
    //  sprintf(file3,"medias/Renorm_Promedio_rasgos_Av_p_l_vs_tpo_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter); 
    //  sprintf(file3,"medias/Renorm_Promedio_rasgos_CC_vs_tpo_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter);
       sprintf(file3,"medias/Renorm_Promedio_rasgos_N_Cl_vs_tpo_N%d_Q%d_f%d_ER_iter%d.dat",N,Q,f,iter);
        
    fich3=fopen(file3,"wt");
    fclose(fich3);
    
//en modo "at" escribo sin borrar lo que ya este escrito
    
    
	
	for(J=0;J<=filas;J++)
	{
	    for(I=0;I<=columnas;I++)
	    {
		M[J][I]=0.;
		M_renorm[J][I]=0.;
	    }       
	}
	
	
	
	
	
	fich2=fopen(file2,"r");     //leo el fichero y lo guardo en un array
	for(J=1;J<=filas;J++)
	{	  
	    for(i=1;i<=columnas;i++)
	    {
		fscanf(fich2,"%f   ",&M[J][i]);		
	    }   
	    	   
	}  
	fclose(fich2);




	for(J=2;J<=filas;J++)   //busco el tpo_cons de la realizacion iter.  (excluyo el primer instante de tpo
	{
	    if( M[J][1]==0 )  
	    {
		tpo_cons[iter]=J-1; //guardo la ultima (no nula, de ahi el -1)linea del archivo 
		lim=tpo_cons[iter];
		//	printf("tpo_cons[%d]=%d\n ",iter,tpo_cons[iter]);
		break;
	    }	 
	}
	
	//getchar();
	


//	printf("fichero leido\n ");				
/*	for (J=1;J<=filas;J++)      //imprimo el propio fichero de entrada como comprobacion
	{ 
	    for (I=1;I<=columnas;I++)
	    {	  		
		printf("%.2f  ",M[J][I]);		
	    }	    
	     getchar();
	    printf("\n ");	    
	}
	getchar();*/
	
	
	
	
	for(I=1;I<=filas;I++)   //copio la inf en otra matriz
	{
	    for(J=1;J<=columnas;J++)
	    {
		M_renorm[I][J]=M[I][J];	
		//printf("%.2f  ",M[J][I]);	
		//getchar();   
		
	    }     
	}
	
	

	for(I=1;I<=filas;I++)   //renormalizo la columna de los tiempo por el tpo cons 
	{
	   
	    M_renorm[I][1]=M_renorm[I][1]/tpo_cons[iter];	
	    // printf("%f  %f  %f  \n",M_renorm[I][1],M_renorm[I][2],M_renorm[I][3]);	
	    // getchar();   
	}
	

	
	
	//imprimo el fichero de las medias con los tiempos renormalizados
	
	fich3=fopen(file3,"wt");                   
	for(j=1;j<=filas;j++)          //calculo la desviacion estandar de la magnitud
	{
	    for(i=1;i<=columnas;i++)
	    {	  
		if(M_renorm[j][2]!=0)
		{  	   
		    fprintf(fich3,"%f  ",M_renorm[j][i]); 	
		}
	
		
	    }
	    if(M_renorm[j][2]!=0)

	    { 
	  	fprintf(fich3,"\n");
	    } 	  
	}
	fclose (fich3);   
	
	

////////////BINES:
    

	for(i=1;i<lim;i++)
	{
	    aux_indice=M_renorm[i][1]/delta_tau;  //la columna 1 contine los tpo renormalizados
	    indice_bines=aux_indice;    //trunco, convirtiendolo en entero

	    Media[indice_bines]+=M_renorm[i][2];//la columna M_renorm[j][2] contine las medias  
	    //printf("M_renorm[%d][2]:%f\n",i,M_renorm[i][2]);     ////OJO!!!!! los indices de  Media[]  van del 0 al Nbines-1

	    //printf("Media[%d]:%f  ",indice_bines,Media[indice_bines]); 
	    //getchar();

	    norma[indice_bines]++;

	}                       



	
	
	
    } //////////////////////fin del bucle a realizaciones


    for(i=0;i<Nbines;i++)
    {
	Media[i]=Media[i]/norma[i];
    }
    
    
    
    //imprimo el fichero de los bines
    
    fich4=fopen(file4,"wt");                   
    for(j=0;j<Nbines;j++)          //calculo la desviacion estandar de la magnitud
    { 	   
	aux_indice=j;
	aux_indice=aux_indice*delta_tau;
	fprintf(fich4,"%f    %f",aux_indice+(delta_tau/2.0),Media[j]); 	
	
	fprintf(fich4,"\n");
	
    }
    fclose (fich4);   
    
    
    
}
