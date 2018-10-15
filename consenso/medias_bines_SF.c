//16-6-09: rutina para hacer las medias sobre rasgos de distintas magnitudes dentro de una iteracion (basado en media_tpos_y_Smax.c)
//17-6-09: renormalizar por el tpo de consenso de cada iter, hacer la media y bienes
//18-6-09: combinacion de lo que hacian los programas medias_rasgos.c y medias_rasgos_iter.c
// es decir, primero hacer la media sobre rasgos de la magnitud que sea en una iteracion concreta, y luego hacer renormalacion por el tpo_consenso de cada iter y luego hacer la media sobre iter y bienes.



#include <stdio.h>
#include <math.h>
#include <stdlib.h>

# define filas 3000     //que contiene el archivo de entrada (estimacion por encima) DEPENDERÁ DEL TAMAÑO DE LA RED y de Q!!!!!!!!!!!!!!! 

# define columnas  14    //que contiene el archivo de entrada (12 el de Ncl,    14 el de S_max, Av_p_l y CCoef)    


# define promediar 5     // el número de la columna del archivo de entrada a partir de la cual empiezan los datos
// de los rasgos sobre los que promediare  (5 para CC, Av, Smax_f, y 3 para N_Cl)

                                                           
# define Niter 3 

#define Nbines 200

# define N 1000


# define Q    10


# define f 10

int J,I,i,j,iter,tpo_cons[Niter+1],indice_bines,lim=0;
float  M[filas+1][columnas+1];
float M_renorm[filas+1][columnas+1];
double Media_bines[Nbines],norma[Nbines],delta_tau=0,aux_indice;
float Media[filas+1];
double aux,aux1,sigma_2,resta_2;



int main()
{ 
    
    FILE *fich1;
    FILE *fich2;
    FILE *fich3;
    FILE *fich4;
    
    
    char  file1[256],file2[256], file3[256], file4[256];
    
    
    
    

    
    ////////// CREACION DEL FICHERO FINAL  ////////////
    
    // sprintf(file4,"medias/Bines%d_Renorm_Promedio_rasgos_Smax_vs_tpo_N%d_Q%d_f%d_SF_%diter.dat",Nbines,N,Q,f,Niter);
    //   sprintf(file4,"medias/Bines%d_Renorm_Promedio_rasgos_Av_p_l_vs_tpo_N%d_Q%d_f%d_SF_%diter.dat",Nbines,N,Q,f,Niter); 
     sprintf(file4,"medias/Bines%d_Renorm_Promedio_rasgos_CC_vs_tpo_N%d_Q%d_f%d_SF_%diter.dat",Nbines,N,Q,f,Niter);

    //  sprintf(file4,"medias/Bines%d_Renorm_Promedio_rasgos_N_Cl_vs_tpo_N%d_Q%d_f%d_SF_%diter.dat",Nbines,N,Q,f,Niter);
    
    fich4=fopen(file4,"wt");
    fclose(fich4);
    



    
    for(i=0;i<Nbines;i++)
    {
	Media_bines[i]=0.;
	norma[i]=0.;
    }
    
    
    
    for(i=0;i<=Niter;i++)
    {
	tpo_cons[i]=0;
    }
    
    
    
    
    
    delta_tau=1.0/Nbines;  //anchura del intervalo de tiempo renormalizado
    
    
    
    
    
    for(iter=1;iter<=Niter;iter++)
    {
	
	
/////////  LECTURA DEL FICHERO  INICIAL////////    
                
//sprintf(file1,"Consenso_N%d_Q%d_f%d_alfa0.0_SF_iter%d.dat",N,Q,f,iter);    
//	   sprintf(file1,"Av_path_length_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);   
sprintf(file1,"CC_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);   
        	     
//	sprintf(file1,"N_Clusters_N%d_Q%d_f%d_alfa0.0_SF_iter%d.dat",N,Q,f,iter);  
	
	
	
	
	
	
	
////////// CREACION DE FICHEROS INTERMEDIOS ////////////
	
//	   sprintf(file2,"medias/Promedio_rasgos_Smax_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
//	   sprintf(file2,"medias/Promedio_rasgos_Av_p_l_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter); 
 sprintf(file2,"medias/Promedio_rasgos_CC_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);

//	   sprintf(file2,"medias/Promedio_rasgos_N_Cl_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
	
 fich2=fopen(file2,"wt");
 fclose(fich2);
 
//en modo "at" escribo sin borrar lo que ya este escrito
	
	
	
	
	
	
	
////////// CREACION DE FICHEROS INTERMEDIOS ////////////
	
//	  sprintf(file3,"medias/Renorm_Promedio_rasgos_Smax_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
//	  sprintf(file3,"medias/Renorm_Promedio_rasgos_Av_p_l_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter); 
	  sprintf(file3,"medias/Renorm_Promedio_rasgos_CC_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);

//	sprintf(file3,"medias/Renorm_Promedio_rasgos_N_Cl_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
        
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
	    Media[J]=0.0;
	}
	
	
	
	
	
	fich1=fopen(file1,"r");     ////////////////////////leo el fichero inicial y lo guardo en una matriz
	for(J=1;J<=filas;J++)
	{	  
	    for(i=1;i<=columnas;i++)
	    {
		fscanf(fich1,"%f   ",&M[J][i]);		
	    }   
	    
	}  
	fclose(fich1);
	
	
	
	
	for(J=2;J<=filas;J++)   /////////busco el tpo_cons de la realizacion iter.  (excluyo el primer instante de tpo)
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
	
	
	
/////////////////////////REALIZACION DE LA MEDIA SOBRE RASGOS///////////////////////	
	
	
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
	
	
	
	
	
	//imprimo el fichero de las medias y la desviacion estandar
	
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
	
	
	
	
	
	
	
	
	
	
	
////////////////////////RENORMALIZACION DE LOS TIEMPOS//////////////////////////
	
	for(I=1;I<=filas;I++)   //copio la inf de las medias sobre rasgos en otra matriz
	{
	    for(J=1;J<=columnas;J++)
	    {
		M_renorm[I][1]=I;
		M_renorm[I][2]=Media[I];	
	
	    }     
	}
	
	
	
	for(I=1;I<=filas;I++)   //renormalizo solo la columna de los tiempo por el tpo cons 
	{	  
	    M_renorm[I][1]=M_renorm[I][1]/tpo_cons[iter];		    
	}
	
	
	
	
	//imprimo el fichero de las medias con los tiempos renormalizados
	
	fich3=fopen(file3,"wt");                   
	for(j=1;j<=filas;j++)          
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
	
	

/////////////        BINES      /////////////////////////////////////////////////////////
    

	for(i=1;i<lim;i++)
	{
	    aux_indice=M_renorm[i][1]/delta_tau;  //la columna 1 contine los tpo renormalizados
	    indice_bines=aux_indice;    //trunco, convirtiendolo en entero

	    Media_bines[indice_bines]+=M_renorm[i][2];//la columna M_renorm[j][2] contine las medias  

	    //printf("M_renorm[%d][2]:%f\n",i,M_renorm[i][2]);     ////OJO!!!!! los indices de  Media_bines[]  van del 0 al Nbines-1

	    //printf("Media_bines[%d]:%f  ",indice_bines,Media_bines[indice_bines]); 
	    //getchar();

	    norma[indice_bines]++;

	}                       



	
	
	
    } //////////////////////fin del bucle a realizaciones


    for(i=0;i<Nbines;i++)
    {
	Media_bines[i]=Media_bines[i]/norma[i];
    }
    
    
    
    //imprimo el fichero de los bines
    
    fich4=fopen(file4,"wt");                   
    for(j=0;j<Nbines;j++)        
    { 	   
	aux_indice=j;
	aux_indice=aux_indice*delta_tau;
	fprintf(fich4,"%f    %f",aux_indice+(delta_tau/2.0),Media_bines[j]); 	 //la posicion en x es el centro de cada intervalo
	
	fprintf(fich4,"\n");
	
    }
    fclose (fich4);   
    
    
    
}
