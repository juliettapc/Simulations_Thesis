//16-6-09: rutina para hacer las medias sobre rasgos de distintas magnitudes dentro de una iteracion (basado en media_tpos_y_Smax.c)
//17-6-09: renormalizar por el tpo de consenso de cada iter, hacer la media y bienes
//18-6-09: combinacion de lo que hacian los programas medias_rasgos.c y medias_rasgos_iter.c
// es decir, primero hacer la media sobre rasgos de la magnitud que sea en una iteracion concreta, y luego hacer renormalacion por el tpo_consenso de cada iter y luego hacer la media sobre iter y bienes.

//22-6-09: añado un printf a mano al final del archivo corresp. con el valor para t/tcons=1
//23-6-09: normalizacion de la magnitud por el valor topologico y bines de: media sobre rasgos, media sobre rasgos normalizada al valor topologico, media del valor global sobre iter, y media del valor global sobre iter normalizado al valor topologico

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

# define filas 3000     //que contiene el archivo de entrada (estimacion por encima) DEPENDERÁ DEL TAMAÑO DE LA RED y de Q!!!!!!!!!!!!!!! 

# define columnas  14    //que contiene el archivo de entrada (12 el de Ncl,    14 el de S_max, Av_p_l y CCoef)    


# define promediar 5     // el número de la columna del archivo de entrada a partir de la cual empiezan los datos
// de los rasgos sobre los que promediare  (5 para CC, Av, Smax_f, y 3 para N_Cl)

                                                           
# define Niter 50

#define Nbines 200

# define N 1000


# define Q    10


# define f 10

int J,I,i,j,iter,tpo_cons[Niter+1],indice_bines,lim=0;
float  M[filas+1][columnas+1];
double M_renorm[filas+1][4];
double Media_bines[Nbines],Magnitud_tot_bines[Nbines+1],norma[Nbines],delta_tau=0,aux_indice;

float Media[filas+1];
double Magnitud_tot_bines_topol[Nbines+1],Media_bines_topol[filas+1];
double valor_final, media_valor_final,valor_final_rasgos, media_valor_final_rasgos;
double aux,aux1,sigma_2,resta_2,norma_magnitud;



int main()
{ 
    
    FILE *fich1;
    FILE *fich2;
    FILE *fich3;
    FILE *fich4;    
    FILE *fich6;
    
    char  file1[256],file2[256], file3[256], file4[256], file6[256];
    
    
    
    
    
    ////////// CREACION DEL FICHERO FINAL  ////////////
    
    // sprintf(file4,"medias/Bines%d_Smax_vs_tpo_N%d_Q%d_f%d_SF_%diter.dat",Nbines,N,Q,f,Niter);
    sprintf(file4,"medias/Bines%d_Av_p_l_vs_tpo_N%d_Q%d_f%d_SF_%diter.dat",Nbines,N,Q,f,Niter); 
    // sprintf(file4,"medias/Bines%d_CC_vs_tpo_N%d_Q%d_f%d_SF_%diter.dat",Nbines,N,Q,f,Niter);
    
    //  sprintf(file4,"medias/Bines%d_N_Cl_vs_tpo_N%d_Q%d_f%d_SF_%diter.dat",Nbines,N,Q,f,Niter);
    
    fich4=fopen(file4,"wt");
    fclose(fich4);
          
    
    
    
    
    
    for(i=0;i<Nbines;i++)
    {
	Media_bines[i]=0.;
	Magnitud_tot_bines[i]=0.;
	norma[i]=0.;
	Magnitud_tot_bines_topol[i]=0.0;
	Media_bines_topol[i]=0.0;
    }
    
    
    
    for(i=0;i<=Niter;i++)
    {
	tpo_cons[i]=0;
    }
    
    
    
    
    
    delta_tau=1.0/Nbines;  //anchura del intervalo de tiempo renormalizado
    
    media_valor_final=0.0;
    media_valor_final_rasgos=0.0;
    
    
    for(iter=1;iter<=Niter;iter++)
    {
	
	valor_final=0.0;	
	valor_final_rasgos=0.0;	
	norma_magnitud=0.0;
	
	
	
/////////  LECTURA DEL FICHERO  INICIAL////////    
	
//	sprintf(file1,"Consenso_N%d_Q%d_f%d_alfa0.0_SF_iter%d.dat",N,Q,f,iter);    
	   sprintf(file1,"Av_path_length_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);   
//	sprintf(file1,"CC_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);   
	
//	sprintf(file1,"N_Clusters_N%d_Q%d_f%d_alfa0.0_SF_iter%d.dat",N,Q,f,iter);  
	
	
	
	
	
	
	
////////// CREACION DE FICHEROS INTERMEDIOS ////////////
	
//	sprintf(file2,"medias/Promedio_rasgos_Smax_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
	   sprintf(file2,"medias/Promedio_rasgos_Av_p_l_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter); 
//	sprintf(file2,"medias/Promedio_rasgos_CC_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
	
//	   sprintf(file2,"medias/Promedio_rasgos_N_Cl_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
	
	fich2=fopen(file2,"wt");
	fclose(fich2);
	
//en modo "at" escribo sin borrar lo que ya este escrito
	
	
	
////////// CREACION DE FICHEROS INTERMEDIOS ////////////
	
//	sprintf(file6,"medias/Renorm_Smax_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
	  sprintf(file6,"medias/Renorm_Av_p_l_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter); 
//	  sprintf(file6,"medias/Renorm_CC_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
	
//	sprintf(file6,"medias/Renorm_Promedio_rasgos_N_Cl_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
        
	fich6=fopen(file6,"wt");
	fclose(fich6);
	
//en modo "at" escribo sin borrar lo que ya este escrito
	
	
	
	
	
////////// CREACION DE FICHEROS INTERMEDIOS ////////////
	
//	sprintf(file3,"medias/Renorm_Promedio_rasgos_Smax_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
	  sprintf(file3,"medias/Renorm_Promedio_rasgos_Av_p_l_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter); 
//	  sprintf(file3,"medias/Renorm_Promedio_rasgos_CC_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
	
//	sprintf(file3,"medias/Renorm_Promedio_rasgos_N_Cl_vs_tpo_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
        
	fich3=fopen(file3,"wt");
	fclose(fich3);
	
//en modo "at" escribo sin borrar lo que ya este escrito
	
	
	

	
	
	for(J=0;J<=filas;J++)
	{
	    for(I=0;I<=columnas;I++)
	    {
		M[J][I]=0.;
	    }   
	   
	    M_renorm[J][1]=0.;
	    M_renorm[J][2]=0.;
	    M_renorm[J][3]=0.;
	    
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
	
	
	
	valor_final=M[lim][promediar-1]; //para añadirlo a mano en el archivo de los bines
	media_valor_final+=valor_final;
	
	norma_magnitud=M[lim][2]; //guarda el valor topológico de la magnitud (para normalizar CC y <l> )
	
 
	//printf("norma%f\n ",	norma_magnitud);
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
	
	
	for(I=1;I<=lim;I++)
	{
	    for(J=promediar;J<=columnas;J++)
	    {
		Media[I]=Media[I]+M[I][J];	  	
	    }     
	}
	
	
	
	for(I=1;I<=lim;I++)
	{
	    Media[I]=Media[I]/f;
	    //printf("Media:%f   (t%d)",Media[I],I); 
	    //getchar();     
	}
	
	
	valor_final_rasgos=Media[lim];//para añadirlo a mano en el archivo de los bines
	media_valor_final_rasgos+=valor_final_rasgos;
	
	
	//imprimo el fichero de las medias y la desviacion estandar
	
	fich2=fopen(file2,"wt");                   
	for(j=1;j<=lim;j++)          //calculo la desviacion estandar de la magnitud
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
	
	for(I=1;I<=lim;I++)   //copio la inf de las medias sobre rasgos en otra matriz
	{
	    
	    M_renorm[I][1]=I;
	    M_renorm[I][2]=Media[I];
	    
	    //   M_renorm[I][3]=Media[I]; //guardo la misma magnitud para luego normalizarla por el valor topol
	    
	}
	
	
	
	for(I=1;I<=lim;I++)   //renormalizo solo la columna de los tiempo por el tpo cons 
	{	  
	    M_renorm[I][1]=M_renorm[I][1]/tpo_cons[iter];

	   
	    M_renorm[I][3]=M_renorm[I][2]/norma_magnitud;      //normalizo la magnitud por el valor topológico (solo CC y <l>)
	    
//printf("%f/%f = %f\n",M_renorm[I][2], norma_magnitud,M_renorm[I][3]); 
	    //getchar();

	}
	
	
	
	
	//imprimo el fichero de: tpos_renormalizados al tpo de consenso, medias sobre capas, medias sobre capas_renormalizadas al valor topol de la red 
	
	fich3=fopen(file3,"wt");                   
	for(j=1;j<=lim;j++)          
	{
	    
	    if(M_renorm[j][2]!=0)
	    {  	   
		fprintf(fich3,"%f   %f   %f\n",M_renorm[j][1],M_renorm[j][2],M_renorm[j][3]); 	
	    }			
	}
	  
	
	fclose (fich3);   
	
	
	
	//imprimo el fichero de la magnitud de consenso global  con los tiempos renormalizados
	
	fich6=fopen(file6,"wt");                   
	for(j=1;j<=lim;j++)          
	{
	    
	    fprintf(fich6,"%f   %f   %f\n",M_renorm[j][1],M[j][promediar-1],M[j][promediar-1]/norma_magnitud); 		  				    
	    
	}
	fclose (fich6);   

	
	
	

	
	
	
	
/////////////        BINES      /////////////////////////////////////////////////////////
	
	
	for(i=1;i<lim;i++)
	{
	    aux_indice=M_renorm[i][1]/delta_tau;  //la columna 1 contine los tpo renormalizados
	    indice_bines=aux_indice;    //trunco, convirtiendolo en entero
	    
	    Media_bines[indice_bines]+=M_renorm[i][2];//la columna M_renorm[j][2] contine las medias  (es como Media[])	    
	    Magnitud_tot_bines[indice_bines]+=M[i][promediar-1];
	    
	    
	    Media_bines_topol[indice_bines]+=M_renorm[i][3];      //M_renorm[j][3] = M_renorm[I][2]/norma_magnitud = Media[I]/norma_magnitud
	    Magnitud_tot_bines_topol[indice_bines]+=M[i][promediar-1]/norma_magnitud;

	    //printf("M_renorm[%d][2]:%f\n",i,M_renorm[i][2]);     ////OJO!!!!! los indices de  Media_bines[]  van del 0 al Nbines-1
	    
	    //printf("Media_bines[%d]:%f  ",indice_bines,Media_bines[indice_bines]); 
	    //getchar();
	    
	    norma[indice_bines]++;
	    
	}                       
	
	
	
	
	
	
    } //////////////////////fin del bucle a realizaciones
    
    
    for(i=0;i<Nbines;i++)
    {
	Media_bines[i]=Media_bines[i]/norma[i];
	Magnitud_tot_bines[i]= Magnitud_tot_bines[i]/norma[i];

	Media_bines_topol[i]=Media_bines_topol[i]/norma[i];
	Magnitud_tot_bines_topol[i]= Magnitud_tot_bines_topol[i]/norma[i];
    }
    
    
    
//imprimo el fichero de los bines
    
    fich4=fopen(file4,"wt");                   
    for(j=0;j<Nbines;j++)        
    { 	   
	aux_indice=j;
	aux_indice=aux_indice*delta_tau;
	fprintf(fich4,"%f    %f    %f    %f    %f",aux_indice+(delta_tau/2.0),Media_bines[j],Magnitud_tot_bines[j],Media_bines_topol[j],Magnitud_tot_bines_topol[j]); 	 //la posicion en x es el centro de cada intervalo	
	fprintf(fich4,"\n");
	
    }
    
    fclose (fich4);   
    
   




 

//imprimo el fichero de los bines
    
/*    fich5=fopen(file5,"wt");                   
    for(j=0;j<Nbines;j++)        
    { 	   
	aux_indice=j;
	aux_indice=aux_indice*delta_tau;
	fprintf(fich5,"%f    %f",aux_indice+(delta_tau/2.0),Magnitud_tot_bines[j]); 	 //la posicion en x es el centro de cada intervalo	
	fprintf(fich5,"\n");
	
    }
    
    fprintf(fich5,"1.0    %f",media_valor_final/Niter);//añado a mano el valor final en el archivo de los bines
    fclose (fich5);  */ 
    
    
    
    
    
}
