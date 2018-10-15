////////////////////////////////////////////////////////////////////
//// Programa para hacer histogramas Pck(t) (leyendo del archivo)////
//// del programa growth_play_tot_Pck.c                   ///////
//// (basado en el programa histograma.c (utilizado en especies.c)/////// 
///////////////////////////////////////////////////////////



#include <stdio.h>
#include <math.h>
#include <stdlib.h>




//#define delta_b 0.2
#define nb 13 

#define  eps 0.99



#define  tauT 10       //(si tauT=10 y tauD=1, pongo 10 nodos y juego una vez)
#define  tauD 1



#define N 1000
#define Niter 100   

#define tiempos 100   // numero de instantes en los que guardo Pck(t): del orden de 100 pq si no gnuplot no podra con ello
                     //  debe ser: jugadas_extra/x=100

# define K_MAX   500 


int J,I,i,j;
float M[K_MAX+1][tiempos];   //en esta matriz meto los datos del archivo "1000tpo_Nck_10-1-b1.90-e0.80_iter6.dat" (tpo y Num coop de cada k)
float Norma [K_MAX+1][Niter];  //en esta matriz meto los datos del archivo "1000Norma_10-1-b0.00-e0.80_10iter.dat" (tpo y Num nodos con cada k)
float Pck_iter[K_MAX+1][tiempos]; //en esta matriz meto la pkc de cada instante de tiempo, y de todo de una sóla iter
float Pck_sumar[K_MAX+1][tiempos]; //en esta matriz voy acumulando la pkc de cada instante de tiempo, y luego lo divido por Niter
int ib, iter;




                        // OJOOOOOOOOOO!!! INDICES DE 0 A N-1 (no queda una posicion vacia)!!!!!!!
int main()
{ 
    
    FILE *fich2;
    FILE *fich3;
    FILE *fich4;
    FILE *fich5;
    
    char  file2[256],file3[256], file4[256], file5[256];
    
    
    double b=1.1;
    double delta_b=0.1;
    
    for(ib=1;ib<=nb;ib++)
    { 
	printf("b=%lf\n",b);    	  
	
	
	sprintf(file4,"%dNorma_%d-%d-b%.2lf-e%.2lf_%diter.dat",N,tauT,tauD,b,eps,Niter); //archivo que leere (uno por cada b)
	
	sprintf(file5,"%dPck_FINAL_%d-%d-b%.2lf-e%.2lf_%diter.txt",N,tauT,tauD,b,eps,Niter); //archivo FINAL que creare (uno por cada b)		
	fich5=fopen(file5,"wt");
	fclose(fich5);
	
	
	
	
	
	for(i=0;i<=K_MAX;i++)
	{
	    for(j=0;j<tiempos;j++)
	    {
		Pck_sumar[i][j]=0.0;		 
	    }       
	}
	


	for(i=0;i<=K_MAX;i++)
	{
	    for(j=0;j<Niter;j++)
	    {
		Norma[i][j]=0.0;
	    }       
	}


	


	fich4=fopen(file4,"r");     //leo el fichero y lo guardo en un array
	for(j=0;j<Niter;j++)       //OJO CON EL INDICE QUE AVANZA PRIMERO Y EL QUE QUEDA FIJO!!!!!!!!!!!!
	{
	    for(i=0;i<=K_MAX;i++)    //(son k_Max +1 posiciones)
	    {
		if(i<K_MAX)
		{
		    fscanf(fich4,"%f     ",&Norma[i][j]);    
		}
		if(i==K_MAX)
		{
		    fscanf(fich4,"%f\n",&Norma[i][j]);     //esto funcionara asi???????????!!!!!!!!!!!!!!!!  (salto de linea)  
		}
	    }
	    //    printf("scanf n:%d hecho\n ",J);
	    
	}  
	fclose(fich4);
	
	








	
	for(iter=1;iter<=Niter;iter++)
	{
	    
	    printf("  iter:%d\n",iter);    	  
	    
	    sprintf(file3,"%dtpo_Nck_%d-%d-b%.2lf-e%.2lf_iter%d.dat",N,tauT,tauD,b,eps,iter); //archivo que leere (uno por cada b e iter)	  	  


	    sprintf(file2,"%dPck_iter_%d-%d-b%.2lf-e%.2lf_%diter.txt",N,tauT,tauD,b,eps,Niter); //archivo aux que creare (uno por cada b e iter)
	    fich2=fopen(file2,"wt");
	    fclose(fich2);
	    
	    
	    for(i=0;i<=K_MAX;i++)
	    {
		for(j=0;j<tiempos;j++)
		{
		    M[i][j]=0.0;
		}       
	    }
	    
	    
	    
	    for(i=0;i<=K_MAX;i++)
	    {
		for(j=0;j<tiempos;j++)
		{
		    Pck_iter[i][j]=0.0;		 
		}       
	    }
	    
	    
	    
	    fich3=fopen(file3,"r");     //leo el fichero y lo guardo en un array
	    for(j=0;j<tiempos;j++)       //OJO CON EL INDICE QUE AVANZA PRIMERO Y EL QUE QUEDA FIJO!!!!!!!!!!!!
	    {
		for(i=0;i<=K_MAX;i++)    //(son k_Max +1 posiciones)
		{
		    if(i<K_MAX)
		    {
			fscanf(fich3,"%f     ",&M[i][j]);    
		    }
		    if(i==K_MAX)
		    {
			fscanf(fich3,"%f\n",&M[i][j]);     //esto funcionara asi???????????!!!!!!!!!!!!!!!!  (salto de linea)  
		    }
		}
		//    printf("scanf n:%d hecho\n ",J);
		
	    }  
	    fclose(fich3);
	    
	    
	 	    
	    
	    
	    
	    for(j=0;j<tiempos;j++)       //calculo la pkc de la iter    
	    {
		//	printf("(tpo:%d)\n\n",j);    
		for(i=1;i<=K_MAX+1;i++)
		{
		    
		    if(Norma[i][iter-1]!=0.0)
		    {
			Pck_iter[i][j]=M[i][j]/Norma[i][iter-1];    //notar que el segundo indice de Norma[][] esta fijo para una misma iter  (y ojo, pq van de 0 a Niter-1  !!!!!!!!!!
			//printf("%lf/%lf = %lf (iter:%d)\n",M[i][j],Norma[i][iter-1],Pck_iter[i][j],iter-1);    	
			//getchar();

		    }
		    else
		    {
			Pck_iter[i][j]=0.0;
		    }
		}       
	    }
	    
	    
	    
	    for(i=1;i<=K_MAX+1;i++)       //voy acumulando la  pkc (de una b fija)
	    {    //notar que la primera columna es la de los tiempos, y no tengo que sumarla
		for(j=0;j<tiempos;j++)
		{
		    Pck_sumar[i][j]+=Pck_iter[i][j];   
		}       
	    }
	    
	   
	    fich2=fopen(file2,"at");     //imprimo el Pck iter (solo como comprobacion) 
	    for(j=0;j<tiempos;j++)     
	    {    
		for(i=0;i<=K_MAX;i++)     //ANTES DE 1  A k_max
		{
		    if(i==0)
		    {
			fprintf(fich2,"%d     ",j);    	  
		    }
		    
		    if(i<K_MAX)
		    {
			fprintf(fich2,"%lf     ",Pck_iter[i][j]);    	  
		    }
		    
		    if(i==K_MAX)
		    {
			fprintf(fich2,"%lf\n",Pck_iter[i][j]);    //(salto de linea)  
		    }
		    
		}		  
	    }
	    fclose(fich2);
	    
	    
	    
	    
	    
	    
	}   ///FIN BUCLE EN ITER
	
	
	
	
	for(j=0;j<tiempos;j++)
	{
	    Pck_sumar[0][j]=j;   //relleno la columna de los tiempos 
	}   
	
	
	
	fich5=fopen(file5,"at");     //normalizo e imprimo el Pck de una b fija  (en bloques de un tiempo fijo serparados por una linea en blanco)
	for(j=0;j<tiempos;j++)     
	{    
	    for(i=1;i<=K_MAX;i++)     //ANTES DE 1  A k_max
	    {	 
		fprintf(fich5,"%d     %d     %lf\n",j,i,Pck_sumar[i][j]/Niter);    	  
	    }
	    fprintf(fich5,"\n");    	  
	}
	fclose(fich5);
	
     
	
	
	b=b+delta_b;           	  
	
	
    }
}








  
 
