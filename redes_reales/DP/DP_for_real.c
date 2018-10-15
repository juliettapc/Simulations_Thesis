// 24-9-09: Programa para implementar el Dilema del Prisionero sobre Redes Reales

//ojo convenio!!!:  e[i]=1 cooperador,  e[i]=0 defector 


////////// para cada red, cambiar  N, filas y los nombre de los archivos de entrada y de salida




# include <stdio.h>
# include <stdlib.h>
# include <math.h>



# define jugadas 50000  //pasos de tiempo para el juego (transitorio)

# define Niter 500       //estadistica  (numero de ci distintas con para implementar la dinamica)



//# define N  1133     // email.dat

# define N 12722      // SCN.dat
//# define N 198      // jazz.dat
//# define N 13259      // ASTRO.dat   //OJO PQ ESTA EMPIEZA A NUMERAR LOS NODOS EN CERO 
                                      //  y ademas tiene tres columnas!!!





//numero de filas del archivo de entrada == numero total links de la red:



//# define filas 5451   // email.dat
# define filas 39967   // SCN.dat
//# define filas 2742   // jazz.dat
//# define filas  123838   // ASTRO.dat





# define R 1.0
# define Pu 0.0
# define Su  0.0


# define b_min    1.10
# define b_max    1.5
# define delta_b  0.01 

# define rho 0.5    //concentracion inicial de cooperadores

int k[N+1],C[N+1][N+1];
double PK[N+1];        //  para la distribucion P(k)
int  e[N+1],e_aux[N+1];                   //estrategia de los nodos
double  ben[N+1];            //beneficios de los nodos 
int  n_coop;
double  c_media,b;  
int iter;

int D[filas+1][3];  //matriz de pares de links



//////////////Generacion de numeros aleatorios//////////////////
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
///////////////////////////////////////////////////////////////




void inicia_rand(int semilla);
void leer_red();
void leer_red_ASTRO();
void construir_C_k();
void histograma_pk();
void juego();
void establecer_ci();

char nombre1[256],nombre2[256],nombre3[256],nombre4[256],nombre5[256];
  
FILE *archivo1;
FILE *archivo2;
FILE *archivo3;
FILE *archivo4;
FILE *archivo5;

int main()
{
    
    
    printf("\n\nPrograma para la implementacion del DP en redes reales.\n  N=%d  b_min=%.2f b_max=%.2f delta_b=%.2f (%d iter sobre c.i.)\n",N,b_min,b_max,delta_b,Niter);
    
    
    
/////////archivo de entrada


//     sprintf(nombre1,"email_1133nodoskmed9.62.dat");
         sprintf(nombre1,"SCN_12722nodos_kmed6.24.dat");
//    sprintf(nombre1,"jazz_198nodos_kmed26.96.dat");
//    sprintf(nombre1,"ASTRO_13259nodos_kmed18.65.dat");    
    
    
    
    
    
    
    
/////////archivos de salida

// sprintf(nombre3,"C_med_vs_b_ASTRO_b%.2f-%.2f_%diter.txt",b_min,b_max,Niter);   

// sprintf(nombre3,"C_med_vs_b_jazz_b%.2f-%.2f_%diter.txt",b_min,b_max,Niter);   

  sprintf(nombre3,"C_med_vs_b_SCN_b%.2f-%.2f_%diter.txt",b_min,b_max,Niter);     
    
// sprintf(nombre3,"C_med_vs_b_email_b%.2f-%.2f_%diter.txt",b_min,b_max,Niter);     
    archivo3=fopen(nombre3,"wt");
    fclose(archivo3);
    
    
    

// sprintf(nombre4,"Pk_ASTRO_%diter.txt",Niter);  
    
    sprintf(nombre4,"Pk_jazz_%diter.txt",Niter);  
   
// sprintf(nombre4,"Pk_SCN_%diter.txt",Niter);  
    
//    sprintf(nombre4,"Pk_email_%diter.txt",Niter);  
    archivo4=fopen(nombre4,"wt");
    fclose(archivo4);
    




    
    inicia_rand(time(0));
    




    
    
      leer_red();          //solo leo la red real una vez, pq la estadistica sera solo sobre c.i.

      // leer_red_ASTRO();      //es distinta funcion pq este archivo tiene 3 columnas, y los demas solo dos

    construir_C_k();
    histograma_pk();
    



    b=b_min;
    while (b<=b_max)
    {

	if(iter<=2)
	{
//	sprintf(nombre5,"Estrategias_ASTRO_%diter.txt",Niter);  
//sprintf(nombre5,"Estrategias_jazz_%diter.txt",Niter);  
  sprintf(nombre5,"Estrategias_SCN_%diter.txt",Niter);  
	//  sprintf(nombre5,"Estrategias_email_b%.2f_%diter.txt",b,Niter);  
    archivo5=fopen(nombre5,"at");    //para que no borre lo que ya este escrito en el
    fclose(archivo5);
	}





	c_media=0.;
	for(iter=1;iter<=Niter;iter++)
	{

	    printf("iter: %d\n",iter);
	  


      /////////archivos para la C(t) 

	    if(iter<=2)
	    {
		//sprintf(nombre2,"Evol_temp_Coop_ASTRO_b%.2lf_iter%d.txt",b,iter);    
		
//		sprintf(nombre2,"Evol_temp_Coop_jazz_b%.2lf_iter%d.txt",b,iter);    
 sprintf(nombre2,"Evol_temp_Coop_SCN_b%.2lf_iter%d.txt",b,iter);    
		// sprintf(nombre2,"Evol_temp_Coop_email_b%.2lf_iter%d.txt",b,iter);    
		archivo2=fopen(nombre2,"wt");
		fclose(archivo2);
		
	    }



	    establecer_ci();	    
	    juego();

	    c_media=c_media+n_coop;

	}      //fin del bucle de estadistica


	c_media=c_media/(double)(N*Niter);




	archivo3=fopen(nombre3,"at");  
	fprintf(archivo3,"%lf  %lf\n",b,c_media);          
	fclose (archivo3);
	


	b=b+delta_b;
    }     //fin dle bucle de barrer en b



} 




/////////////////////////////
//////////////////////////////////
///////////////////////////////////////


void inicia_rand(int semilla)
{

int i;
int dummy;

srand((unsigned)semilla);
for(i=0;i<111;i++)
  rand();
ip=128;
ip1=ip-24;
ip2=ip-55;
ip3=ip-61;
for(i=0;i<256;i++)
  ira[i]=(unsigned)rand()+(unsigned)rand();
for(i=0;i<1111;i++)
  dummy=RANDOM;
}






/////////////////////////////
//////////////////////////////////
//////////////////////////////////////







void leer_red()
{

    int i,c1,c2;
 
    


    for(i=1;i<=filas;i++)
    {
	D[i][1]=0;
	D[i][2]=0;
    }
   

    archivo1=fopen(nombre1,"r");  
    for(i=1;i<=filas;i++)
    {
	fscanf(archivo1,"%d %d\n",&c1,&c2);

	D[i][1]=c1;
	D[i][2]=c2;
    }
    fclose(archivo1);


    /*  for(i=1;i<=filas;i++)   //comprobacion
    {
	printf("%d-%d\n",D[i][1],D[i][2]);
    }
    getchar();*/
}



//////////////////////////////
///////////////////////////////////////
/////////////////////////////////



void establecer_ci()      //ojo convenio!!!:  e[i]=1 cooperador,  e[i]=0 defector 
{
    int i;
    double r;



    for(i=1;i<=N;i++)    // todos defectores 
    {
	e[i]=0;
    }
    

    n_coop=0;
    for(i=1;i<=N;i++)
    {
	r=FRANDOM;

	if(r<rho)
	{
	    e[i]=1;
	    n_coop++;
	}
    }

    if(iter<=2)     //solo imprimo un par de ejemplos de evolucion temporal C(t)
    {
	archivo2=fopen(nombre2,"at");
	fprintf(archivo2,"0  %d\n",n_coop);      
	fclose(archivo2);
	
    }


}




//////////////////////////////
///////////////////////////////////////
/////////////////////////////////


void juego ()    //ojo!!!!  el convenio es 0=defect y 1=coop!!!!!!!!
{
    int i,j,tpo,y,w,t,k_m;
    double r,p;



    for(i=1;i<=N;i++)
    {
	ben[i]=0;

    }



    tpo=0;
    while(tpo<jugadas)
    {

	for(i=1;i<=N;i++)   //pongo a cero los beneficios de todos
	{
	    ben[i]=0;
	}

	for(i=1;i<=N;i++)   //cada nodo juega con todos sus vecinos
	{
	    for(j=1;j<=k[i];j++)
	    {
		y=C[i][j];
		if(y>i)  //para no repetir parejas de nodos
		{
		    if(e[i]==1 && e[y]==1)   //ambos coop
		    {
			ben[i]+=R;
			ben[y]+=R;
		    }  
		    
		    if(e[i]==1 && e[y]==0)   // i coop, y defect
		    {
			ben[i]+=Su;
			ben[y]+=b;
		    }  
		    
		    if(e[i]==0 && e[y]==1)   //i defect, y coop
		    {
			ben[i]+=b;
			ben[y]+=Su;
		    }  
		    
		    
		    if(e[i]==0 && e[y]==0)   //ambos defect
		    {
			ben[i]+=Pu;
			ben[y]+=Pu;
		    }  


		}

	    }    //i acaba de jugar con sus vecinos
	    
	    	    
	}  // fin de las partidas de todos




	for(i=1;i<=N;i++)   //copio la estrategia de antes de la comparcion
	{
	    e_aux[i]=e[i];  
	}


	for(i=1;i<=N;i++)   //comparan beneficios
	{
	    
	    r=FRANDOM;
	    r=r*k[i];
	    w=(int)r+1;
	    
	    
	    t=C[i][w];    //el vecino elegido al azar
	    
	    if(k[i]>k[t])
	    {k_m=k[i];}
	    else
	    {k_m=k[t];}
	    
	    //printf("juego1 ben[i]=%lf ben[t]=%lf\n",ben[i],ben[t]);//CONTROL
	    
	    
	    if(ben[i]<ben[t])
	    {
		p=(ben[t]-ben[i]);
		
		//CORRECCION en la normalización: diferencia máxima de coef de la matriz de payoff*max k (antes: b*k_m siempre)
		
		if(Su>=0)
		{
		    p=p/(b*(double)k_m);
		}
		else
		{
		    p=p/((b-Su)*(double)k_m);
		}
		
		
		r=FRANDOM;
		if(r<p)
		{
		    e_aux[i]=e[t];
		    
		}
	    }
	    
	    
	 

	}// fin de las comparaciones


	for(i=1;i<=N;i++)   //actualizaciones en paralelo
	{
	    e[i]=e_aux[i];
	}


	n_coop=0;
	for(i=1;i<=N;i++)
	{
	    if(e[i]==1)	    
		n_coop++;	      
	}


	if(iter<=2)     //solo imprimo un par de ejemplos de evolucion temporal C(t)
	{
	    archivo2=fopen(nombre2,"at");
	    fprintf(archivo2,"%d  %d\n",tpo,n_coop);      
	    fclose(archivo2);
	}
	
//	printf("jugada:%d   n_coop:%d\n",tpo,n_coop);	



	if(n_coop==0 || n_coop==N )
	{
	    break;    //me salto el resto de las jugadas.
	}

	tpo++;

    }  //fin bucle while

	printf("jugada:%d   n_coop:%d\n",tpo,n_coop);	

	if(iter<=2)
	{
	    archivo5=fopen(nombre5,"wt"); 
	    for(i=1;i<=N;i++)             //escribo la estrategia de cada nodo al final de cada iter
	    { 
		fprintf(archivo5,"%d  %d\n",i,e[i]);          
	    }
	    fprintf(archivo5,"\n"); //para separar iters  
	    fclose (archivo5);
	}    
	




    
}


////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////



void construir_C_k()      //parto de la matriz de pares de links D[][]
{

    int i,j,d1,d2;

    for(i=1;i<=N;i++)
    {
	k[i]=0;
	for(j=1;j<=N;j++)
	{
	    C[i][j]=0;
	}
    }



    for(i=1;i<=filas;i++)
    {
//	printf("fila%d: ",i);

	d1=D[i][1];    //nodo 1 del link
	d2=D[i][2];    //nodo 2 del link

//	printf("%d-%d\n",d1,d2);


	k[d1]++;
	k[d2]++;
	
	/*printf("k1:%d,  k2:%d\n",k[d1],k[d2]);

	getchar();*/

	C[d1][k[d1]]=d2;
	C[d2][k[d2]]=d1;

    }
    
    /* for(i=1;i<=N;i++)   //comprobacion
    {

	printf("\n%d (k:%d):",i,k[i]);  
	for(j=1;j<=k[i];j++)
	{
	    printf("%d   ",C[i][j]);    
	}
	getchar();
	}*/
}


/////////////////////////////////
/////////////////////////////////////////////
/////////////////////////////////////


void histograma_pk()            //costruccion de P(k)
{
    
    int i;
    
    for(i=0;i<=N;i++)           //inicializo
    {
	PK[i]=0.;
    }
    

    for(i=1; i<=N;i++)            //recuento
    {
	PK[k[i]]++;
    }
    
    



    for(i=1;i<=N;i++)         //normalizo
    {
	PK[i]=PK[i]/N;
	
    }

 
    // Escribo k y PK en un archivo 





    archivo4=fopen(nombre4,"at"); 
    for(i=1; i<N; i++)
    {
	fprintf(archivo4,"%d  %lf\n", i, PK[i]);          
    }
    fclose (archivo4);
    


/*for(i=1; i<N; i++)
  {fprintf(escribe,"%d  %lf\n", i, PK[i]);}*/




}        


/////////////////////////////
//////////////////////////////////
//////////////////////////////////////







void leer_red_ASTRO()    //el fichero ASTRO tiene tres columnas, las dos primeras son las que me interesan
{                      // ademas, empieza a numerar los nodos en 0, no en 1 !!!!!!

    int i,c1,c2,c0;
 
    


    for(i=1;i<=filas;i++)
    {
	D[i][1]=0;
	D[i][1]=0;
	D[i][2]=0;
    }
   

    archivo1=fopen(nombre1,"r");  
    for(i=1;i<=filas;i++)
    {
	fscanf(archivo1,"%d %d %d\n",&c1,&c2,&c0);

	D[i][0]=c0;
	D[i][1]=c1+1;
	D[i][2]=c2+1;
    }
    fclose(archivo1);


    /*  for(i=1;i<=filas;i++)   //comprobacion
    {
	printf("%d-%d\n",D[i][1],D[i][2]);
    }
    getchar();*/
}
