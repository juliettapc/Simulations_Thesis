////////////////////////////////////////////////////////////////////////////////////
///  Programita para implementar varias ecuaciones de la dinámica evolutiva:   ////
// cuasiespecies, replicador y replicador-mutador                             ////
// utilizando un algoritmo de Runge-Kutta (varias dimensiones= nº especies)  ////
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
# include <stdlib.h>



# define T 1.35
#define R 1.0
#define Pu 0.0
#define Su -0.3

#define especies 2



//Generacion de numeros aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];



char nombre [256];
FILE *escribe;


void inicia_rand(int semilla);
double funcion(double);


double fitness[especies], mutaciones[especies][especies],fit_media, x_old[especies];
 double     k_1[especies], k_2[especies], k_3[especies], k_4[especies];
int i;

 main()
{
    
    sprintf(nombre,"evolucion.dat");
    escribe=fopen(nombre,"wt");

    
   
    int j; 
    double delta_t, t_final, t;
    double norma,suma;  //para normalizar la x_old inicial (frecuencias)
    double norma_m[especies];  //para normalizar la matriz de mutaciones

    
    inicia_rand(time(0));
    
    
    t=0;
    delta_t=0.001;
    t_final=10;
    
    for(i=0;i<especies;i++)     //inicializaciones a cero
    {
	printf("i=%d\n",i);
	k_1[i]=0;
	k_2[i]=0;
	k_3[i]=0;
	k_4[i]=0;

	fitness[i]=0.0;
	norma_m[i]=0.0;
	x_old[i]=0.0;

	for(j=0;j<especies;j++)
	{
	    mutaciones[i][j]=0.0; 
	   
	}
    }
    fit_media=0.0;


    norma=0.0;

    for(i=0;i<especies;i++)    // frecuencias iniciales aleat.
    {
    x_old[i]=FRANDOM;
    norma+=x_old[i];
    }
   
    suma=0.0;
    for(i=0;i<especies;i++)    // normalizo
    {
	x_old[i]=x_old[i]/norma;
	suma+=x_old[i];
    }
    printf("x1=%lf  x2=%lf  suma = %lf",x_old[0],x_old[1],suma);
    getchar();


    
    for(i=0;i<especies;i++)
    {
	fitness[i]=0.5;    //por ahora, landscape cte
    }
    
    
    
    for(i=0;i<especies;i++)       //matriz de mutaciones (tiene que ser estocastica????????????????????)
    { 

	norma_m[i]=0.0;
	    
	for(j=0;j<especies;j++)
	{
	   mutaciones[i][j]=FRANDOM;
	   norma_m[i]+= mutaciones[i][j];
	}

	for(j=0;j<especies;j++)   //normalizo (matriz estocástica)
	{
	    mutaciones[i][j]=mutaciones[i][j]/norma_m[i];
	}


    } 
    

    for(i=0;i<especies;i++)       //calculo la fitness media de la poblacion
    {	
	fit_media+=fitness[i]+x_old[i];
	    
    }

    while(t<=t_final)    //inicio del algoritmo R-K
    { 
	
//printf("k_1=%lf k_2=%lf k_3=%lf k_4=%lf, t=%lf\n\n",k_1, k_2, k_3, k_4, t);
	
  for(i=0;i<especies;i++)
     {
	
     k_1[i]=funcion(x_old[i])*delta_t;


     k_2[i]=funcion(x_old[i]+0.5*k_1[i])*delta_t;

	
     k_3[i]=funcion(x_old[i]+0.5*k_2[i])*delta_t;


     k_4[i]=funcion(x_old[i]+k_3[i])*delta_t;
     
	 //printf("x(t=%lf)=%lf\n",t,x_old[i]);

 
     x_old[i]=x_old[i] + (1.0/6.0)*(k_1[i]+2*k_2[i]+2*k_3[i]+k_4[i]);

     printf("x(t=%lf)=%lf\n",t,x_old[i]);
     fprintf(escribe,"%lf  ",x_old[i]); 


     if (i==(especies-1))
	 fprintf(escribe,"%lf  ",x_old[0]+x_old[1]); 


     printf("x1=%lf  x2=%lf  suma = %lf \n",x_old[0],x_old[1],x_old[0]+x_old[1]);  //comprobacion de la condicion de cierre
     }
      //ojo con la fraccion, pq 1/6=0   !!!!!

 fprintf(escribe,"\n"); 

   t+=delta_t;

/* for(i=0;i<especies;i++)
    {
	printf("k_1=%lf k_2=%lf k_3=%lf k_4=%lf\n",k_1[i], k_2[i], k_3[i], k_4[i]);
	}	*/


    }  //fin del algoritmo R-K
    
    fclose(escribe);
}



////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

double funcion(double x)       //corresponde a la ec. cuasiespecies
{

 double resultado,primero;
 int j;

 
 if(i==0)     // especie 0
 {
     
     primero=0.0;
     
     for(j=0;j<especies;j++)       
     {
	 primero+=x_old[j]*fitness[j]*mutaciones[j][i];
     }
     
     resultado=primero - x_old[i]*fit_media;
     
     
     return (resultado); 
 }
 else 
  {
     if(i==1)      // especie 1 

     {
     primero=0.0;
     
     for(j=0;j<especies;j++)       
     {
	 primero+=x_old[j]*fitness[j]*mutaciones[j][i];
     }
     
     resultado=primero - x_old[i]*fit_media;
     
     return (resultado); 
     }
     else
      {
	 if(i==2)      // especie 2
	 {
	  primero=0.0;
     
     for(j=0;j<especies;j++)      
     {
	 primero+=x_old[j]*fitness[j]*mutaciones[j][i];
     }
     
     resultado=primero - x_old[i]*fit_media;
     
	  return (resultado); 
	 }
      }
  }
 
//printf("x=%lf f(x)=%lf\n",x,valor);

}


////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


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
