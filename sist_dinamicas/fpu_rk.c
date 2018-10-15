#define DEBUG

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

 
# define  PI 3.141592653589793
# define  N 32
# define  alfa 0.250


double x[N+1],x_old[N+1], v[N+1],v_old[N+1], a[N+1];
double  delta_v, t, d_tiempo, E[N+1], Energia, E_media, Q[N+1], Q_punto[N+1], posicion , aux, omega_2[N+1];
int i,j,k ,k_max_interno, j_max_externo, q;
char nombre [256];


double aceleracion(int);
double funcion_Q(int), funcion_Q_punto(int);


FILE *escribe;
FILE *escribe2;


int main ()
{

    FILE *fevol;
   

    k_max_interno=200;
    j_max_externo=10000;
   


    sprintf(nombre,"cadena_fpu%.3lf.dat",alfa);
    escribe=fopen(nombre,"wt");
    escribe2=fopen("frecuencias_nat_fpu.dat","wt");


    for(i=1;i<=N;i++)        // condic. iniciales: excitado con un seno de lanbda=2N
    {                     //   y  extremos fijos
	x[i]=0.0;
	v[i]=sin(PI*(i-1.0)/(N-1.0));
    }
    
    
    
    for(q=1;q<=N;q++)
    {
	omega_2[q]=2.0*sin(PI*q/(2.0*(N+1.0)));
	omega_2[q]=omega_2[q]*omega_2[q];
	
	fprintf(escribe2,"\n w(%d)= %f ",q,omega_2[q]);
    }    
    fclose(escribe2);
    


    d_tiempo=2.0*PI/(k_max_interno*sqrt(omega_2[N]));
  
    t=0.0;
    j=0;
    
    while(j<=j_max_externo)          //bucle temporal externo (escribire cada vez)
    {
	printf("\n t= %f ",t);

	
	Energia=E_media=0.0; 
	for(q=1;q<=N;q++)         //bucle a los N modos normales
	{
	    Q[q]=funcion_Q(q);       //llamo a las funciones que calculan cada cosa
	    Q_punto[q]=funcion_Q_punto(q);
	    	    
	    E[q]=0.5*(Q_punto[q]*Q_punto[q]+omega_2[q]*Q[q]*Q[q]);   //calculo la energia del modo q
	    Energia+=E[q];  
	}

	fprintf(escribe,"\n%f   %f",t, Energia); 
 	for(q=1;q<=N;q++)          
	{
	    fprintf(escribe,"%f   ",E[q]); 
	    fflush(escribe);
	}
	fprintf(escribe,"\n"); 	

	printf("\n Energia(%f)= %f ",t,Energia);
	
       
	for(k=1;k<=4;k++)
	{
	    
	    //algoritmo rk
	    for(i=1;i<=N;i++)         
	    {
		a[i]=aceleracion(i);	
		v_old[i]=v[i];
		v[i]+=0.5*d_tiempo*a[i];
	    }
	    
	    for(i=1;i<=N;i++)         
	    {			
		x[i]+=0.5*d_tiempo*v_old[i];
	    }
	    
	    t+=d_tiempo*0.5;     //este tiempo solo se usara  para escribirlo en el fichero	
	    
	}
	
	
	printf("\n j= %d\n ",j);
	j++;
	
    }   
    
    
    fclose(escribe);
}




double aceleracion (n)       //al llamar a la funcion, le paso el numero de la particula en la que estoy 
{
    double a;

 
    if( n==1 || n==N)         //extremos fijos
	{   
	    a=0.0;
	}
    else                 //calculo la aceleracion de la particula correspondiente
    {
	a=(x[n+1]-2.0*x[n]+x[n-1])+alfa*((x[n+1]-x[n])*(x[n+1]-x[n])-(x[n]-x[n-1])*(x[n]-x[n-1]));
    }    


    return a;    //devuelvo al programa principal el valor obtenido
}




double funcion_Q (modo)        //al llamar a la funcion, le paso el numero del modo
{
    double resultado;
    int m;

    resultado=0.0;

    for(m=1;m<=N;m++)
    {
	resultado+=(x[m]*sin(PI*modo*m/(N+1.0)));  //calculo el sumatorio extendido a todas las particulas
    }

    resultado=resultado*sqrt(2.0/(N+1.0));    //lo multiplico por la raiz

    return resultado;         //devuelvo al programa principal el valor obtenido
}




double funcion_Q_punto (modo)          //al llamar a la funcion, le paso el numero del modo
{
    double resultado1;
    int m;

    resultado1=0.0;

    for(m=1;m<=N;m++)
    {
	resultado1+=v[m]*sin(PI*modo*m/(N+1.0));   //calculo el sumatorio extendido a todas las particulas
	
    }
    resultado1=resultado1*sqrt(2.0/(N+1.0));      //lo multiplico por la raiz

    return resultado1;               //devuelvo al programa principal el valor obtenido
}




