#define DEBUG

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

 
# define  PI 3.141592653589793
# define  N 32
# define  alfa 0.250


double x[N+1], v[N+1], a[N+1];
double  delta_v, t, d_tiempo, E[N+1], Energia, Q[N+1], Q_punto[N+1], posicion , aux, omega_2[N+1];
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
    j_max_externo=4000;
   


    sprintf(nombre,"cadena_fpu%.3lf.dat",alfa);
    escribe=fopen(nombre,"wt");
    escribe2=fopen("frecuencias_nat_fpu.dat","wt");


    for(i=1;i<=N;i++)        // condic. iniciales: excitado con un seno de lanbda=2N
    {                     //   y  extremos fijos
	    x[i]=0.0;
	    posicion=i;
	    v[i]=sin(PI*(posicion-1.0)/(N-1.0));
	    //printf("v(%d)=%f \n",i,v[i]);
	    //printf("x(%d)=%f \n",i,x[i]);
    }

 
    printf("v[1]=%f,v[%d]=%f\n",v[1],N,v[N]);

    for(q=1;q<=N;q++)
    {
	omega_2[q]=2.0*sin(PI*q/(2.0*(N+1.0)));
	omega_2[q]=omega_2[q]*omega_2[q];

	fprintf(escribe2,"\n w(%d)= %f ",q,omega_2[q]);
    }

    fclose(escribe2);
    
    d_tiempo=2.0*PI/(k_max_interno*sqrt(omega_2[N]));

    //printf("d_tiempo=%f",d_tiempo);
    
    
    t=0.0;
    j=0;
    
    while(j<=j_max_externo)          //bucle temporal externo (escribire cada vez)
    {
	printf("\n t= %f ",t);

	
	Energia=0.0; 
	for(q=1;q<=N;q++)         //bucle a los N modos normales
	{
	    //printf("\n  x(%d)= %lf ",q,x[q] );
	    //printf("\n  v(%d)= %lf ",q,v[q] );
	    
	    Q[q]=funcion_Q(q);       //llamo a las funciones que calculan cada cosa
	    Q_punto[q]=funcion_Q_punto(q);

	    	    
	    E[q]=(Q_punto[q]*Q_punto[q]+omega_2[q]*Q[q]*Q[q])/2.0;   //calculo la energia del modo q

	    

	    //E[q]=0.5*v[q]*v[q]+0.5*x[q]*x[q];
	    Energia= Energia+E[q];  

	}
	

/*	for(q=1;q<=N;q++)           //suma a todos los modos normales para obtener la energia total
	{
	    printf("\n E(%d)= %f, Q=%f,Q_punto=%f x=%f v=%f ",q,E[q],Q[q],Q_punto[q],x[q],v[q]);
	}
*/	
//	printf("\n Energia(%lf)= %f ",t,Energia);
	fprintf(escribe,"\n%f   %f",t, Energia); 

	for(q=1;q<=N;q++)          
	{
	    fprintf(escribe,"%f   ",E[q]/Energia); 
	    fflush(escribe);
	}
	fprintf(escribe,"\n"); 	



/*#ifdef DEBUG
	fevol=fopen("evol.plt","wt");
	#endif*/

	for(k=1;k<=k_max_interno;k++)         //bucle temporal interno (solo calculo, no escribo)
	{
	    
	    //algoritmo euler
	    for(i=1;i<=N;i++) 
	    {
		a[i]=aceleracion(i);
		delta_v=a[i]*d_tiempo;
		v[i]+=delta_v;
	    }
	    for(i=1;i<=N;i++)       //bucle a las N particulas          
	    {
		x[i]+=d_tiempo*v[i];
		
		
/*#ifdef DEBUG
		if(k==1)
		fprintf(fevol,"\n  %d %lf %lf %lf",i,x[i],v[i],a[i] );
#endif
		//printf("\n  v(%d)= %lf ",i,v[i] );
		//printf("\n  a(%d)= %lf ",i,a[i] );*/

	    }
	    t+=d_tiempo;     //este tiempo solo se usara  para escribirlo en el fichero	
    
	   
	}
	
	   
/*#ifdef DEBUG
	fclose(fevol);
	getchar();
	#endif*/

	printf("\n j= %d\n ",j);
	j++;
 
    }   
    
    
    fclose(escribe);
}




double aceleracion (n)       //al llamar a la funcion, le paso el numero de la particula en la que estoy 
{

    double a;

    //printf("\n  x(%d)= %lf ",i,x[i] );

    if( n==1 || n==N)         //extremos fijos
	{   
	    a=0.0;
	}
    else                 //calculo la aceleracion de la particula correspondiente
    {
	a=(x[n+1]-2.0*x[n]+x[n-1]);//+alfa*((x[n+1]-x[n])*(x[n+1]-x[n]) -(x[n]-x[n-1])*(x[n]-x[n-1]));
    }    

//printf("\n  a(%d)= %lf ",n,a);

    return a;    //devuelvo al programa principal el valor obtenido
}




double funcion_Q (modo)        //al llamar a la funcion, le paso el numero del modo
{
    double resultado;
    int m;

    resultado=0.0;

//printf("\n  x(%d)= %lf ",q,x[q] );

    for(m=1;m<=N;m++)
    {
	resultado+=(x[m]*sin(PI*modo*m/(N+1.0)));  //calculo el sumatorio extendido a todas las particulas
    }

    resultado=resultado*sqrt(2.0/(N+1.0));    //lo multiplico por la raiz

    return resultado;         //devuelvo al programa principal el valor obtenido
}




double funcion_Q_punto (modo)          //al llamar a la funcion, le paso el numero del modo
{
    double resultado;
    int m;

    resultado=0.0;

//printf("\n  v(%d)= %lf ",q,v[q] );
    for(m=1;m<=N;m++)
    {
	

	resultado+=v[m]*sin(PI*modo*m/(N+1.0));   //calculo el sumatorio extendido a todas las particulas
	
    }
    resultado=resultado*sqrt(2.0/(N+1.0));      //lo multiplico por la raiz

    return resultado;               //devuelvo al programa principal el valor obtenido
}




