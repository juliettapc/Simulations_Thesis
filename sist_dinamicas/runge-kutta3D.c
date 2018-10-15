			//programa para implementar el metodo de Runge-Kutta de resolucion
			//de ecuaciones diferenciales (en tres dimensiones)



# include <stdio.h>
# include <stdlib.h>
# include <math.h>
//#define valor_real 0,367879441171442    // valor de 1/e

# define N 3


double x_old[N],  k_1[N], k_2[N], k_3[N], k_4[N];
double  t_final, t, t_o,argumento, delta_t;
double funcion(double );
int i,j;



FILE *escribe;


int main()
{

 escribe=fopen("runge.dat", "wt");

for(i=0;i<N;i++)
    {
    x_old[i]=1;
    }        //c.i.
t=0;

delta_t=0.001;
t_final=1;
printf("condiciones iniciales\n");


    for(i=0;i<N;i++)
    {
	printf("i=%d\n",i);
	k_1[i]=0;
	k_2[i]=0;
	k_3[i]=0;
	k_4[i]=0;
    }
  
	
    printf("t_final=%lf  delta_t=%lf\n",t_final,delta_t);
	    
while(t<=t_final)
  {

//printf("k_1=%lf k_2=%lf k_3=%lf k_4=%lf, t=%lf\n\n",k_1, k_2, k_3, k_4, t);

  for(i=0;i<N;i++)
     {
	
     k_1[i]=funcion(x_old[i])*delta_t;


     k_2[i]=funcion(x_old[i]+0.5*k_1[i])*delta_t;

	
     k_3[i]=funcion(x_old[i]+0.5*k_2[i])*delta_t;


     k_4[i]=funcion(x_old[i]+k_3[i])*delta_t;
     
	 //printf("x(t=%lf)=%lf\n",t,x_old[i]);

 
     x_old[i]=x_old[i] + (1.0/6.0)*(k_1[i]+2*k_2[i]+2*k_3[i]+k_4[i]);

     printf("x(t=%lf)=%lf\n",t,x_old[i]);
     fprintf(escribe,"%lf  ",x_old[i]); 
     }
      //ojo con la fraccion, pq 1/6=0!!!!!

 fprintf(escribe,"\n"); 

   t+=delta_t;

 for(i=0;i<N;i++)
    {
	printf("k_1=%lf k_2=%lf k_3=%lf k_4=%lf\n",k_1[i], k_2[i], k_3[i], k_4[i]);
    }

 }            //fin del algoritmo


printf("\nel resultado es\nx(%lf)=(",t);
for(i=0;i<N;i++)
     {
	 printf("%lf,",x_old[i]);

     }
printf(")\n"); 
//getch();


 fclose(escribe);
}


double funcion(double x)             //corresponde a la ec.dif.: f(x)=dx/dt=-x
{
    double resultado;
 
 if(i==0)     //coord. x
 {                    
 resultado=-x;
 return (resultado); 
 }
 else 
  {
     if(i==1)      //coord. y
     {
     resultado=-x;
     return (resultado); 
     }
     else
      {
	 if(i==2)      //coord. z
	 {
	  resultado=-x;
	  return (resultado); 
	 }
      }
  }
 
//printf("x=%lf f(x)=%lf\n",x,valor);

}
