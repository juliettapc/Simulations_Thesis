			//programa para resolver las ecuaciones de Lorenz y obtencion del
			//atactor extraño en un fichero .dat (tres columnas de datos: (x,y,z) )


# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# define N 3


double x_old[N],  k_1[N], k_2[N], k_3[N], k_4[N];
double  t_final, t, t_o,argumento, delta_t;
double funcion(double );
int i,j;
char nombre[256];



FILE *escribe;

double sigma, r, b;
  
 

int main()
{
  sigma=10;          //con estos valores para los parametros, es caótico
  r=28;
  b=8.0/3.0;
  
  delta_t=0.005;
  t_final=50;                


 sprintf(nombre,"lorenz_sigma%.0lf_r%.0lf_b%.2lf.dat",sigma,r,b);
 escribe=fopen(nombre, "wt");

for(i=0;i<N;i++)
    {
    x_old[i]=1.0;
    }        //c.i.
t=0;



    for(i=0;i<N;i++)
    {

	k_1[i]=0;
	k_2[i]=0;
	k_3[i]=0;
	k_4[i]=0;
    }
  
	
  
while(t<=t_final)
  {


  for(i=0;i<N;i++)
     {
	
     k_1[i]=funcion(x_old[i])*delta_t;


     k_2[i]=funcion(x_old[i]+0.5*k_1[i])*delta_t;

	
     k_3[i]=funcion(x_old[i]+0.5*k_2[i])*delta_t;


     k_4[i]=funcion(x_old[i]+k_3[i])*delta_t;
     

     x_old[i]=x_old[i] + (1.0/6.0)*(k_1[i]+2*k_2[i]+2*k_3[i]+k_4[i]);

    
     fprintf(escribe,"%lf  ",x_old[i]); 
     }
      //ojo con la fraccion, pq 1/6=0!!!!!

 fprintf(escribe,"\n"); 

   t+=delta_t;



 }            //fin del algoritmo



 fclose(escribe);
}


double funcion(double x)             //corresponde a las ec. de lorenz
{
  
double resultado;

 if(i==0)     //coord. x
 {                    
 resultado=sigma*(x_old[1]-x_old[0]);
 return (resultado); 
 }
 else 
  {
     if(i==1)      //coord. y
     {
     resultado=r*x_old[0]-x_old[1]-x_old[0]*x_old[2];
     return (resultado); 
     }
     else
      {
	 if(i==2)      //coord. z
	 {
	  resultado=x_old[0]*x_old[1]-b*x_old[2];
	  return (resultado); 
	 }
      }
  }
 


}
