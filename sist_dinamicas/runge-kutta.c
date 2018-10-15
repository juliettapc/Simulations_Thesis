			//programa para implementar el metodo de Runge-Kutta de resolucion
			//de ecuaciones diferenciales (en una dimension)


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
//#define valor_real 0,367879441171442    // valor de 1/e


double x_old,  k_1, k_2, k_3, k_4,esto;
double delta_t, t_final, t, t_o;
double funcion(double);


int main()
{
/*printf("introduce la condicion inicial: (t x) y  el paso de tiempo\n");
scanf("%lf %lf %lf",&x_old,&t_o,&delta_t);
printf("introduce el tiempo para el que quieres la solucion\n");
scanf("%lf ",&t_final);
*/

x_old=1;     //c.i.
t=0;


delta_t=0.001;
t_final=1;

k_1=0;
k_2=0;
k_3=0;
k_4=0;

printf("t_final=%lf  delta_t=%lf\n",t_final,delta_t);

while(t<=t_final)
 {

  printf("k_1=%lf k_2=%lf k_3=%lf k_4=%lf, t=%lf\n\n",k_1, k_2, k_3, k_4, t);
 
  k_1=funcion(x_old)*delta_t;
  k_2=funcion(x_old+0.5*k_1)*delta_t;
  k_3=funcion(x_old+0.5*k_2)*delta_t;
  k_4=funcion(x_old+k_3)*delta_t;

  printf("x(t=%lf)=%lf\n",t,x_old);

  x_old=x_old + (1.0/6.0)*(k_1+2*k_2+2*k_3+k_4);

      //ojo con la fraccion, pq 1/6=0!!!!!

  t+=delta_t;


  printf("k_1=%lf k_2=%lf k_3=%lf k_4=%lf, t=%lf\n",k_1, k_2, k_3, k_4, t);
 }  

printf("\nel resultado es\n x(%lf)=%lf\n",t,x_old);
//getch();
}

double funcion(double x)
{
double valor;                    //corresponde a la ec.dif.: f(x)=dx/dt=-x
valor=-x;
//printf("x=%lf f(x)=%lf\n",x,valor);
return (valor);
}
