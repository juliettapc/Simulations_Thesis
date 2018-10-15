		//programa para implementar el metodo de Euler de resolucion
		//de ecuaciones diferenciales

# include <stdio.h>
# include <stdlib.h>
# include <math.h>


double delta_t, x_old, t_final, t, t_o;
double funcion(double);


int main()
{
/*printf("introduce la condicion inicial: (t x) y  el paso de tiempo\n");
scanf("%lf %lf %lf",&x_old,&t_o,&delta_t);
printf("introduce el tiempo para el que quieres la solucion\n");
scanf("%lf ",&t_final);
*/

x_old=1;     //c.i.


delta_t=t=0.0001;
t_final=1;

printf("t_final=%lf  delta_t=%lf\n",t_final,delta_t);

while(t<=t_final)
 {
 x_old=x_old + funcion(x_old)*delta_t;
 printf("x(t=%lf)=%lf\n",t,x_old);
 t+=delta_t;

 }

printf("el resultado es\n x(%lf)=%lf\n",t,x_old );
getch();
}

double funcion(double x)
{
double valor;                    //corresponde a la ec.dif.: f(x)=dx/dt=-x
valor=-x;
printf("x=%lf f(x)=%lf\n",x,valor);
return (valor);
}
