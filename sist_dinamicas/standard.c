
                      //simulacion numerica del map standard

# include <stdio.h>
# include <stdlib.h>
# include <math.h>


#define Npuntos 100

//Generación de números aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])   //número aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];

char nombre [256];


void inicia_rand(int semilla);

FILE *escribe;

double p,q,k;
double t, delta_t, t_fin;
double x,PI;
int n;

int main()
{
inicia_rand(time(0));

    k=0.3;


    sprintf(nombre,"standard_k%.2lf.dat",k);
    escribe=fopen(nombre,"wt");

    PI=3.141592653589793;
    delta_t=0.1;
    t_fin=500;

    for (n=1;n<=Npuntos;n++)
    {    
	x=FRANDOM*2.0*PI;
	p=x;
	x=FRANDOM*2.0*PI;
	q=x;
	//printf ("%lf  %lf\n",q,p);
	fprintf(escribe,"%lf  %lf \n",p,q);
  
	t=0;
    while(t<=t_fin)
    {
	p=p+k*sin(q);
	q=q+p;
	while (p>2*PI)
	{
	    p-=2*PI;
	}
        while (q>2*PI)
	{
	    q-=2*PI;
	}
        while (p<0)
	{
	    p+=2*PI;
	}
        while (q<0)
	{
	    q+=2*PI;
	}
	fprintf(escribe,"%lf  %lf\n",q,p);
	//printf ("%lf  %lf\n",q,p);
	t+=delta_t;
    }
  }

fclose(escribe);

}





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
