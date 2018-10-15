// Programa para hacer el diagrama de error de una serie temporal correspondiente a las 
// variaciones relativas de la bolsa, a partir de un fichero ya filtrado con un umbral, y convertido en
// 0's  y 1's. tb hace distribucion de probabilidad de intervalos entre sucesos.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>



#define filas 2000

#define landa 0.01   //la inv de la media de la distrib de prob de la serie temporal    NO ES LA INVERSA DE LA ANCHURA?!!!1



//Generacion de numeros aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];




int i,j;
int V[filas+1];  //serie  temporal 
double p_dt[filas+1], norma; // histograma de tiempos entre eventos
int dt_max,dt_min,interv,tau,tau_ini, tau_fin;
int flat;
double fa, fe;
int n_sucesos,cont;

int t, aux;
double v,w;
double dobleprec,dobleprec1,dobleprec2,dobleprec3;

void inicia_rand(int semilla);




int main()
{ 
    FILE *fich1;
    FILE *fich2;
    FILE *fich3;
    FILE *fich4;

    char file1[256],file2[256],file3[256],file4[256];
			  
			  
    inicia_rand(time(0));
			  
 
			  
			 
			  
  sprintf(file1,"serie_temporal_gauss.dat");  
  fich1=fopen(file1,"wt");
  fclose(fich1);


			  
  /* sprintf(file2,"TCvsK0.00-bini%.3lf_binning25.dat",b);   //falta por modificar!!	  
  fich2=fopen(file2,"wt");
  fclose(fich2);*/

  sprintf(file3,"p_dt_gauss_%d_%.2lf.dat",filas,landa);
  fich3=fopen(file3,"wt");
  fclose(fich3);


  sprintf(file4,"diagrama_error_gauss_%d_%.2lf.dat",filas,landa);
  fich4=fopen(file4,"wt");
  fclose(fich4);






  for(i=0;i<=filas;i++)
    {
	V[i]=0;
	p_dt[i]=0;
    }



  //////////   creo la serie temporal de prueba       /////////////

  t=0;
  while(t<=filas)
  {
      w=FRANDOM;       //genero un numero aleat entre (0,1)        
       
      v=-log(1-w);   //para convertir la distribuc plana en una poisson     (ojo!: log es el NEPERIANO)
      v=v/landa;
                  
      t+=v;         //lo convierto a entero   y lo sumo a la posicion anterior
      V[t]=1;

      
  }


  fich1=fopen(file1,"at");
  for (i=1;i<=filas;i++)
  {
      //  printf("%d\n",V[i]);
      fprintf(fich1,"%d \n",V[i]);     
  }
  
  fclose(fich1);
  

    
////////////////////////////////





/*            //  tomo la serie temporal de un archivo

  fich2=open(file2,"r");

  for(i=1;i<=filas;i++)
    {
      fscanf(fich2,"%d\n",&V[i]);
    }
  
  fclose(fich2);
  
  

   for (i=1;i<=filas;i++)
    {
      printf("%d\n",V[i]);
    }
  
*/ 
  
  
  //calculo dt_max y dt_min
  
  dt_max=0;
  dt_min=filas;
                         //EL INTERVALO ENTRE DOS SUCESOS CONSECUTIVOS ES UNO!!!
  interv=0;     
  for(i=1;i<=filas;i++)
    {
	if(V[i]==0)    //sin no hay suceso
	{
	    interv++;	  	 
	}
      else      //   V[i]==1  si lo hay
      {	  

	  if(interv > dt_max)
	  {
	      dt_max=interv;
	      
	      //printf("encontrado mejor max %d\n", k_max);
	  }
	  if( interv < dt_min)       // (        && (interv > 0))
	  {
	      dt_min=interv;
	      
	      //printf("encontrado mejor max %d\n", k_max);
	  }


	  p_dt[interv]++;   //recuento para el histograma de intervalos entre sucesos
	  norma++;    //////normalizo sobre el numero de intervalos 

	  interv=0;       //EL INTERVALO ENTRE DOS SUCESOS CONSECUTIVOS ES CERO!!!   por eso
      }
    }
 
  printf("dt_max: %d  dt_min: %d\n",dt_max,dt_min);

  tau_ini=dt_min;      
  tau_fin=dt_max;

  // getchar();

  dobleprec=norma;
  for(i=0;i<=filas;i++)
    {
	p_dt[i]= p_dt[i]/dobleprec;   ////////normalizo sobre el numero de intervalos 
    }
  
 
  
  
  fich3=fopen(file3,"wt");      //escribo la p(dt)

  for(i=1;i<=filas;i++)
    {    
	fprintf(fich3,"%d  %lf\n",i,p_dt[i]);   
	//printf("%d  %lf\n",i,p_dt[i]);  
    }
  
  fclose(fich3);



 
  tau=0;//tau_ini;
  // printf("tau_ini:%d  tau_fin:%d\n",tau_ini,tau_fin);  
  //getchar();
  
			  
  while (tau<=dt_max+1)       //bucle para obtener los distintos puntos del diagramas de error
  {
      fa=fe=0.;
      n_sucesos=0;
      cont=1;   
      flat=0;      //alarma en off

      for(i=1;i<=filas;i++)
      {
	 
	  if((cont==tau) && (tau>0))     //casos habituales
	  {
	      flat=1;      //pongo la alarma
	  }

	  if(tau==0)     //parche para el caso especial
	  {
	      flat=1;      //pongo la alarma siempre
	  }

	  
	  if(V[i]==0)     //si no ocurre un suceso 
	  {
	      if(flat==1)    //y la alarma estaba puesta
	      {
		  fa++;
	      }
	  }

	  if(V[i]==1)     //si ocurre un suceso 
	  {
	      n_sucesos++;
	      
	      
	      if(flat==1)       //ojo!!!!!: el tiempo de alarma tambien hay que sumarlo en este caso!!!!!!!! 
	      {
		  fa++;  
	      }

	      if(flat==0)   //y la alarma no estaba puesta
	      {
		  fe++;
	      }

	      flat=0;    //despues de un suceso, quito la alarma
	      cont=0;   //y pongo a cero el contador de tiempo a esperar antes de conectar la alarma  
	      if(tau==0)
	      {
		  flat=1;       //parche para el caso especial: no espero ni un paso de tiempo
	      }
	  }

	 cont++;
	  
      }                 

      
      // printf("fa:%d  (filas:%d)   fe:%d  (n_sucesos:%d)\n",fa,filas,fe,n_sucesos);
      
      fa=fa/filas;       //normalizo al tiempo total
      fe=fe/n_sucesos;
      
      
      printf("fe:%lf    fa:%lf   tau:%d\n\n",fe,fa,tau);
      
// getchar();
      
//escribo el punto del diagrama de error
      
      fich4=fopen(file4,"at");     
      fprintf(fich4,"%lf   %lf    %d\n",fe,fa,tau);     
      fclose(fich4);
      
      
      tau++;
  }
  







  
  
  
}



/////////////////////////////////////////////////
/////////////////////////////////////////////////


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



/////////////////////////////////////////////
/////////////////////////////////////////////   




  
 
