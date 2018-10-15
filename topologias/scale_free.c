
        /*Creación de una red Scale-Free por el método de Barabasi y Albert*/


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# define n_ini 3
# define N_tot 1000
# define num_links 3
#define iter 20         //repetir el proceso para tener estadística suficiente


//Generación de números aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//número aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];

void inicia_rand(int semilla);

	int k[N_tot+1], conectado[num_links+1], M[N_tot+1][num_links+1], P[N_tot+1];

main()
{

	int i, j, n, cont, x, sum;
	float  r;
	double PK[N_tot],PK_tot[N_tot];

FILE *escribe;
FILE *escribe2;

	if( (escribe=fopen("conectividad_correg_estad.txt","wt")) == NULL)
    	{
        printf("error en la apertura del archivo conectividad.txt");
        getchar();
        exit(0);
        }
     
    if( (escribe2=fopen("datos_correg_estad.txt","wt")) == NULL)
   	{
        printf("error en la apertura del archivo datos.txt");
        getchar();
        exit(0);
        }

inicia_rand(time(0));

     for(i = 0; i < N_tot; i++)
		PK_tot[i] = 0;

for (n=0 ; n<iter; n++)              //bucle para hacer estadística
    {
	for(i = 0; i <= num_links; i++)
		conectado[i] = 0;
    
	/* inicializo k, M y P, para n_ini nodos conectados entre sí */
	
	for(i = 0; i <=N_tot; i++)
	  {
	    k[i] = 0;
	    P[i] = 0;
	    for(j = 0; j <= num_links; j++)
	      M[i][j]=0;
	  }


	for(i = 1; i <= n_ini; i++)         //los n_ini primeros nodos conectados entre si
	  {
	    k[i] = n_ini - 1;
	    for(j = 1; j < i; j++)
	      M[i][j] = j;    
	  }
    	
	for(i = 1; i <= N_tot; i++) //las prob  (nulas salvo para los n_ini nodos)
	  P[i] = P[i-1] + k[i];
	
	
	/* lanzo los m links */
	
   for(i = n_ini+1; i <= N_tot; i++)      //para los siguientes nodos, a partir del n_ini+1
	  {
	    for(j = 1; j <= num_links; j++)
	      {
		r = FRANDOM;
		r = r*P[N_tot];
		
		for (cont = 1; ((cont < N_tot) && (P[cont] <= r)); cont++);     //busco el primer nodo que cumple r<Prob(i)
		conectado[j] = cont;

		for(x = conectado[j]; x <= N_tot; x++)     //quito el nodo i para no volver a contarlo
		  P[x] = P[x] - k[conectado[j]];
	      }
      

   	    for(j = 1; j <= num_links; j++)        //reescribo el vector k
	      k[conectado[j]]++;
	    

	    k[i] = num_links;
			

	    for(j = 1; j <= num_links; j++)       //reescribo el vector p
	      {
		for(cont = conectado[j]; cont <= N_tot; cont++)
		  P[cont]+= k[conectado[j]];
	      }
	    
	    for(cont = i; cont <= N_tot; cont++)
	      P[cont] = P[i] + num_links;
	    

   	    for(j = 1; j <= num_links; j++)                //reescribo la matriz M
	      M[n_ini + i][j] = conectado[j];
	    
	  }		
		
	printf("salgo del bucle en el que hago los links\n");

	/* Escribo las conectividades */
	for(i = 0; i <= N_tot; i++)
		fprintf(escribe2, "%d %d\n", i, k[i]);
	
	/* Calculo la distribución de conectividad */
	/* primero inicializo el vector PK a 0 */
	
	for(i = 0; i < N_tot; i++)
		PK[i] = 0;
	
	for(i =1; i <= N_tot; i++)
        PK[k[i]]++;

  	

   	
	/* Normalizo PK */
  	for( i = 1; i < N_tot; i++)
   		{
          PK[i] = PK[i]/N_tot;
          PK_tot[i]+=PK[i];
        }



 }     //fin del bucle en iter


	/* Normalizo PK_tot */
  	for( i = 1; i < N_tot; i++)
   		PK_tot[i] = PK_tot[i]/iter;


	/* Escribo k y PK en un archivo */

 for(i = num_links; i <= N_tot - 1; i++)
		fprintf(escribe, "%d %lf\n", i, PK_tot[i]);
		
 	fclose(escribe);
	fclose(escribe2);
 	
 	exit(1);

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
