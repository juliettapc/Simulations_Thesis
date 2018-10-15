
        /*Craación de una red aleatoria*/
        /* y calculo del " average path length"*/


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# define p 0.004        //probabilidad de conexion de una pareja cualquiera de nodos
# define N 1000      //tamaño de la red
# define iter 1    //estadistica


//Generación de números aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//número aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
	

void inicia_rand(int semilla);
 
 int k[N+1], D[N+1][N+1], C[N+1][N+1], M[N+1][N+1], k[N+1];
 double r, PK[N+1],PK_tot[N+1];;
char file32[256];
	
int main()
{
	
	int number, i, j, h, l, ultimo, label[N+1], analizar[N+1], cont, n, m;
	float D_med;
	FILE *escribe;
	FILE *escribe2;
	FILE *fich32;    
	if( (escribe2=fopen("distancia_correg_rand_estad.txt","wt")) == NULL)
    		{
	        printf("error en la apertura del archivo distancia.txt");
	        getchar();
	        exit(0);
	        }
	
	if( (escribe=fopen("PK_correg_rand_estad.txt","wt")) == NULL)
    		{
	        printf("error en la apertura del archivo PK.txt");
	        getchar();
	        exit(0);
	        }

D_med=0.0;

  for(i = 0; i <= N ; i++)
		PK_tot[i] = 0;


 for (n=0;n<iter;n++)     //bucle para hacer estadística
 {




      sprintf(file32,"%dNet_p%f.NET",N,p);   
      fich32=fopen(file32,"wt");            //pintar las redes
      fclose(fich32);
      
	    





	/* Inicializo C */
	for(i = 0; i <= N; ++i)
	  {
	    k[i]=0;
	    for(j = 0; j <= N; ++j)
	      {
		C[i][j]=0;        //matriz de conectividad: la fila i contiene todos los nodos a los que esta unido el  nodo i
	      }
	  }	
	
	/* Construyo la red */
	
    inicia_rand(time(0));

	for(i = 1; i <= N; i++)    //bucle a todos los nodos
	  {
	    ultimo = 1;
	    for(j = i + 1; j <= N; j++)
	      {
		r = FRANDOM;

		   if(r < p)
		     {
   		      k[i]++;
		      k[j]++;
		      C[i][k[i]]=j;
		      C[j][k[j]]=i;
		     }
	       }
	   }




	printf("He construido la red\n");
	
	/* Calculo la distribución de conectividad */
	/* primero inicializo el vector PK a 0 */
	
	for(i = 0; i < N; i++)
		PK[i] = 0;
	
	
	for(j = 1; j <= N; j++)
        PK[k[j]]++;


		  
   	
   	/* Normalizo PK */
  	for( i = 1; i <=N; i++)
	  {
       PK[i] = PK[i]/N;
       PK_tot[i]+=PK[i];
       }
   	


	  

	printf("He calculado la distribucion de conectividad\n");
	
	/* Escribo k y PK en un archivo */

	for(i = 1; i<=N; i++)
	  fprintf(escribe, "%d %lf\n", i, PK[i]);     
		
		      

	/* Calculo el average path length */
	for(i = 1; i <= N; i++)
		{
		for(j = 1; j <= N; j++)
			label[j] = 2*N;

		label[i] = 0;
		number = 1;
		analizar[1] = i;

		for(j = 1; j <= number; j++)
		  {
		    for(h = 1; h <= k[analizar[j]]; h++)
		      {
			if(label[C[analizar[j]][h]] > N+1)
			  {
			    number++;
			    analizar[number] = C[analizar[j]][h];
			    label[C[analizar[j]][h]] = label[analizar[j]] + 1;
			  }
		      }
		  }

		for(j = 1; j <= N; j++)
		  {
		    D[i][j] = label[j];
		  }
		}
		



	for(i = 1; i <= N; i++)
	  {
	    for(j = i+1; j <=N; j++)
	      D_med+= D[i][j];
	  }



  
      fich32=fopen(file32,"at");                                //guardo los pares de links para pintar la red
      fprintf(fich32,"*Vertices      %d\n*Edges\n",N);      
      fclose (fich32);
      
      fich32=fopen(file32,"at");
      for(i = 1; i<N; i++)         
	{
	  for(j = 1; j<N; j++)         
	    {
	      if(C[i][j]!=0)
		{
		  fprintf(fich32,"%d     %d\n",i,C[i][j]);      
		}
	    }
	
	}
      fclose (fich32);      








}  //fin del bucle en iter


/* Normalizo PK_tot */
  	for( i = 1; i <= N; i++)
   		PK_tot[i] = PK_tot[i]/iter;


 for(i = 0; i <= N; i++)
		fprintf(escribe, "%d %lf\n", i, PK_tot[i]);

	D_med = (2.*D_med)/((N+1)*N*iter);


	
	printf("He calculado el camino libre medio\n");
	
    fprintf(escribe2, "%f\n", D_med);
	printf("D: %f\n", D_med);
	


	fclose(escribe);
	fclose(escribe2);
	//    getch();
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
		
