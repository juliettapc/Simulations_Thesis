     /*Implementacion del Dilema del Prisionero en una red gardeñes-moreno*/
                   //contiene ademas un bucle de estadística y otro para barrer en b y obtener el n_coop en le equil. (==despues de muchas jugadas)


# include <stdio.h>
# include <stdlib.h>
# include <math.h>


# define N 4000 //tamaño de la red
# define m_o  2      //nodos inicialmente unidos
# define m   2      //nodos nuevos añadidos a cada paso de tiempo
# define alfa 1.0

# define jugadas 30000  //pasos de tiempo para el juego
# define Niter 10   //estadistica

# define K_MAX   300          //para el tamaño de C[N][K_MAX]

//Generación de números aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])   //número aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
char nombre [256],nombre3 [256];

 
int k[N+1],k_PA[N+1],A[N+1];                    //conectividad topológica y de Preferential Attachment (PA)
int   M[N+1][m+1], unido[m+1],Cuidado[N+1];              //  unido[] guarda a quién se ha unido i , y cuidado[] quien le ha lanzado ya un link
double P[N+1],P_prov[N+1];        //P_prov[] es la provisional para manipularla
double PK[N+1],PK_tot[N+1];        //  para la distribucion P(k)
double r,v;                                  //para guardar los aleatorios
int  i, j, jj, w,g,d, q, x[N+1],y,z, C[N+1][K_MAX+1],n, ValorA,steps,s,flat;
int norma, norma2, norma_aleat,tipo[m_o+1];
double dnorma_aleat;
int si[N+1], cont,cont2;

int n_coop, e[N+1],  w, t, k_m;       // variables para el dilema
double  ben[N+1],p;

double b, b_max, b_min, delta_b, c_media;                 //barrer en b y estadistica
int iter;
	
void inicia_rand(int semilla);
void construir_red();
void histograma();
void comprobar_red();
void juego();

FILE *escribe;
FILE *escribe2;
FILE *escribe3;


int main()
{


b_min=1.0;
b_max=2.0;
delta_b=0.1;

/*escribe=fopen(nombre,"wt");
escribe2=fopen("comprobar.dat","wt");*/
sprintf(nombre3,"c_med_frente_b%.2lf-%.2lf_N%d_alfa%.3f.dat",b_min,b_max,N,alfa);

//sprintf(nombre3,"tiempo_gm_N%d_alfa%.3f.dat",N,alfa);
escribe3=fopen(nombre3,"wt");
//fclose(escribe3);    


printf("Red de N=%d, con alfa=%lf \n",N,alfa);
printf("jugando desde: b_min=%lf hasta: b_max=%lf \n",b_min,b_max);

inicia_rand(time(0));



b=b_min;
while(b<=b_max)
{
 printf("b=%lf\n",b);

 c_media=0;
 for(iter=1;iter<=Niter;iter++)
  {
   printf("iteracion:%d\n",iter);
   construir_red();
    //comprobar_red();
    //histograma();
   juego();
   c_media+=n_coop;
  }
  
 c_media=c_media/(Niter*N);          //normalizo y tanto por uno

printf("c_media:%lf\n",c_media);


fprintf(escribe3,"%lf  %lf\n", b,c_media);
fflush(escribe3);                  // bolcar el buffer

 b=b+delta_b;

}

fclose(escribe3);
//fclose(escribe2);
//fclose(escribe);


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



void construir_red()
{

 ValorA=3;
 steps=N-m_o;

 for(i = 0; i <=N ; i++)
     PK_tot[i] = 0;


	for(i = 0; i<=N; i++)
	  {
	    k[i] = 0;
            k_PA[i] = 0;
	    A[i] = 0;

        for(j = 0; j <= m_o; j++)
	        M[i][j] =0;
      }


    for(j = 1; j <=m_o ; j++)                    // inicializo los m_o primeros nodos conectados entre si
	   A[j]=ValorA;

    for(i = 2; i <=m_o ; i++)
	   {
	    for(j = 1; j <=i-1 ; j++)
	      {
           M[i][j] = j;
           //printf("M[%d][%d] = %d \n",i,j,M[i][j]);
          }
	   }


	for(i = 0; i <= N; i++)
	  {
      P[i] = 0;
      P_prov[i]=0;
      }

   for(i = 1; i<= N; i++)                    //prob iniciales
      {
       for(j=1; j<=i; j++)
	    { P[i] = P[i] + k_PA[j] + A[j];}

       P_prov[i]=P[i];
      }

   norma=P[N];
   norma2=norma;

    for(s = 1; s <= steps; s++)    //bucle a todos los  demas nodos (paso de tiempo)
      {

      //printf("\nnodo: %d \n",s+m);

      for(i=1;i<=N;i++)
        si[i]=i;
      for(i=m+s;i<N;i++)   //me quito yo
       si[i]=si[i+1];

      si[N]=0;


      for(i=0; i<=m; i++)
         {unido[i]=0;}

      for(i=1; i<=N; i++)
         {P_prov[i]=P[i];}

      norma2=norma;
      P_prov[m_o + s]=0.;

      for(i=m_o + s + 1; i<=N; i++)
         {P_prov[i]=P_prov[i] - (k_PA[m_o+s]+A[m_o+s]); }

      norma2=norma2-(k_PA[m_o+s]+A[m_o+s]);


      //Quitar los posibles enlaces lanzados a mi anteriormente

      jj=0;
      cont2=0;
      for(q = 1; q<m+s; q++)      //yo soy el nodo m+s
	    {

	    for(i=1;i<=m;i++)
	      {
	      if(M[q][i]==m+s)
		   {
		   jj++;
		   Cuidado[jj]=q;     //guarda los que ya me han elegido antes
		   P_prov[q]=0.;

		  for(j=q+1;j<=N;j++)
		    { 
		      P_prov[j]=P_prov[j] - (k_PA[q]+A[q]); 
		    }

		  norma2=norma2-(k_PA[q]+A[q]);
		  }
	    }
	  }

  cont=0;
  for(i=1;i<=jj;i++)         //quito los nodos que me han lanzado un link
    {
    d=Cuidado[i];

    for(j=d-cont;j<N;j++)
      {
      si[j]=si[j+1];
      }
      cont++;
    }

   cont2=cont;
   for(q = 1; q<=m; q++)        //bucle a los m links nuevos que voy a añadir
       {
         r = FRANDOM;

         if (r>alfa)
            {flat=1;}
         else
            {flat=0;}

         if (flat==1)    //scale-free
            {
            tipo[q]=1;
            v=FRANDOM*norma2;
            for(j=1; j<=N; j++)
               {
               if(v<P_prov[j])        //elijo el nodo
                 {
                  unido[q]=j;
                  break;
                 }
               }
             g=unido[q];
             P_prov[g]=0;    //lo guardo y anulo su prob para no volver a cogerlo

             for(j=g+1; j<=N; j++)
	            {P_prov[j] = P_prov[j] -( k_PA[g] + A[g]);}

             norma2=norma2-(k_PA[g]+A[g]);


             for(i=g-cont2;i<N;i++)
               si[i]=si[i+1];


            }

         else   //aleatoria
            {
             tipo[q]=0;
             dnorma_aleat=N-jj-1-q;
	         norma_aleat=dnorma_aleat;

	         r=FRANDOM;
	         r=r*dnorma_aleat;
	         w=(int)r+1;

             unido[q]=si[w];     //me uno al aleatorio corresp. del los posibles
             g=si[w];     //lo guardo y anulo su prob para no volver a cogerlo
             P_prov[g]=0;

             for(j=g+1; j<=N; j++)
                {P_prov[j] = P_prov[j] -( k_PA[g] + A[g]);}

             norma2=norma2-(k_PA[g]+A[g]);

             //actualizo:
             for(i=w;i<N;i++)
                si[i]=si[i+1];
           }

        cont2++;

       }    //fin del bucle sobre los m links lanzados

       for(i=1; i<=m; i++)     //bucle sobre los m nuevos links que acabo de lanzar
          {
           g=unido[i];

           if(A[g]==ValorA)   //si ya habia recibido algun link  antes
             {
             k[g]++;   //aumento su conectividad topológica

             if(tipo[i]==1)  //link por PA
               {
               k_PA[i]++;
               for(j=g; j<=N; j++)
                 {P[j]++;}
               norma++;
               }
             }

          else     //si no habia recibido links antes
             {
             A[g]=ValorA;
             k[g]++;

             if(tipo[i]==1)  //link por PA
               {
               k_PA[i]++;

               for(j=g; j<=N; j++)
                 {P[j]=P[j]+ValorA+1;}
               norma=norma+ValorA+1;
               }
             else
              {
               if(tipo[i]==0)  //link aleatorio
                 {
                 k_PA[g]=k_PA[g];

                 for(j=g; j<=N; j++)
                   {P[j]=P[j]+ValorA;}
                 norma=norma+ValorA;
                 }
               }
             }
          M[m+s][i]=g;

          }       //fin dle bucle sobre los m links nuevos


       if(A[m+s]==0)
         {
         A[m+s]=ValorA;

         for(j=m+s; j<=N; j++)
            {P[j]=P[j]+A[m+s];}
         }

      norma=P[N];


 }  

   //fin del bucle a los nodos de la red
		
    for(i=1;i<=m_o;i++)
    {k[i]=k[i]+m-1;}
    
    for(i=m_o+1; i<=N; i++)
    {k[i]=k[i]+ m;}


 //obtencion de la matriz C a partir de la M


     for(i=0;i<=N;i++)
       {
	   x[i]=0;
	   
	   for(j=0;j<=N;j++)
	   {
	       C[i][j]=0;
	   }
       }

     for(i=1;i<=N;i++)
        {
	    for(j=1;j<=m;j++)
            { 
		if(M[i][j]!=0)
		{
		    x[i]++;
		    z=M[i][j];
		    //printf("%i  %i\n",i,z);
		    C[i][x[i]]=z;
		    x[z]++;
		    C[z][x[z]]=i;	 
		}
	    }
        }

//printf("red construida\n");

}  //fin de la funcion "red()"



void comprobar_red()
{

 for(i=1;i<=N;i++)
    {
    for(j=1;j<=N;j++)
       {
       if(C[i][j]==i)       //si esta conectado a si mismo
          {
          fprintf(escribe2,"error: conexion consigo mismo, nodo %d\n", i );
          }
       }

    for(j=2;j<=k[i];j++)
       {
       for(y=1; y<=j-1; y++)
         {
          if(C[i][j-y]!=0  &&  C[i][j-y]==C[i][j])          // si existen enlaces dobles
             {
              fprintf(escribe2,"error: enlace doble, nodo %d\n", i );
             }
         }
       }
    }

printf("red comprobada\n");
}         //fin de la funcion "comprobar_red()"



void histograma()            //costruccion de P(k)
{

 for(i = 0; i <= N; i++)           //inicializo
	PK[i] = 0;
	
 for(i =1; i <= N; i++)            //recuento
    PK[k[i]]++;

 for(i = 1; i <= N; i++)         //normalizo
	{
	 //PK[i] = PK[i]/N;
     PK_tot[i]+=PK[i];
     }


	/* Escribo k y PK en un archivo */

for(i = 1; i<N; i++)
    {fprintf(escribe,"%d  %lf\n", i, PK[i]);}



printf("histograma hecho\n");

}          //fin de la funcion "histograma()"




void juego()            //implementacion del dilema del prisionero:
{

//printf("the game begins \n");

n_coop = p = k_m = 0;


for(i=1;i<=N;i++)
   {
    ben[i]=0;
    e[i]=1;
   }

/*for(i=1;i<=N/2;i++)
   {
   r=FRANDOM;
   r=r*N;
   w=(int)r+1;
   e[w]=0;
   n_coop++;
   }
*/

for(i=1;i<=N;i++)           // establezco inicialmente el 50% de cooperantes
   {
   r=FRANDOM;

   if(r<0.5)
     {
     e[i]=0;
     n_coop++;
     }
   }



for(n=1;n<=jugadas;n++)            //pasos de tiempo
   {
    n_coop=0;

    for(i=1;i<=N;i++)          // cada nodo juega 1partida con cada vecino
    {
	ben[i]=0;
    }
    
    for(i=1;i<=N;i++)          // cada nodo juega 1partida con cada vecino
    {
	for(j=1;j<=k[i];j++)
         {
          y=C[i][j];            //para no repetir
          if(y>i)
          {  
          
	     if(e[i]==0 && e[y]==0)
	     {
		 ben[i]++;
		 ben[y]++;
	     
		 //printf("%d : %d  empate coop \n", i,y);
	     }

	     if(e[i]==0 && e[y]==1)
	     {
		 //ben[i]=ben[i];
		 ben[y]=ben[y]+b;
		 //printf("%d : %d  coop contra defect\n", i,y);
	     }

	     if(e[i]==1 && e[y]==1)
	     {
		 // ben[i]=ben[i];
		 //ben[y]=ben[y];
		 //printf("%d : %d  empate defect \n", i,y);
	     }

	     if(e[i]==1 && e[y]==0)
	     {
		 ben[i]=ben[i]+b;
		 //ben[y]=ben[y];
		 //printf("%d : %d  empate defect \n", i,y);
	     }
	 }

         }
    }

  for(i=1;i<=N;i++)        //comparo estrategias
     {
       r=FRANDOM;
       r=r*k[i];
       w=(int)r+1;

    
       t=C[i][w];    //el vecino elegido al azar
       //if(t<1)
       //{printf("ERROR-----------------%i  %i  %i\n",i,w,k[i]);}

       if(k[i]>k[t])
       {k_m=k[i];}
       else
       {k_m=k[t];}

       if(ben[i]<ben[t])
       {
	   p=(ben[t]-ben[i]);
	   p=p/(b*(double)k_m);

          //printf("p= %f  \n",p);
          v=FRANDOM;
          if(v<p)
	  {
	      e[i]=e[t];
	      //printf("cambio a estrategia: %d  \n",e[t]);
	  }
       }

      }

    for(i=1;i<=N;i++)     //recuento de cooperantes
      {
      if(e[i]==0)
        n_coop++;
      }
   //printf("num coop.: %d  \n",  n_coop);

   }


 
}
