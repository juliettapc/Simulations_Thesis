        /*Craación de una red gardeñes-moreno*/


# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# define N 50    //tamaño de la red
# define m_o  3      //nodos inicialmente unidos
# define m   3      //nodos nuevos añadidos a cada paso de tiempo
# define iter 100     //estadistic
# define alfa 0.0

# define K_MAX   300          //para el tamaño de C[N][K_MAX]

//Generación de números aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//número aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
char nombre [256];


int k[N+1],k_PA[N+1],A[N+1];                    //conectividad topológica y de Preferential Attachment (PA)
int   M[N+1][m+1], unido[m+1],Cuidado[N+1];              //  unido[] guarda a quién se ha unido i , y cuidado[] quien le ha lanzado ya un link
double P[N+1],P_prov[N+1];        //P_prov[] es la provisional para manipularla
double PK[N+1],PK_tot[N+1];        //  para la distribucion P(k)
double r,v;                                  //para guardar los aleatorios
int  i, j, jj, w,g,d, q,x[N+1],y,z, C[N+1][K_MAX+1],n, ValorA,steps,s,flat;
int norma, norma2, norma_aleat,tipo[m_o+1];
double dnorma_aleat;
int si[N+1], cont,cont2;
	
void inicia_rand(int semilla);


int main()
{
	
FILE *escribe;
FILE *escribe2;

sprintf(nombre,"pk_gm_N%d_alfa%.3f.dat",N,alfa);
escribe=fopen(nombre,"wt");
escribe2=fopen("comprobar.dat","wt");
 


ValorA=m;
steps=N-m_o;

for(i = 0; i <=N ; i++)
    PK_tot[i] = 0;


for (n=0; n < iter; n++)     //bucle para hacer estadística

  {
    printf("iteracion numero %d \n",n);
    fprintf(escribe2, "iteracion numero %d \n",n);
    
    inicia_rand(time(0));


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
    //printf("Cuidado[%d]=%d\n",i,d);

    for(j=d-cont;j<N;j++)
      {
      //printf("j=%d\n",j);
      si[j]=si[j+1];
      //printf("desplazando el %d al puesto del %d: %d -> %d \n",j+1,j,si[j],si[j+1]);
      }
      cont++;
    }

    /*printf("\ninicialmente puedo unirme a %d nodos:\n ",N-jj-1);  //todos menos yo y aquellos que me hayan cogido
     for(j=1;j<=N;j++)
        printf("%d  ",si[j]);
     printf("\n \n");
     getch();*/

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
             //printf("me uno (p-a) al nodo %d\n",g);
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
             dnorma_aleat=N-jj-1-q;           //esta bien?????
	         norma_aleat=dnorma_aleat;

	         r=FRANDOM;
	         r=r*dnorma_aleat;
	         w=(int)r+1;
             //printf("aleatorio w=%d\n",w);

             unido[q]=si[w];     //me uno al aleatorio corresp. del los posibles
             g=si[w];     //lo guardo y anulo su prob para no volver a cogerlo
             P_prov[g]=0;
             //printf("me uno (aleat) al nodo %d\n\n",g);

             for(j=g+1; j<=N; j++)
                {P_prov[j] = P_prov[j] -( k_PA[g] + A[g]);}

             norma2=norma2-(k_PA[g]+A[g]);

             //actualizo:
             for(i=w;i<N;i++)       // i=w esta bien???
                si[i]=si[i+1];
           }

        cont2++;
       /* printf("ahora puedo a:\n ");
         for(i=1;i<=N;i++)
            printf("%d  ",si[i]);

        printf("\n");        */

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
         // printf("nuevo enlace %d con %d \n",m+s,g);

          }       //fin dle bucle sobre los m links nuevos

  /* for(i = 1; i <=N ; i++)
	   {
	    for(j = 1; j<=m ; j++)
           printf("M[%d][%d] = %d \n",i,j,M[i][j]);
	   }
         */

       if(A[m+s]==0)
         {
         A[m+s]=ValorA;

         for(j=m+s; j<=N; j++)
            {P[j]=P[j]+A[m+s];}
         }

      norma=P[N];


 }  

   //fin del bucle a los nodos de la red
		
	for(i=1; i<=N; i++)
          {k[i]=k[i]+ m;}

	/*
	for(i=1;i<=N;i++)
	  {
	    printf(" \n %i      ", i);
	    for(j=1;j<=m;j++)
	      {
		printf("%i  ", M[i][j]);
	      }
	  }

	//getch();
	printf("\n  \n");
	*/

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
             C[i][j]=M[i][j];  //en las m primeras posiciones de C, copio M

         for(j=1;j<=m;j++)
            {                            //a partir de la cuarta posicion, voy añadiendo
             z=M[i][j];
	         x[z]++;
             C[z][x[z]+m]=i;
            }
        }


 /* for(i=1;i<=N;i++)
     {
     printf(" %d      ",i);
      for(j=1;j<=N;j++)
         {
         printf(" %d  ",C[i][j]);
         }
      printf("\n");
     }

  */



 //comprobaciones de:

 for(i=1;i<=N;i++)
    {
    for(j=1;j<=N;j++)
       {
       if(C[i][j]==i)       //si esta conectado a si mismo
          {
          fprintf(escribe2,"error: conexion consigo mismo, iteracion numero %d  nodo %d\n",n, i );
          //printf("error de conexion tipo 1, iteracion numero %d  nodo %d\n",n, i );
          }
       }

    for(j=2;j<=k[i];j++)
       {
       for(y=1; y<=j-1; y++)
         {
          if(C[i][j-y]!=0  &&  C[i][j-y]==C[i][j])          // si existen enlaces dobles
             {
              fprintf(escribe2,"error: enlace doble, iteracion numero %d  nodo %d\n",n, i );
              //printf("error de conexion tipo 2, iteracion numero %d  nodo %d\n",n, i );
             }
         }
       }
    }




    //costruccion del histograma P(k)

    for(i = 0; i <= N; i++)           //inicializo
		PK[i] = 0;
	
	for(i =1; i <= N; i++)            //recuento
        PK[k[i]]++;

  	for(i = 1; i <= N; i++)         //normalizo
   		{
		  //PK[i] = PK[i]/N;
                  PK_tot[i]+=PK[i];
        }


  } //    fin del bucle en iter    (estadistica)



	// Normalizo PK_tot
for( i = 0; i < N; i++)
   {PK_tot[i] = PK_tot[i]/iter;}



	/* Escribo k y PK en un archivo */

for(i = 1; i<N; i++)
   { fprintf(escribe,"%d  %lf\n", i, PK_tot[i]);}


fclose(escribe);
fclose(escribe2);

getch();
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
		
