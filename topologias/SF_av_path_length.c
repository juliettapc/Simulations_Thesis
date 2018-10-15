
// creacion de una red SF 

# include <stdio.h>
# include <stdlib.h>
# include <math.h>


# define N 4500 //tamaño de la red
# define m_o  2      //nodos inicialmente unidos
# define m   2     //nodos nuevos añadidos a cada paso de tiempo



# define Niter 20        //estadistica

# define K_MAX   500          //para el tamaño de C[N][K_MAX]     ojo pq depende de N en SF!!!



//Generacion de numeros aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];



 
int k[N+1];                    //conectividad topológica 
int   M[N+1][m+1], unido[m+1], Cuidado[N+1];         //  unido[] guarda a quién se ha unido i , y cuidado[] quien le ha lanzado ya un link
double P[N+1], P_prov[N+1];        //P_prov[] es la provisional para manipularla
double PK[N+1], PK_tot[N+1];        //  para la distribucion P(k)
double r,v;                                  //para guardar los aleatorios
int  i,ii, j, jj,  w,g,gg,  q, x[N+1],y,z, C[N+1][K_MAX+1],n, ValorA, steps,s,flat;
int norma, norma2, norma_aleat,tipo[m_o+1];
double dnorma_aleat;
int si[N+1], cont,cont2;

//variables para el dilema

int  e[N+1],  w, t, k_m; 
double  ben[N+1],p;
double ro, ro_min, ro_max, delta_ro, c_media;      //media (sobre Niter realizaciones) de cooperadores (sin discriminar)
double cp_media,dp_media,fl_media;    //id pero sobre las Niter realizacines
double n_coop_equil;   //media de cooperadores en la ventana del equilibrio
int  n_coop; //numero de coop. instantaneos -ultimo paso de tiemp antes del equil-
int  cp,dp,fl;       //cooperadores puros, defectores puros y fluctuantes   (en una realizacion)
int S[N+1],Estado[N+1][2];      //etiqueta para discriminar entre puros y fluctuantes
int iter;
int marca;    //marcador para saber si ha salido n_coop=0
double b;
int K_max;

int actual[N+1],d[N+1],dist,h;
double dist_tot,D_med;

void inicia_rand(int semilla);


char file7[256];    
FILE *fich7;         




int main()
{
  
  
  
             //siempre dentro del main!!!
  
  
    sprintf(file7,"Pk_%d_SF.dat",N);        
  fich7=fopen(file7,"wt");
  fclose(fich7);
  
   inicia_rand(time(0));
 

   D_med=0.;

printf("\n\nN:%d\n",N);
  for(n=1;n<=Niter;n++)
    {
      
      
      for(i =1; i <= N; i++)
	PK[k[i]]=0;
      
      
      steps=N-m_o;       //numero de nodos que se añadiran en total
      
      
      

      
      
      for(i = 0; i<=N; i++)
	{
	  k[i] = 0;     
	  
	  for(j = 0; j <= m_o; j++)
	    {
	  M[i][j] =0;
	    }        
	}
      
      
      for(j = 1; j <=m_o ; j++)         // inicializo los m_o primeros nodos conectados entre si
	k[j]=m_o-1;
      
      
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
      
      
      for(j=1; j<=N; j++)
	{ 
	  P[j] = P[j-1] + k[j];
	  //printf("P[%d] = %lf \n",j,P[j]);
	}
      
      
      
      norma=P[m_o];     
      norma2=norma;
      
      
      //empieza la construccion de la red:
      
     
      for(s = 1; s <= steps; s++)    //bucle a todos los demas nodos (paso de tiempo)
	{
	  
	  //printf("\nnodo: %d \n",s+m);
	  
	  
	  for(i=0; i<=m; i++)
	    {unido[i]=0;}
	  
	  for(i=1; i<=N; i++)
	    {P_prov[i]=P[i];}        //guardo la prob en un vector aux para manipularlo
	  
	  norma2=norma;
	  
	  
	  
	  
	  for(q = 1; q<=m; q++)        //bucle a los m links nuevos que voy a añadir
	    {
	      
	      v=FRANDOM*norma2;
	
	      //getchar();
	      
	      for(j=1; j<=m_o+s; j++)
		{
		  if(v<P_prov[j])        //elijo el nodo
		    {
		      unido[q]=j;
		      //printf(" %f   %d   %f nodo:%d\n",v,norma2,FRANDOM, j);
		      break;
		    }
		}
	      
	      g=unido[q];
	      P_prov[g]=0;    //lo guardo y anulo su prob para no volver a cogerlo
	      
	      for(j=g+1; j<=m_o+s; j++)   
		{
		  P_prov[j] = P_prov[j] - k[g] ;
		}
	      
	      norma2=norma2-k[g];
	      
	      //printf(" norma:   %d \n",norma2);
	    }    //fin del bucle sobre los m links lanzados
	  
	  
	  
	  for(j=1;j<=m_o;j++)
	    {
	      g=unido[j];
	      k[g]=k[g]+1;     
	      M[m_o+s][j]=g;
	      
	      for(ii=g;ii<=m_o+s;ii++)
		{
		  P[ii]=P[ii]+1;    ////////antes la prob prov
		}
	      
	    }
	  
	  k[m+s]=m;
	  norma=norma+2.0*m;      /////////////antes +m_o
	  P[m+s]=norma;
	  
	  
	  
	}  
      
      //fin del bucle a los nodos de la red
      
      
      
      
      
      //obtencion de la matriz C a partir de la M
      
      
      for(i=0;i<=N;i++)
	{
	  x[i]=0;
	  
	  for(j=0;j<=K_MAX;j++)
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
      

   	//red	de prueba:

      /*	C[1][1]=2;
	C[2][1]=1;
	C[2][2]=4;
	C[3][1]=4;
	C[4][1]=2;
	C[4][2]=3;
	C[4][3]=5;
	C[4][4]=6;
	C[5][1]=4;
	C[6][1]=4;
	C[6][2]=7;
	C[6][3]=9;
	C[7][1]=6;
	C[8][1]=6;
	C[9][1]=6;
	C[9][2]=10;
	C[9][3]=11;      ////////  debe salir: <l>= 2.57
	C[10][1]=9;
	C[11][1]=9;


	C[1][1]=2;
	C[2][1]=1;
	C[2][2]=3;
	C[3][1]=2;
	C[3][2]=4;
	C[4][1]=2;
	C[4][2]=3;
	C[5][1]=4;
	C[5][2]=6;
	C[6][1]=5;
	C[6][2]=7;
	C[7][1]=6;
	C[7][2]=8;
	C[8][1]=7;
	C[8][2]=9;
	C[9][1]=8;*/
      

    ////////  debe salir: <l>= 3.33333


  



      //calculo de average path length:   

      
      dist_tot=0.;
      
      for (i=1;i<=N;i++)     // para todos los nodos de la red
	{
	  dist=0;    //guardara la suma de las distancias de un nodo a todos los demas

	  for (j=1;j<=N;j++)
	    {
	      d[j]=-100;   //inicializo las banderas de los nodos como "no visitado"
	    }	               //guardara la distancia desde el nodo i a todos los demas


	  cont=1;
	  actual[cont]=i;     //guardo el nodo que estoy mirando en este momento (empiezo por mi)
	  d[i]=0;   //para que si parto de mi, no vuelva a mi mismo de nuevo

	  for(j=1;j<=cont;j++)
	    {	   

	      for(h = 1; h <= k[actual[j]]; h++)     //bucle sobre la conectividad del nodo que estoy mirando
		{
		  if(d[C[actual[j]][h]]==-100)   //si no he visitado el nodo
		    {
		      cont++;
		      actual[cont] = C[actual[j]][h]; 
		      d[C[actual[j]][h]]=d[actual[j]]+1;   //lo marco con 1+ la distancia del nodo actual
		    }
		}
	    }
	  


	  for (j=1;j<=N;j++)
	    {
	      if (d[j]>0)
		{
		  dist+=d[j];
		  //printf("dist:%d\n",dist);
		}
	    }
	  
	  dist_tot+=dist;
	  //printf("dist_tot:%lf\n",dist_tot);
	  
	}// fin del bucle a todos los nodos de la red
      


      dist_tot=dist_tot/((N-1)*N);




printf("dist_tot:%lf\n",dist_tot);






      
      /* Calculo la distribucion de conectividad */
      
      
      for(i =1; i <= N; i++)
	PK[k[i]]++;
      
      
      
  
      /* Normalizo PK */
      for( i = 1; i < N; i++)
	{
	  PK[i] = PK[i]/N;
	  PK_tot[i]+=PK[i];
	}
      
      
      
      
      //calculo la maxima conectividad:
      K_max=0;
      for(i=1;i<=N;i++)
	{
	  if(k[i]>K_max)
	    K_max= k[i];
	}
      
      
      printf("conectividad max: %d en iter %d\n",K_max,n);
      
      
      D_med+=dist_tot;
      
    }    //fin bucle estadistica
  
  //escribo la p(k):
  

 for( i = 1; i < N; i++)
   {
     PK_tot[i] = PK_tot[i]/Niter;
   }
 
 D_med=D_med/Niter;
 printf("en %d iteraciones, d_med=%lf\n",Niter,D_med);
  
   fich7=fopen(file7,"at");
  for(i = 2; i <= N; i++)
    
    fprintf(fich7,"%d   %lf\n",  i, PK_tot[i]); 
  
  
  
    fclose(fich7);
  
  
}  //fin de la funcion "red()"






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





