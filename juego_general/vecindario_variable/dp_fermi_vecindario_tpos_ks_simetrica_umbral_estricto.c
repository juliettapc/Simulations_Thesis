// 16/9/09:     Programa para la implementacion del DP en redes estáticas pero 
//     jugando una ronda con k_star vecinos, de entre todos los que tenga un nodo, en cada ronda seran distintos
//    basado en el juego_general_k-new.c   pero con


//  matriz de payoff dependiente del cociente de  b (beneficio) y c (coste) :  b/c=ratio  es ahora de 
// la forma:
//                  (    ratio-1     -1   )     ( R   S )                                 
//             M=   (                     )  =  (       )            (declarare igualmente R, T, S y P )
//                  (      ratio      0   )     ( T   P )                    
//



// y con la prob de cambio de estrategia dada por la regla de Fermi


// OJO! CON LA PROB DE CAMBIO DE ESTRAT DE FERMI, EL SISTEMA ACABARA O BIEN EN EL ESTADO C=0 O C=1, NO EN UNO INTERMIEDIO MEDIO!!!! por tanto, no tiene sentido medir cp, dp y fl


//ademas, en cada ronda, cada nodo elige de entre sus vecinos topologicos, solo a k_star de ellos, con los que jugara
//  y con uno de los cuales tb luego se comparara

// 7-10-09:  añado calculo de los tiempos de fijacion al estado absorbente (c=0 o c=1) del sistema (por ahora sin discriminar entre ambos estados)

//  8-10-09: estudio de la kefectiva (kin + kout) vs ktopologica

// 24-10-09: modificacion para que las interacciones sean simetricas i con j -->> y j con i

// 16-12-09: añado umbral para que si coop>0.95N o <(1-0.95)N, acabe la simu
// para evitar los problemas de lentitud de convergencia al estado absorbente
//OJO!! PERO  CON ESTE PROGRAMA NO CALCULARE POR TANTO TIEPOS DE FIJACION!!

// 13-1-10: como el umbral del 5% o 95% aun sigue haciendo simus muy lentas
// dejo evolucionar el sistema 5000 pasos temporales y miro si coop < o > del 50
// y le asigno un valor final de 0 o 1, respectivamente!!
//   ( a efectos de simulacion, ahora el umbral=0.5)

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>


# define Niter   50    //estadistica


# define N    4000       //tamaño de la red
# define m_o  2      //nodos inicialmente unidos
# define m   2      //nodos nuevos añadidos a cada paso de tiempo
# define  ro 0.5 
# define alfa 0.0


# define ratio_min    6.51    //cociente   beneficio/coste  (siempre >1)
# define ratio_max    10.01
# define delta_ratio  0.1


# define seleccion    1.0




 
# define KSTAR     10  // numero de vecinos con los que juego/me comparo (entre m y N)



# define jugadas 500000  //pasos de tiempo para el juego (transitorio)



# define umbral   0.5    // tanto por uno de cooperadores que ya considero coop total (o, 1-umbral, defecc tot)


# define jugadas_umbral   2000   //momento en el cual evaluo si la coop es < o > del 50%

/////////////////////////Generacion de numeros aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
////////////////////////////////////////////////////////



char nombre3[256],nombre2[256];

 
int k[N+1],k_PA[N+1],A[N+1];                    //conectividad topológica y de Preferential Attachment (PA)
int   M[N+1][m+1], unido[m+1], Cuidado[N+1];         //  unido[] guarda a quién se ha unido i , y cuidado[] quien le ha lanzado ya un link
double P[N+1], P_prov[N+1];        //P_prov[] es la provisional para manipularla
double PK[N+1], PK_tot[N+1];        //  para la distribucion P(k)
double r,v;                                  //para guardar los aleatorios
int  i, j, jj,  w,g, d, q, x[N+1],y,z, C[N+1][N+1],n, ValorA, steps,s,flat;
int norma, norma2, norma_aleat,tipo[m_o+1];
double dnorma_aleat;
int si[N+1], cont,cont2;
int k_star[N+1],C_star[N+1][N+1],C_aux[N+1][N+1];

//variables para el dilema

int  e[N+1], e_aux[N+1],  w, t, k_m; 
double  ben[N+1],p;
double  c_media;      //media (sobre Niter realizaciones) de cooperadores (sin discriminar)
int  n_coop; //numero de coop. instantaneos -ultimo paso de tiemp antes del equil-
int iter;
int marca;    //marcador para saber si ha salido n_coop=0
double ratio;

            //para estudiar los clusters de coop y los tiempos



double R, Pu,Su,T;      // los coef de la matriz de pagos tendran valores como funcion del ratio

int k_efect[N+1];

	

void inicia_rand(int semilla);
void construir_red();
void juego();
void selecciono_vecindario();


FILE *escribe3,*escribe2;






int main()
{
  
    //  inicia_rand(26);

     inicia_rand(time(0));
    
    
    printf("\n\n\nPrograma PD con vecindad variable   --umbral--\n\n"); 
    printf("Red de N=%d    alfa=%lf  (%d iteraciones) seleccion=%.2lf  K_star:%d\n",N,alfa,Niter,seleccion,KSTAR);
    
    
    
    //Guarda <c>, CP,DP y  F : un solo fichero (pq ahora no barrermos en rho)
    sprintf(nombre3,"v__C_med_ratio%.2lf-%.2lf_alfa%.1lf_selecc%.2lf_kstar%d_N%d_%diter_sim_umbr.dat",ratio_min,ratio_max,alfa,seleccion,KSTAR,N,Niter); 
    escribe3=fopen(nombre3,"wt");  
    fclose(escribe3);
    
    
    
    
    
    
    
    
    ratio=ratio_min;
    
    while(ratio<=ratio_max)    //bucle externo para barrer en b
    {
	printf("\n\nb/c=%lf\n",ratio);  
	
	

	R=ratio-1.0;        //doy valores a los parametros de la matriz de pagos en funcion de b/c
	Su=-1.0;
	T=ratio;
	Pu=0.0;
	



	
	sprintf(nombre2,"v__Evoluc_coop_ratio%.2lf_alfa%.1lf_selecc%.2lf_kstar%d_N%d_%diter_umbr.dat",ratio,alfa,seleccion,KSTAR,N,Niter); 
	escribe2=fopen(nombre2,"wt");  
	fclose(escribe2);



	  
	
	


	c_media=0.0;



	for(iter=1;iter<=Niter;iter++)
	{
	    
	     printf("\niteracion:%d\n",iter);
	    
	     if(iter<=5)
	       {
		 escribe2=fopen(nombre2,"at");
		 fprintf(escribe2,"\n");    
		 fclose(escribe2);
	       }

	    
	    construir_red();

	    //printf("red construida\n");  

	    juego();    //incluye las "jugadas"  

	    // printf("juego hecho\n");  
	    
	    c_media+=n_coop;

	    //  printf("n_coop=%d\n",n_coop);

	 
	 	    
	}       //fin bucle estadística
	
	

	 c_media=c_media/(double)(Niter*N);            //normalizar, y tanto por uno	
	
	




	printf("\nb/c=%lf c_media:%f  \n",ratio ,c_media); 
	




	escribe3=fopen(nombre3,"at");
	fprintf(escribe3,"%f  %f\n",ratio,c_media);      
	fclose(escribe3);
	
	
	
	
  	ratio+=delta_ratio;    
	
	
    }     //fin del bucle en b/c
    
    
    
    
    

}



/////////////////////////////
////////////////////////////
//////////////////////


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




/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////




void juego ()    //ojo!!!! aqui el convenio es 0=defect y 1=coop!!!!!!!!
{
    int i,j,tpo,y,w;
    double r,p;



    for(i=1;i<=N;i++)
    {
	ben[i]=0;
	e[i]=0;  //todo defectores
    }


    n_coop=0;
    for(i=1;i<=N;i++)       //c.i.
    {
	r=FRANDOM;
	if(r<ro)
	{
	    e[i]=1;
	    n_coop++;
	}
    }
    
    if(iter<=5)
      {	
	escribe2=fopen(nombre2,"at");
	fprintf(escribe2,"0  %d\n",n_coop);      
	fclose(escribe2);
      }



    tpo=0;
    while(tpo<jugadas)
    {

//printf("tpo:%d\n",tpo);  
	selecciono_vecindario();
	//printf("   seleccinado vecindario\n");  

	for(i=1;i<=N;i++)   //pongo a cero los beneficios de todos
	{
	    ben[i]=0;
	}

	for(i=1;i<=N;i++)   //cada nodo juega con todos sus vecinos
	{
	    for(j=1;j<=k_star[i];j++)
	    {
		y=C_star[i][j];
		if(y>i)  //para no repetir parejas de nodos
		{
		    if(e[i]==1 && e[y]==1)   //ambos coop
		    {
			ben[i]+=R;
			ben[y]+=R;
		    }  
		    
		    if(e[i]==1 && e[y]==0)   // i coop, y defect
		    {
			ben[i]+=Su;
			ben[y]+=T;
		    }  
		    
		    if(e[i]==0 && e[y]==1)   //i defect, y coop
		    {
			ben[i]+=T;
			ben[y]+=Su;
		    }  
		    
		    
		    if(e[i]==0 && e[y]==0)   //ambos defect
		    {
			ben[i]+=Pu;
			ben[y]+=Pu;
		    }  


		}

	    }    //i acaba de jugar con sus vecinos
	    
	    	    
	}  // fin de las partidas de todos

	//printf("   fin de la partida\n");  


	for(i=1;i<=N;i++)   //copio la estrategia de antes de la comparcion
	{
	    e_aux[i]=e[i];  
	}


	for(i=1;i<=N;i++)   //comparan beneficios
	{
	    
	    r=FRANDOM; //elijo un vecino al azar
	    r=r*k_star[i];
	    w=(int)r+1;
	    y=C_star[i][w];


	    p=seleccion*(ben[i]-ben[y]);  //prob de que i imite a  y
	    p=exp(p);
	    p=1.0/(1.0+p);



	    r=FRANDOM;
	    if(r<p)
	    {
		e_aux[i]=e[y];
	    }



	}// fin de las comparaciones


	for(i=1;i<=N;i++)   //actualizaciones en paralelo
	{
	    e[i]=e_aux[i];
	}


	n_coop=0;
	for(i=1;i<=N;i++)
	{
	    if(e[i]==1)	    
		n_coop++;	      
	}

	
	if(iter<=5)
	  {
	    escribe2=fopen(nombre2,"at");
	    fprintf(escribe2,"%d  %d\n",tpo,n_coop);      
	    fclose(escribe2);
	  }
	
	//	printf("jugada:%d   n_coop:%d\n",tpo,n_coop);

        if(n_coop==0)     //me salto el resto de las jugadas!!!
        {
            break;
        }

        if(n_coop==N)
        {
            break;
        }


	
	if(tpo>=jugadas_umbral)
	  {
	    
	    if(n_coop >= (umbral*N))
	      {
		n_coop=N;
		break;  
	      }
	    else
	      {
		n_coop=0;
		break;  
	      }
	  }
	
	

/*if( (tpo>=jugadas_umbral) && (n_coop>umbral*N) )
	{
	    n_coop=N;
	    break;  
	}
	else
	{
	    if( (tpo>=jugadas_umbral) && ( n_coop<(1.0-umbral)*N))
	    {
		n_coop=0;
		break;  	
	    }
	    }*/




	tpo++;

    }  //fin bucle while (jugadas)

    
    /*if(tpo==jugadas &&  n_coop!=0  &&  n_coop!=N)   //si se acaban las jugadas y no he alacanzado el estado absorbente
    {
	printf("\n-------no he alacanzado el estado absorbente");
	tpo_f[iter]=jugadas;

	media_tpo_fijac+=jugadas;
	
	}*/
  
    
    
}


////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////





void selecciono_vecindario()
{

    int i, j, w,cont,vecino;   // el tamaño de C_aux[][] da problemas!!!!!!!!!!!!!!!! pero no si lo declaro en el main
    double  r;

	

    for(i=1;i<=N;i++)
    {
	for(j=1;j<=N;j++)
	{	   
	    C_aux[i][j]=C[i][j]; //copio la matriz para hacerle perrerias
	
	    C_star[i][j]=0;



	}
	k_star[i]=0;

	k_efect[i]=0;   //kin +kout
    }


   
    for(i=1;i<=N;i++)        //recorro toda la red   AHORA LAS INTERACCIONES SON SIMETRICAS!! 
    {

//	printf("nodo:%d\n",i);

	if(k[i]<KSTAR)
	{
	    k_star[i]=k[i];

	    for(j=1;j<=k[i];j++)
	    {

		k_efect[i]++;
		C_star[i][k_efect[i]]=C[i][j];
			

		vecino=C[i][j];
		k_efect[vecino]++;
		C_star[vecino][k_efect[vecino]]=i;
	    }
	}
	else
	{
	    k_star[i]=KSTAR;

	    cont=1;
	    while(cont<=KSTAR)
	    {
		r=FRANDOM;             //elijo un vecino al azar
		r=r*k[i];
		w=(int)r+1;

		if(C_aux[i][w]!=0)     //si i no ha elegido a w ya antes
     		{

		    k_efect[i]++;
		    C_star[i][k_efect[i]]=C[i][w];		   
		    C_aux[i][w]=0;

		    vecino=C[i][w];
		    k_efect[vecino]++;

		    C_star[vecino][k_efect[vecino]]=i;


		    cont++;
		}

		
	    }
	}


    }      //fin bucle nodos de la red





    
}
