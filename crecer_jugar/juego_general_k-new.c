 //Implementacion de un juego con matriz payoff general en
//una red gardeñes-moreno
//obtencion de la dependencia de la concentracion de cooperantes puros,
//defectores puros y fluct.en funcion de la concentracion inicial de cooperantes
//  esta version incluye un bucle externo para barrer en b, ademas del de barrer en rho

// calculo tb  de la densidad de coop puros, defect puros y fluct en funcion de su  conectividad


//CORRECCION en la normalización: diferencia máxima de coef de la matriz de payoff*max k (antes: b*k_m siempre)


// estudio de las fluctuaciones (desviacion estandar) tanto de C como de los beneficio (suma sobre N de los payoff), para 
//   ver si existe una divergencia en torno a b=1.4  (encontrada en el modelo dipolar)

# include <stdio.h>
# include <stdlib.h>
# include <math.h>


# define N 1000 //tamaño de la red
# define m_o  3      //nodos inicialmente unidos
# define m  3      //nodos nuevos añadidos a cada paso de tiempo
 
# define alfa 0.0



#define nb 28    //numero de valores de b


#define R 1.0
#define Pu 0.0


#define Su  0.0


# define jugadas 25000  //pasos de tiempo para el juego (transitorio)
# define jugadas_equil 30000    //pasos de tiempo para hacer la media y discriminar
# define Niter  1000      //estadistica

# define K_MAX   500          //para el tamaño de C[N][K_MAX]     ojo pq depende de N en SF!!!

//Generacion de numeros aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];

char nombre[256],nombre3[256],file7[256],file8[256],file9[256], nombre4[256];

 
int k[N+1],k_PA[N+1],A[N+1];                    //conectividad topológica y de Preferential Attachment (PA)
int   M[N+1][m+1], unido[m+1], Cuidado[N+1];         //  unido[] guarda a quién se ha unido i , y cuidado[] quien le ha lanzado ya un link
double P[N+1], P_prov[N+1];        //P_prov[] es la provisional para manipularla
double PK[N+1], PK_tot[N+1];        //  para la distribucion P(k)
double r,v;                                  //para guardar los aleatorios
int  i, j, jj,  w,g,gg, d, q, x[N+1],y,z, C[N+1][K_MAX+1],n, ValorA, steps,s,flat;
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

            //para estudiar los clusters de coop y los tiempos
double NCLUSTERSMED,NCLUSTERSMEDd,GCMED, GCMEDd,GCompd,GComp, CLUSTERS,CLUSTERSd,da,Tiempos[jugadas_equil+1],VentanasCoop[jugadas_equil+1];

int NN[N+1],Phase[N+1][K_MAX+1],warning[N+1],analyse[N+1],semilla[N+1],TC[N+1],numvar[N+1];
int NonNull,NonNullCP,NonNullDP,NCLUSTERS,NCLUSTERSd;

double normatiempos,normaVentanasCoop,FREC[N+1], FRECMED,frecmed[K_MAX+1] , max;

int Validas,Salir, contador,size,tiempo, ib;

double frecfinal,frecPuros,frecDPuros,dobleprec;

double NORMAFREC,normafrec[K_MAX];

int  kk;   

double  rho_C_k[K_MAX+1], rho_D_k[K_MAX+1], rho_Fl_k[K_MAX+1] ,rho_C_tot[K_MAX+1], rho_D_tot[K_MAX+1], rho_Fl_tot[K_MAX+1],iter_k[K_MAX+1];     

double  media_ben, sigmaC,sigmaB, media_sigmaB,media_sigmaC;
double cooperacion [jugadas_equil+1], suma_ben[jugadas_equil+1];
	

void inicia_rand(int semilla);
void construir_red();
void histograma();
void comprobar_red();
void juego();
void juego_equil();
void TopologiaCP();
void TopologiaDP();
void histograma_rho_k();

FILE *escribe;
FILE *escribe2;
FILE *escribe3;
FILE *fich7,*fich8,*fich9,*fich10,*fich11, *escribe4;



double b=1.0,  delta_b=0.1;

int main()
{
  
  
  
  printf("\n\n\nPrograma juego general (corregida normalizacion prob. cambio estrat.)\n\n"); 
  printf("Red de N=%d   alfa=%lf  (%d iteraciones)\n",N,alfa,Niter);
  printf("S=%lf P=%lf \n", Su, Pu);  


  //Guarda los Clusters de Cooperantes etc.
  sprintf(file7,"Clusters_b%.2lf-%.2lf_alfa%.1lf_S%.2lf_P%.2lf_N%d_%diter_c.dat",b,b+nb*delta_b,alfa,Su,Pu,N,Niter);
  fich7=fopen(file7,"wt");
  fclose(fich7);
  
  //Guarda <c>, CP,DP y  F : un solo fichero (pq ahora no barrermos en rho)
  sprintf(nombre3,"fluct_b%.2lf-%.2lf_alfa%.1lf_S%.2lf_P%.2lf_N%d_%diter_c.dat",b,b+(nb-1)*delta_b,alfa,Su,Pu,N,Niter);    
  escribe3=fopen(nombre3,"wt");  
  fclose(escribe3);

 
 //Guarda la densidad de coop puros, defect puros y fluct en funcion de la conectividad : un solo fichero 
  sprintf(nombre4,"rho_vs_k_b%.2lf-%.2lf_alfa%.1lf_S%.2lf_P%.2lf_N%d_%diter_c.dat",b,b+(nb-1)*delta_b,alfa,Su,Pu,N,Niter);    
  escribe4=fopen(nombre4,"wt");  
  fclose(escribe4);




  printf("b_ini=%lf nb=%d delta_b=%lf\n",b,nb, delta_b);
  
  for(ib=0;ib<nb;ib++)     //bucle externo para barrer en b
    {
      
      ro_min=0.5;
      ro_max=0.51;
      delta_ro=0.1;
      
      printf("ro_min=%lf hasta: ro_max=%lf \n\n",ro_min,ro_max);
      printf("b=%lf\n\n",b);
      
      /*sprintf(nombre,"pk_N%d_alfa%.3lf.dat",N,alfa);
	escribe=fopen(nombre,"wt");
	escribe2=fopen("comprobar.dat","wt");*/
      
      
      
      
      
      inicia_rand(time(0));
      /* comprobar=time(0);
	 printf("time0= %f \n",comprobar);*/
      
      
      ro=ro_min;
      while(ro<=ro_max)     //bucle para barrer en rho
	{
	  
	  //Guarda la Distribucion de Tiempos de Cooperacion : un histogr. para cada valor de rho y de b
	  sprintf(file8,"Histogr_Times_b%.2lf_alfa%.1lf_S%.2lf_P%.2lf_N%d_%diter_c.dat",b,alfa,Su,Pu,N,Niter);
	  fich8=fopen(file8,"wt");
	  fclose(fich8);	
	  
          //Guarda la Distribucion de Frecuencias versus K: un histogr. para cada valor de rho y de b
	  sprintf(file9,"Frecs_k_b%.2lf_alfa%.1lf_S%.2lf_P%.2lf_N%d_%diter_c.dat",b,alfa,Su,Pu,N,Niter);
	  fich9=fopen(file9,"wt");
	  fclose(fich9);
	  
	  
	  //inicializo variables glabales (mediaran sobre estadística)	    
	  
	  for(i=0;i<=jugadas_equil;i++)
	    {
	      Tiempos[i]=0.; //Guarda la distribucion de tiempos de Coop. tot de los fluct
	      VentanasCoop[i]=0.;   //histogr ventanas tiempo coop de los fluct
	    }
	  
	  normatiempos=0.; //Norma para la distribucion de tiempos de Coop.
	  normaVentanasCoop=0.;
	  FRECMED=0.; //Frecuencia media de cambio.
	  NORMAFREC=0.; //Numero de tios que han cambiado.
	  
	  for(i=0;i<K_MAX;i++)
	    {
	      frecmed[i]=0.;  //Distrib. de Frecs. en funcion de K.
	      normafrec[i]=0.;  //Numero de oscilantes de conectividad K
	    }
	  
	  
	  NonNullCP=0;
	  NonNullDP=0;
	  NCLUSTERSMED=0.;
	  NCLUSTERSMEDd=0.;
	  GCMED=0.;
	  GCMEDd=0.;
	  Validas=0;

	  
	  for(i=1;i<=K_MAX;i++)    //para los histogr de dens. cp dp y fl en funcion de k
	    {
	      rho_C_tot[i]=0.;
	      rho_D_tot[i]=0.;
	      rho_Fl_tot[i]=0.;
	      iter_k[i]=0.;
	    }
	  

	  
	  
	  c_media=cp_media=dp_media=fl_media=0;
	  media_ben=media_sigmaC=media_sigmaB=0.0;
	 
	  for(iter=1;iter<=Niter;iter++)
	    {
	      
	      printf("iteracion:%d\n",iter);
	      
	      for(i=1;i<=N;i++)      //inicializo variables locales (mediaran en una iteracion)
		{
		  TC[i]=0;         //Guarda el Tiempo de Coop del nodo i.
		  FREC[i]=0.;      //Guarda la Frec media (=tiempo medio de coop) del nodo i.
		  numvar[i]=0;     //Numero de veces que cambia i.
		  Estado[i][0]=1;  //los inicializo a defectores
		  Estado[i][1]=0;	
		}
	      
	      
	      construir_red();

	      //comprobar_red();

	      histograma();    //lo calcula pero no lo imprime ni lo normaliza  a 1 (solo recuento)
	      
	      juego();     //  transitorio
	      
	      
	      //guardo la estrategia de la ultima jugada antes del equilibrio,	
	      // para luego  usarla para discriminar   ( S[i] = Marcador[i] ):
	      
	      for(i=1;i<=N;i++)    //mi convenio: 0=CP, 1=DP, -1=Fl.
		{                  		    
		  //		S[i]=-1;   //todos como fluct  	(AUN NO HAY FLUCTUANTES)
		  
		  if(e[i]==0)      //si es coop
		    {
		      Estado[i][0]=0;  
		      Estado[i][1]=0;
		      S[i]=0;        //lo marco como coop  puro
		    }  
		  else
		    {
		      S[i]=1;    //y si es defect.,  como defect. puro
		      //Estado[i][0]=1; 
		    }         
		  
		} 
	      
	      if(marca==1)
		{
		  juego_equil();    //  jugadas_equil 
		  Validas++;        //recuento numero de realizaciones no nulas
		}
	      
	      //printf("n_coop=%d\n",n_coop);
	      c_media+=n_coop_equil;
	      cp_media+=cp;
	      dp_media+=dp;
	      fl_media+=fl;
	      
	      
	      //   ( calculos de histogramas tiempos al final de la funcion juego_equil )
	      
	      
	      TopologiaCP();
	      TopologiaDP();
	      

	      histograma_rho_k();      

	      for(i=1;i<=K_MAX;i++)     //acumulo
		{
		  rho_C_tot[i]+=rho_C_k[i];
		  rho_D_tot[i]+=rho_D_k[i];
		  rho_Fl_tot[i]+=rho_Fl_k[i];
		}
	      
	      
	      
	      
//calculo la desviacion de C y de los beneficios en una iteracion

	      sigmaC=0.0;
	      sigmaB=0.0;	      

	      for(i=1;i<=jugadas_equil;i++)
	      {
		  sigmaC+=(cooperacion[i]-n_coop_equil)*(cooperacion[i]-n_coop_equil);
		  sigmaB+=(suma_ben[i]-media_ben)*(suma_ben[i]-media_ben);
	      }
	      
	      sigmaC=sigmaC/(jugadas_equil-1.0);
	      sigmaB=sigmaB/(jugadas_equil-1.0);

	      sigmaC=sqrt(sigmaC);
	      sigmaB=sqrt(sigmaB);


	      media_sigmaC+=sigmaC;     //para hacer la media sobre iteraciones
	      media_sigmaB+=sigmaB;


	    }       //fin bucle estadística
	  
	  
	  
	  c_media=c_media/(double)(Niter*N);            //normalizar, y tanto por uno
	  
	  cp_media=cp_media/(double)(Niter*N); 
	  dp_media=dp_media/(double)(Niter*N); 
	  fl_media=fl_media/(double)(Niter*N);


	  media_sigmaC=media_sigmaC/(double)Niter;     //para hacer la media sobre iteraciones
	  media_sigmaB=media_sigmaB/(double)Niter;
	  


	  
	  FRECMED=FRECMED/NORMAFREC;
	  
	  for(i=1;i<=K_MAX;i++)
	    {
	      if (normafrec[i]!=0)
		frecmed[i]=frecmed[i]/normafrec[i];
	      else 
		normafrec[i]=0;
	      
	    }
	  
	  
	  for(i=1;i<=K_MAX;i++)      //normalizo al numero de iteraciones
	{
	  rho_C_tot[i]= rho_C_tot[i]/iter_k[i];
	  rho_D_tot[i]= rho_D_tot[i]/iter_k[i];
	  rho_Fl_tot[i]=rho_Fl_tot[i]/iter_k[i];
	}
      
 






	 
      
      printf("\nb=%lf ro= %f  c_media:%f  cp_media:%f  dp_media:%f  fl_media:%f\n", b,ro,c_media,cp_media,dp_media,fl_media); 
      
      escribe3=fopen(nombre3,"at");
      fprintf(escribe3,"%f   %f   %f  %f  %f  %f  %f  %f  %f\n", Su,b, ro,c_media,cp_media,dp_media,fl_media,media_sigmaC,media_sigmaB);
      
      fclose(escribe3);
      
      //escribo ficheros tiempos etc
      
      
      
      fich9=fopen(file9,"at");
      
      for(i=1;i<=K_MAX;i++)          //histograma tiempos coop. en funcion de la conectividad
	{
	  fprintf(fich9,"%i   %lf    %lf   %lf   %lf   %i\n",
		  i,frecmed[i],Su,b, ro,Validas);
	  //printf("frecmed[%d]= %f  \n",i, frecmed[i]);
	}
      fprintf(fich9,"\n");
      fclose(fich9);
      
      
      
      fich8=fopen(file8,"at");
      //getchar();
      for(i=1;i<=jugadas_equil;i++)              
	{
	  if(normatiempos!=0.)
	    Tiempos[i]=Tiempos[i]/normatiempos;   //histograma tiempos coop. tot
	  else 
	    Tiempos[i]=0;
	  
	  if(normaVentanasCoop!=0.)
	    VentanasCoop[i]=VentanasCoop[i]/normaVentanasCoop;   //histograma ventanas de coop.
	  else	   
	    VentanasCoop[i]=0;
	  
	  fprintf(fich8,"%i   %lf   %lf   %lf   %lf   %lf\n",i,Tiempos[i],VentanasCoop[i],Su,b,ro);
	  //printf("Tiempos[%d]= %f    VentanasCoop[i]=%lf\n",i, Tiempos[i],VentanasCoop[i]);
	}
      
      
      //getchar();
      
      
      fprintf(fich8,"\n");
      fclose(fich8);	      
      
      
      
      //promedio clusters:
      if (NonNullCP>0)
	{
	  NCLUSTERSMED=NCLUSTERSMED/NonNullCP;
	  GCMED=GCMED/NonNullCP;
	}
      
      else 
	{
	  NCLUSTERSMED=0.;
	  GCMED=0.;
	}
      
      
      
      if (NonNullDP>0)
	{
	  NCLUSTERSMEDd=NCLUSTERSMEDd/NonNullDP;
	  GCMEDd=GCMEDd/NonNullDP;
	}
      else 
	{
	  NCLUSTERSMEDd=0.;
	  GCMEDd=0.;
	}	    
      
      
      fich7=fopen(file7,"at");
      fprintf(fich7,"%lf   %lf   %lf   %lf   %lf   %lf   %i   %i   %i   %lf\n",
	      Su,b,GCMED,GCMEDd,NCLUSTERSMED,NCLUSTERSMEDd,Validas,
	      NonNullCP,NonNullDP,ro);
      
      fclose(fich7);
      
      
      printf("GCMED=%.2lf   GCMEDd=%.2lf   NCLUSTERSMED=%.2lf   NCLUSTERSMEDd=%.2lf   NonNullCP=%i   NonNullDP=%i\n",GCMED,GCMEDd,NCLUSTERSMED,NCLUSTERSMEDd,NonNullCP,NonNullDP); 
      
      
      
      escribe4=fopen(nombre4,"at");
      for (i=1;i<=K_MAX;i++)
	{
	  fprintf(escribe4,"%lf   %lf   %d   %lf   %lf   %lf  \n",
		             Su,   b,   i,rho_C_k[i], rho_D_k[i], rho_Fl_k[i]);
	}
      fprintf(escribe4,"\n");
      fclose(escribe4);	
      
      
      
      ro=ro+delta_ro;
    }        //fin del bucle en rho
  
  
  
  
  b+=delta_b;    
  //fclose(escribe);
  //fclose(escribe2);
  
}     //fin del bucle en b


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


    /*   K_MAX=0;
    for(i=1;i<=N;i++)
    {
	if(k[i]>K_MAX)
	    K_MAX= k[i];
    }


    printf("conectividad max: %d\n",K_MAX);*/



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

}  //fin de la funcion "red()"



void comprobar_red()
{

 for(i=1;i<=N;i++)
    {
    for(j=1;j<=K_MAX;j++)                //MODIFICACION!!!  (antes: =0; <=N)  
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
     PK[i]=0.;
	
 for(i =1; i <= N; i++)            //recuento
     PK[k[i]]++;



 //
 //ahora no normalizo a 1, pq me interesa el Nº de nodos con cada conectividad, para normalizar las densid.
 //


 /*for(i = 1; i <= N; i++)         //normalizo
	{
	 PK[i] = PK[i]/N;
     PK_tot[i]+=PK[i];
     }

 */
	/* Escribo k y PK en un archivo 

for(i = 1; i<N; i++)
  {fprintf(escribe,"%d  %lf\n", i, PK[i]);}*/




}        





void juego()            //implementacion del dilema del prisionero:  periodo transitorio
{

//printf("the game begins \n");

n_coop = p = k_m = 0;


for(i=1;i<=N;i++)
   {
    ben[i]=0;
    e[i]=1;         //todo defectores
   }



for(i=1;i<=N;i++)           // establezco el % de coop. iniciales 
   {
   r=FRANDOM;
   //printf("r=%f \n",r);

   if(r<ro)
     {
     e[i]=0;
     n_coop++;
     }
   }


//printf("num coop.: %d  \n",  n_coop);

    marca=1;
for(n=1;n<=jugadas;n++)            //pasos de tiempo
   {
    n_coop=0;

    for(i=1;i<=N;i++)         
    {
	ben[i]=0;
    }
    
    for(i=1;i<=N;i++)          // cada nodo juega 1partida con cada vecino
    {
	for(j=1;j<=k[i];j++)
         {
          y=C[i][j];            
          if(y>i)            //para no repetir vecinos en una partida
          {  
          
	     if(e[i]==0 && e[y]==0)
	     {
		 ben[i]+=R;
		 ben[y]+=R;
	     
		 //printf("%d : %d  empate coop \n", i,y);
	     }

	     if(e[i]==0 && e[y]==1)
	     {
		 ben[i]+=Su;
		 ben[y]+=b;
		 //printf("%d : %d  coop contra defect\n", i,y);
	     }

	     if(e[i]==1 && e[y]==1)
	     {
		 ben[i]+=Pu;
		 ben[y]+=Pu;

		 //printf("%d : %d  empate defect \n", i,y);
	     }

	     if(e[i]==1 && e[y]==0)
	     {
		 ben[i]+=b;
		 ben[y]+=Su;
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
      
       if(k[i]>k[t])
       {k_m=k[i];}
       else
       {k_m=k[t];}

       //printf("juego1 ben[i]=%lf ben[t]=%lf\n",ben[i],ben[t]);//CONTROL


       if(ben[i]<ben[t])
       {
	   p=(ben[t]-ben[i]);

	   //CORRECCION en la normalización: diferencia máxima de coef de la matriz de payoff*max k (antes: b*k_m siempre)

	   if(Su>=0)
	   {
	       p=p/(b*(double)k_m);
	   }
	   else
	   {
	       p=p/((b-Su)*(double)k_m);
	   }
	   
	   v=FRANDOM;
	   if(v<p)
	   {
	      e[i]=e[t];
	    
	  }
       }

      }

    for(i=1;i<=N;i++)     //recuento de cooperantes
      {
      if(e[i]==0)
        n_coop++;
      }
    //printf("num coop.: %d  \n",  n_coop);


    if (n_coop==0) 
    {
	marca=0;
	cp=0.0;
	dp=N;
	fl=0.0;
	n_coop_equil=0;

	//printf("marca=0 \n");
	break;      //se saltara los pasos de tiempo que falten, pq n_coop no va a aumentar
    }
    
  
   }    //fin de los pasos de tiempo "jugadas"
                 
 
}


void juego_equil()     // jugar en el equilibrio: parto de la configuracion anterior
{                               //Mi convenio: -1 fluct, 0 coop, +1 defect
    
    
    
    cp=dp=fl=n_coop_equil=media_ben =0.0;    
    for(n=1;n<=jugadas_equil;n++)            //pasos de tiempo
    {
	n_coop=0.0;
	
	for(i=1;i<=N;i++)         
	{
	    ben[i]=0;
	}
	
	for(i=1;i<=N;i++)          // cada nodo juega 1partida con sus vecinos
	{
	    for(j=1;j<=k[i];j++)
	    {
		y=C[i][j];            //para no repetir
		if(y>i)
		{  
		    
		    if(e[i]==0 && e[y]==0)
		    {
			ben[i]+=R;
			ben[y]+=R;
			
			//printf("%d : %d  empate coop \n", i,y);
		    }
		    
		    if(e[i]==0 && e[y]==1)
		    {
			ben[i]+=Su;
			ben[y]+=b;
			
			//printf("%d : %d  coop contra defect\n", i,y);
		    }
		    
		    if(e[i]==1 && e[y]==1)
		    {
			ben[i]+=Pu;
			ben[y]+=Pu;
			
			//printf("%d : %d  empate defect \n", i,y);
		    }
		    
		    if(e[i]==1 && e[y]==0)
		    {
			ben[i]+=b;
			ben[y]+=Su;
			
			//printf("%d : %d  empate defect \n", i,y);
		    }
		}
		
	    }
	}                  //fin de la partida

	
	for(i=1;i<=N;i++)        //comparo estrategias
	{
	    r=FRANDOM;
	    r=r*k[i];
	    w=(int)r+1;
	    
	      
	    t=C[i][w];    //el vecino elegido al azar
	    
	    if(k[i]>k[t])
	    {k_m=k[i];}
	    else
	    {k_m=k[t];}


	    //printf("juego2 ben[i]=%lf\n",ben[i]);//CONTROL
	    
	    if(ben[i]<ben[t])
	    {
		p=(ben[t]-ben[i]);
		
//CORRECCION en la normalización: diferencia máxima de coef de la matriz de payoff*max k (antes: b*k_m siempre)
		
		if(Su>=0)
		{
		    p=p/(b*(double)k_m);
		}
		else
		{
		    p=p/((b-Su)*(double)k_m);
		}
		
		v=FRANDOM;
		if(v<p)
		{
		    e[i]=e[t];
		    
		}
	    }
	    
	}                   //realizados cambios de estrategia
	
	for(i=1;i<=N;i++)     //recuento de cooperantes instantaneos (medios, sin discriminar)
	{
	    if(e[i]==0)	    
		n_coop++;
	    
	}
	//printf("num coop.: %d  \n",  n_coop);
	
	n_coop_equil+=n_coop;
	


	cooperacion[n]=n_coop;           /////guardo C(t)
	cooperacion[n]=cooperacion[n]/N;

	for (i=1;i<=N;i++)      
	{
	    suma_ben[n]+=ben[i];        //guardo la suma instantanea de los beneficios
	}


	media_ben+=suma_ben[n];     // voy acumulando los beneficios para hacer la media temporal
 //discriminacion:

	for(i=1;i<=N;i++)         
	{
	    
	    if(e[i] != Estado[i][0])    //si ha cambiado en la ultima jugada
	    {
		if(Estado[i][0]==1)     //y si era defect
		{
		    Estado[i][0]=0;    //ahora es coop
		    Estado[i][1]=n;
		    S[i]=-1;
		}
		else             //y si era coop
		{
		    Estado[i][0]=1;       //ahora es defect
		    numvar[i]++;
		    FREC[i]=FREC[i]+(n-Estado[i][1]);
		   
		    S[i]=-1;
		    tiempo=(n-Estado[i][1]);
		    VentanasCoop[tiempo]=VentanasCoop[tiempo]+1.;
		    normaVentanasCoop=normaVentanasCoop+1.;
		    Estado[i][1]=0;
		}
	    }
	    
	    if(e[i]==0)         //si es coop
	    {	
		TC[i]=TC[i]+1;
	    }
	    
	}
	
	
    }  //fin pasos de tiempo
    


    n_coop_equil=n_coop_equil/jugadas_equil;      //media temporal de numero de cooperadores (sin discriminar)
    
    media_ben=media_ben/jugadas_equil;         //media temporal de los beneficios
    
    
    //printf("num coop. equil: %f  \n",  n_coop_equil);


   

    for(i=1;i<=N;i++)  //recuento numero cooperadores puros etc
	{
	    if(S[i]==0)
	    {
		cp++;
	    }
	    if(S[i]==1)
	    {
		dp++;
	    }
	    if(S[i]==-1)
	    {
		fl++;	    
	    }
	}
    //printf("cp:%f  dp:%f  fl:%f\n",cp,dp,fl);    
    
    for(i=1;i<=N;i++)     
    {
	
	if(TC[i]>0)                           //histograma tiempos cooperac. total
	{
	    if(TC[i]<jugadas_equil)   //evito los coop puros
	    {	
		Tiempos[TC[i]]=Tiempos[TC[i]]+1.;
		normatiempos=normatiempos+1.;
	    }
	}
    }
    
    /*for(i=1;i<=jugadas_equil;i++)     
    {
	printf ("Tiempos[%i]=%f \n",i,Tiempos[i]); 
	}*/


    for(i=1;i<=N;i++)
    {
	if(numvar[i]>0) 
	{
	    //printf ("numvar[%i] distinto de cero\n",i);
	    dobleprec=numvar[i];
	    FRECMED=FRECMED+FREC[i]/dobleprec;        // tiempo medio intervalos coop. (de los fluct)
	    NORMAFREC++;
	    frecmed[k[i]]=frecmed[k[i]]+FREC[i]/dobleprec;        //media tiempos coop. (de los fluct) segun su conectividad
	    normafrec[k[i]]++;	  
	}
    }

}



/////////////////////////////////////
/////////////////////////////////////
///////////////////////////////////////
//////////////////////////////////////
 
void TopologiaCP()    //Mi convenio: -1 fluct, 0 coop, +1 defect
{
     int J, I, IIII, ii;

    for(J=1;J<=N;J++)
    {
	NN[J]=0;  //NN[] es k[] pero solo enlaces entre cooperantes 
	for(IIII=1;IIII<=K_MAX;IIII++)
	{
	    Phase[J][IIII]=0;  //Phase[][] es C[][] pero solo entre cooperantes
	}
    }

    for(J=1;J<=N;J++)
    {
	if(S[J]==0) //Si es Coop. Puro
	{	
	    for(IIII=1;IIII<=k[J];IIII++)
	    {
		if(C[J][IIII]>J)
		{
		    if(S[C[J][IIII]]==0)
		    {
			NN[J]++;
			NN[C[J][IIII]]++;
			Phase[J][NN[J]]=C[J][IIII];
			Phase[C[J][IIII]][NN[C[J][IIII]]]=J;
		    }
		}
	    }
	}
	else
	{
	    NN[J]=0;
	}
    }

// Concluida la Matriz de Cooperantes Puros

    
    for(I=0;I<N;I++)
    {warning[I]=0;}   //Indica que se ha mirado si el nodo está en un cluster 
                      //de cooperadores con un 1 y que no se ha mirado con un 0
    GComp=0.;
    NCLUSTERS=0;

    //////

  
    for(I=1;I<=N;I++)
    {
	kk=0;

	for(J=1;J<=N;J++)
	{analyse[J]=0;}              //Indica los nodos que forman parte del 
	                             //cluster que se mira en ese momento

	if(NN[I]>0)
	{
	    if(warning[I]==0)
	    {	
		kk=1;
		warning[I]=1;
		analyse[kk]=I;
		for(J=1;J<=kk;J++)
		{ 
		    g=analyse[J]; 
		    for(ii=1;ii<=NN[g];ii++)
		    { 
			gg=Phase[g][ii];
			if(warning[gg]==0)
			{
			    warning[gg]=1;
			    kk=kk+1;
			    analyse[kk]=gg;
			} 
		    }
		}
		kk=kk+1; 
	    }
	}
	else
	{warning[I]=1;}

	if(kk>1)        //si ha encontrado un  cluster
	{
	    max=kk;
	    //printf("max (CP):%lf  GComp:%lf\n",max,GComp);
	    NCLUSTERS++;

	    if (max>GComp)    //guardo el tamaño de la componente gigante
		GComp=max;
	  	    
	}	

    }
   
     
    if (NCLUSTERS>0)  
	{
	    NonNullCP++;
	    NCLUSTERSMED=NCLUSTERSMED+NCLUSTERS;
	    GCMED=GCMED+GComp;
	}

    
}

    

//////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////

 
void TopologiaDP()
{
     int J, I, IIII,ii;
for(J=1;J<=N;J++)
    {
	NN[J]=0;  //NN[] es k[] pero solo enlaces entre defectores 
	for(IIII=1;IIII<=K_MAX;IIII++)
	{
	    Phase[J][IIII]=0;  //Phase[][] es C[][] pero sólo entre defectores
	}
    }

    for(J=1;J<=N;J++)
    {
	if(S[J]==1) //Si es Def. Puro
	{	
	    for(IIII=1;IIII<=k[J];IIII++)
	    {
		if(C[J][IIII]>J)
		{
		    if(S[C[J][IIII]]==1)
		    {
			NN[J]++;
			NN[C[J][IIII]]++;
			Phase[J][NN[J]]=C[J][IIII];
			Phase[C[J][IIII]][NN[C[J][IIII]]]=J;
		    }
		}
	    }
	}

	else
	{
	    NN[J]=0;
	}
    }

// Concluida la Matriz de Defectores Puros

  
    
    for(I=0;I<N;I++)
    {warning[I]=0;}   //Indica que se ha mirado si el nodo está en un cluster 
                      //de defectores con un 1 y que no se ha mirado con un 0
    
    GCompd=0.;
    NCLUSTERSd=0;

    //////

   

    for(I=1;I<=N;I++)
    {
	kk=0;
	for(J=1;J<=N;J++)
	{analyse[J]=0;}              //Indica los nodos que forman parte del 
	                             //cluster que se mira en ese momento

	if(NN[I]>0)
	{
	    if(warning[I]==0)
	    {
		kk=1;
		warning[I]=1;
		analyse[kk]=I;
		for(J=1;J<=kk;J++)
		{ 
		    g=analyse[J]; 
		    for(ii=1;ii<=NN[g];ii++)
		    { 
			gg=Phase[g][ii];
			if(warning[gg]==0)
			{
			    warning[gg]=1;
			    kk=kk+1;
			    analyse[kk]=gg;
			} 
		    }
		}
		kk=kk+1; 
	    }
	}
	else
	{warning[I]=1;}

	if(kk>1)
	{
	    max=kk;
	    //printf("max (DP):%lf  GCompd:%lf\n",max,GCompd);
	    NCLUSTERSd++;

	  if (max>GCompd)    //guardo el tamaño de la componente gigante
		GCompd=max;
	    
	   
	}	
    }
   
     
    if (NCLUSTERSd>0)  
	{
	    NonNullDP++;
	    NCLUSTERSMEDd=NCLUSTERSMEDd+NCLUSTERSd;
	    GCMEDd=GCMEDd+GCompd;
	}
    
    //printf("GC: %lf\n NClusters: %i\n",GComp*(float)NODOS,NCLUSTERS);
    
    

   
    
    
    //printf("GCompDP:%lf  NClustersd: %i\n",GCompd,NCLUSTERSd);  
 
    
}



void histograma_rho_k ()
{
  int i;
  
  for(i=1;i<=K_MAX;i++)    
    {
      rho_C_k[i]=0.;
      rho_D_k[i]=0.;
      rho_Fl_k[i]=0.;
    }
  
  for(i=1;i<=N;i++)
    {
      
      if(S[i]==0)   //si es coop puro
	{
	  rho_C_k[k[i]]++;	
	}
      else 
	if(S[i]==1)  //si es defect puro
	  {
	    rho_D_k[k[i]]++;	 
	  }
	else     //si es fluct
	  {
	    rho_Fl_k[k[i]]++;	   
	  }
    }
  
  for(i=1;i<=K_MAX;i++)    // lo normalizo al numero de nodos con esa conectividad en esa realizacion
    {
      if(PK[i]>0.1)
	{
	  rho_C_k[i]=rho_C_k[i]/PK[i];
	  rho_D_k[i]= rho_D_k[i]/PK[i];
	  rho_Fl_k[i]= rho_Fl_k[i]/PK[i];
	  iter_k[i]++;
	}
      else 
	{
	  rho_C_k[i]=0.;
	  rho_D_k[i]= 0.;
	  rho_Fl_k[i]=0.;
	}
    }
  
}
