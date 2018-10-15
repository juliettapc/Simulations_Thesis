//Implementacion de una dinamica de juego combinada con el crecimiento de la red por pref.-attach

// bucle de estadística y de barrido en b, ademas de computar para el p(k) solo las realizaciones no nulas   

//actualizacion paralela de las estrategias

//calculo de average path lenth y de clustering coeficient


//(basado en el growth&play-v2_validas.c)


// continuacion de la dinamica del juego una vez acabada de crecer la red.

//CORRECCION en la normalización: diferencia máxima de coef de la matriz de payoff*max k (antes: b*k_m siempre)
		



# include <stdio.h>
# include <stdlib.h>
# include <math.h>


# define N 1000 //tamaño de la red
# define m_o  2      //nodos inicialmente unidos
# define m   2      //links nuevos añadidos a cada paso de tiempo

# define Niter 50  //estadistica
# define K_MAX   1000      //para el tamaño de C[N][K_MAX]     ojo pq depende de N en SF!!!
                          //tambien calculable en la funcion matriz_Conect

# define  eps 0.2

# define  tauT 1      //(si =10, pongo 10 nodos y juego una vez)
# define  tauD 1

# define  ro 0.5
# define  R 1.0
# define  Pu 0.0
# define  Su 0.0

# define b_min 1.5
# define b_max 1.55

# define delta_b 0.2


# define jugadas_extra   10000       //despues de acabar de crecer la red

# define jugadas_equil   20000       //para discriminar CP,DP y Fl



//Generacion de numeros aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];




 
int k[N+1],k_PA[N+1],A[N+1];  //conectividad topológica y de Preferential Attachment (PA)
int   M[N+1][m+1], unido[m+1], Cuidado[N+1]; //  unido[] guarda a quién se ha unido, y cuidado[] quien le ha lanzadolinks

double CC[K_MAX+1],NC[K_MAX+1],PC[K_MAX+1];

double P[N+1], P_prov[N+1];     //P_prov[] es la provisional para manipularla
double PK[N+1], PK_tot[N+1], PKa[N+1], PKa_tot[N+1], PK_tot_val[N+1], PKa_tot_val[N+1];    //  para la distribucion P(k)
double r,v;                     //para guardar los aleatorios
int  II,i,j,jj,g,gg, d, q,y,z, x[N+1], C[N+1][K_MAX+1], steps,iter;
double norma,norma2,Normandia;
int si[N+1], cont,cont2;

//Variables para el dilema

//double b,R,Pu,Su;
//int tauT,tauD;
int tiempo;

int  e[N+1],e_aux[N+1], t, w, k_m; 
double  ben[N+1],p;
int  n_coop; //numero de coop. instantaneos
int iter,counter,n;

int s, pasos,ii;


double fit[N+1],c_media;       
double b;
double aux;
int validas,Validas;


    //para el av_path_l
double dist_tot,D_med;
    	
int D[m*N+1][m+1];   //matriz de pares de links

int D_prov[m*N+1][m+1];          //para el clust_coef

double clust_coef_tot;   

double n_coop_equil;
int cp,dp,fl;
double cp_media,dp_media,fl_media;    //id pero sobre las Niter realizacines
int S[N+1],Estado[N+1][2];      //etiqueta para discriminar entre puros y fluctuantes

double NCLUSTERSMED,NCLUSTERSMEDd,GCMED, GCMEDd,GCompd,GComp, CLUSTERS,CLUSTERSd;

int NN[N+1],Phase[N+1][K_MAX+1],warning[N+1],analyse[N+1],semilla[N+1];
int NonNull,NonNullCP,NonNullDP,NCLUSTERS,NCLUSTERSd;


void inicia_rand(int semilla);
void nucleo_inicial();
void jugada();
void discriminar();
void nuevo_nodo();
void histograma_pk();
void matriz_Conect();
void av_path_length();
void clustering_coef();
void matriz_D();
void TopologiaCP();
void TopologiaDP();


char nombre[256],nombre3[256],file8[256] ,file9[256], nombre4[256];
char file7[256]; 
char file[256], file1[256], file4[256],file44[256];  
char file2[256],file3[256],file77[256];

FILE *fich7;                   
FILE *fich; 
FILE *fich2; 
FILE *fich3; 
FILE *fich1; 
FILE *fich4;
FILE *fich77;      
FILE *fich8; 
FILE *fich44; 




int main()
{
  
/*  sprintf(file3,"%dC_med%d-%d-b_ini%.2lf-e%.2lf_%diter.dat",N,tauT,tauD,b_min,eps,Niter);   //un solo archivo     
  fich3=fopen(file3,"wt");
  fclose(fich3);*/
  

  sprintf(file4,"%dProp-top_%d-%d-b_ini%.2lf-e%.2lf_%diter_s.dat",N,tauT,tauD,b_min,eps,Niter);   //un solo archivo     
  fich4=fopen(file4,"wt");
  fclose(fich4);
  


  //Guarda los Clusters de Cooperantes etc.
  sprintf(file8,"%dClusters_%d-%d-b_ini%.2lf-e%.2lf_%diter_s.dat",N,tauT,tauD,b_min,eps,Niter);
  fich8=fopen(file8,"wt");
  fclose(fich8);







  /* sprintf(file2,"%dLength%d-%d-b%.2lf-e%.2lf_pay.dat",N,tauT,tauD,b_min,eps);        
     fich2=fopen(file2,"wt");
     fclose(fich2);*/
  
  	  inicia_rand(time(0));
  
  printf("\n\nPrograma Growth&Play + Seguir_jugando (%d iteraciones)\nN=%d  b_ini=%lf b_fin=%lf eps=%lf  TauT=%d  TauD=%d\n", Niter,N,b_min,b_max,eps,tauT,tauD);
  
  
  b=b_min;
  
  while(b<=b_max)     //bucle externo para barrer en b
    {
      printf("\nb=%lf\n",b);
      


  
      sprintf(file7,"%dPk%d-%d-b%.2lf-e%.2lf_%diter_s.dat",N,tauT,tauD,b,eps,Niter);   //un p(k) para cada b     
      fich7=fopen(file7,"wt");
      fclose(fich7);

      sprintf(file77,"%dCk%d-%d-b%.2lf-e%.2lf_%diter_s.dat",N,tauT,tauD,b,eps,Niter);   //un p(k) para cada b     
      fich77=fopen(file77,"wt");
      fclose(fich77);

      sprintf(file,"%dEvolCoop%d-%d-b%.2lf-e%.2lf_%diter_s.dat",N,tauT,tauD,b,eps,Niter);    // id.    
      fich=fopen(file,"wt");
      fclose(fich);

      sprintf(file1,"%dnum_coop%d-%d-b%.2lf-e%.2lf_%diter_s.dat",N,tauT,tauD,b,eps,Niter);    // id.    
      fich1=fopen(file1,"wt");
      fclose(fich1);




      c_media=0.;
      Normandia=0.;

      for (i=0;i<=N;i++)
	{
	  PK_tot[i]=0.;
	  PKa_tot[i]=0.;
	  PK_tot_val[i]=0.;
	  PKa_tot_val[i]=0.;
	}

      for(i=0;i<=K_MAX;i++)
	{
	  CC[i]=0.;
	  NC[i]=0.;
	  PC[i]=0.;
	}
      
      
      
      NonNullCP=0;
      NonNullDP=0;
      NCLUSTERSMED=0.;
      NCLUSTERSMEDd=0.;
      GCMED=0.;
      GCMEDd=0.;
      Validas=0;
      
      
      clust_coef_tot=0.0;
      D_med=0.;
      validas=0;
      c_media=cp_media=dp_media=fl_media=0;

      for(iter=1;iter<=Niter;iter++)         //inicio bucle de estadistica
	{
	  
	  printf("\niter:%d\n",iter);
	  
	  
	  
	  // INICIAliZO:
	  
	  steps=N-m_o;       //numero de nodos que se añadiran en total
	  s=m_o;         //tamagno actual de la red: s=N(t)
	  tiempo=1;       //contador de tiempos para tauT y tauD...	
	  n_coop= 0;  
	  
	  for(i=1;i<=N;i++)      
	    {
	      ben[i]=0;
	      e[i]=1;     //TODOS DEFECTORES      
	    }
	  

	  
	  nucleo_inicial();
	  //    matriz_Conect();  //hay que actualizarla cada vez para poder jugar
	  //printf("num coop.: %d  \n",  n_coop);
	  
	  
	  do{
	    
	    for(II=1;II<=tauD;II++) 
	      {
		jugada();
		
		n_coop=0;
		for(i=1;i<=s;i++)
		  {
		    if(e[i]==0)	    
		      n_coop++;	      
		  }
		fich=fopen(file,"at");
		fprintf(fich,"%d   %d\n",s,n_coop);      
		fclose (fich);
		
	      }
	    
	    //Recalculo Probabilidad de attachment: (SOLO AQUI PARA OPTIMIZAR)
	    
	    for(i=1;i<=N;i++)
	      {P[i]=0.;}
	    
	    for(i=1;i<=s;i++)
	      {
		counter=counter+k[i];
	    for(j=1;j<=i;j++)  
	      {
		P[i]=P[i]+fit[j]; 
	      }
	      }
	    norma=P[s];
	    
	    //ATTACHMENT:
	    
	    for(II=1;II<=tauT;II++)
	    {
		if(s<N)               //para evitar un segmentation fault, en s=1002 cuando tauT=10!!!
		{
		    s++;	    
		    nuevo_nodo();	      
		    //printf("%i   %lf\n",s,norma);	
		    //	    matriz_Conect();
		}
		
	    }
	    
	    //printf("num coop.: %d  \n",  n_coop);
	    
	    tiempo=tiempo+1; 
	    
	  }while(s<N);        //acaba de crecer la red
	  
	  
	  
	  matriz_D();                   //calculo av_path_leng y clustering  (acumulacion del valor al final de las funciones)
	  av_path_length();
	  clustering_coef();
	  

	  
	  if(n_coop>4*m) //para la norma del p(k)
	  {
	      validas++;   
	  }
	  
	  histograma_pk(); 
	  
	  
	  
	  for(ii=1;ii<=jugadas_extra;ii++)     //bucle para continuar jugando con la red ya crecida
	  {
	      jugada();
	      n_coop=0;
	      for(i=1;i<=s;i++)
	      {
		  if(e[i]==0)	    
		      n_coop++;	      
	      }
	      fich=fopen(file,"at");
	      fprintf(fich,"%d   %d\n",s+ii,n_coop);      
	      fclose (fich);
	      
	  }

	  
//	  printf("n_coop=%i\n",n_coop);
	  


	  
	  //guardo el estado tras el transitorio como referencia para discriminar
	  

	  for(i=1;i<=N;i++)      //inicializo variables locales (mediaran en una iteracion)
	  {
	      
	      Estado[i][0]=1;  //los inicializo a defectores
	      Estado[i][1]=0;	
	  }
	  
	  
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
	  
	 
	  //  printf("fin del transitorio, comienzan jugadas de equilibrio\n");
	  
	  cp=dp=fl=n_coop_equil=0.0;    
	  
	  for(ii=1;ii<=jugadas_equil;ii++)     //bucle para alcanzar el equilibrio y discriminar
	  {
	      jugada();
	      n_coop=0;
	      for(i=1;i<=s;i++)
	      {
		  if(e[i]==0)	    //recuento de cooperantes instantaneos (medios, sin discriminar)
		      n_coop++;	      
	      }
	      fich=fopen(file,"at");
	      fprintf(fich,"%d   %d\n",s+jugadas_extra+ii,n_coop);      
	      fclose (fich);
	      
	      
	      n_coop_equil+=n_coop;	
	      
	      
	      
	      // printf("n_coop_equil=%lf\n",n_coop_equil);
	      
	      
	      discriminar();
	      
	      
	      
	  }                                      //fin jugadas equil
	  
	  n_coop_equil=n_coop_equil/jugadas_equil;
	  
	  
	  for(i=1;i<=N;i++)  //recuento numero CP, DP y Fl
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
	  
	  
	  printf("cp:%d  dp:%d  fl:%d\n",cp,dp,fl);  
	  
	  
	  c_media+=n_coop_equil;
	  cp_media+=cp;
	  dp_media+=dp;
	  fl_media+=fl;
	  
	  
	  
	  TopologiaCP();
	  TopologiaDP();



	  
	  aux=n_coop/(double)N;
	  fich1=fopen(file1,"at");	  
	  fprintf(fich1,"%d  %lf\n", iter,aux);     //guardo la coop_media de cada iter	  	  
	  fclose (fich1);      
	  
	  
	}                               //fin bucle estadistica
      
      
      
      
      c_media=c_media/(double)(Niter*N);            //normalizar, y tanto por uno
      
      cp_media=cp_media/(double)(Niter*N); 
      dp_media=dp_media/(double)(Niter*N); 
      fl_media=fl_media/(double)(Niter*N);
      
      
      clust_coef_tot=clust_coef_tot/Niter;
      D_med=D_med/Niter;
      
      printf("\n\nb:%lf   clust_coef:%lf   Av-path_l:%lf   \n c_med:%lf  cp_med:%lf  dp_med:%lf   fl_med:%lf \n",b, clust_coef_tot, D_med, c_media,cp_media,dp_media,fl_media);
      
      //  Normalizo pk y pk_acumulada
      
      for (i=1;i<=N;i++)
      {
	  PK_tot[i]=PK_tot[i]/Niter;
	  PKa_tot[i]=PKa_tot[i]/Niter;
	  
	  if(validas!=0)
	  {	  
	      PK_tot_val[i]=PK_tot_val[i]/validas;
	      PKa_tot_val[i]=PKa_tot_val[i]/validas;
	  }
	  
	  else
	  {
	      PK_tot_val[i]=0;
	      PKa_tot_val[i]=0;
	  }
	  
      }
      
      printf("Niter:%d   Validas:%d\n",Niter, validas);
      
      // Escribo PK  
      
      fich7=fopen(file7,"at");
  
      for(i = 1; i<N; i++)
	{
	  fprintf(fich7,"%d  %lf  %lf  %lf  %lf\n", i, PK_tot[i],PKa_tot[i], PK_tot_val[i],PKa_tot_val[i]);
	}                   //normalizo sobre iters y tb sobre realizaciones con coop no nula
      
      fclose (fich7);
      


      fich77=fopen(file77,"at");
  
      for(i = 1; i<=K_MAX; i++)
	{

	  CC[i]=CC[i]/NC[i];
	  PC[i]=PC[i]/NC[i];

	  fprintf(fich77,"%d  %lf  %lf\n", i, CC[i],PC[i]);
	}                   //normalizo sobre iters y tb sobre realizaciones con coop no nula
      
      fclose (fich77);
      
      
      /*escribo la c_media
      
      fich3=fopen(file3,"at");
      fprintf(fich3,"%lf   %lf   %lf   %d   %d\n",b,c_media,eps,tauT, tauD);      
      fclose (fich3);*/


//escribo  el clust-coef , el Av-path_length y la c_med

      fich4=fopen(file4,"at");
      fich4=fopen(file4,"at");
      fprintf(fich4,"%lf   %lf   %lf   %lf     %lf   %lf   %lf   %lf   %d   %d\n",b, clust_coef_tot, D_med, c_media, cp_media,dp_media,fl_media,eps, tauT, tauD);      
      fclose (fich4);

      
      
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
      
      
      fich8=fopen(file8,"at");
      fprintf(fich8,"%lf   %lf   %lf   %lf   %lf   %i   %i   %i     %lf   %i   %i\n",
	      b,GCMED,GCMEDd,NCLUSTERSMED,NCLUSTERSMEDd,Validas,
	      NonNullCP,NonNullDP,eps,tauT,tauD);
      
      fclose(fich8);
      
      
      printf("GCMED=%.2lf   GCMEDd=%.2lf   NCLUSTERSMED=%.2lf   NCLUSTERSMEDd=%.2lf   NonNullCP=%i   NonNullDP=%i\n",GCMED,GCMEDd,NCLUSTERSMED,NCLUSTERSMEDd,NonNullCP,NonNullDP); 















      b+=delta_b;    
  //fclose(escribe);
  //fclose(escribe2);
      
    }     //fin del bucle en b
  
  


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


void nucleo_inicial()    //conecto los m_o nodos iniciales todos entre si
{


    for(i = 0; i<=N; i++)
    {
	k[i] = 0;     
      
      for(j = 0; j <= m_o; j++)
      {M[i][j] =0;}
 
      for(j=0;j<=K_MAX;j++)
      {C[i][j]=0;}
    }
  
  
    for(j = 1; j <=m_o ; j++) 
	k[j]=m_o-1;

    //  C[1][1]=2;
    //  C[2][1]=1;

    for(i = 1; i <=m_o ; i++)
    {
	counter=0;
	for(j = 1; j <=m_o ; j++)
	{
	    if(i!=j)
	    {
		counter++;
		C[i][counter]=j;
	    }
	}
    }
  

    for(i = 2; i <=m_o ; i++)
    {
	for(j = 1; j <=i-1 ; j++)
	{
	  M[i][j] = j;
	  //printf("M[%d][%d] = %d \n",i,j,M[i][j]);
	}
    }
  
    for(i=1;i<=m_o;i++)                 
    { 
	fit[i]=1.-eps;  
    }

  for(i = 0; i <= N; i++)
  {
      P[i] = 0;
      P_prov[i]=0;
  }
  
  for(i = 1; i<= m_o; i++)                    //prob de los nodos iniciales
  { 
      
      for(j=1; j<=i; j++)
      {                              
	  P[i] = P[i] + fit[j];        //para los m_o fit[i]=1
      }
      
      P_prov[i]=P[i];
  }
  
  norma=P[m_o];
  norma2=norma;

 // establezco como cooperadores a los  nodos iniciales:
  

  for(i=1;i<=m_o;i++)
    {     
      e[i]=0;       //cooperador
      n_coop++;
    }
  
  
}


//////////////////////////////////////
//////////////////////////////////////



void nuevo_nodo() 
{
  
// establezco el carcter del nuevo nodo:
  
  r=FRANDOM;
  
  if(r<ro)
    {                
      e[s]=0;       //cooperador
    }


  for(i=0; i<=m; i++)
    {unido[i]=0;}
  
  for(i=1; i<=N; i++)
    {P_prov[i]=P[i];}   //guardo la prob en un vector aux para manipularlo
  
  norma2=norma;
  P[s]=0.;  
  P_prov[s]=0.;  
  
  
  for(q = 1; q<=m; q++)        //bucle a los m links nuevos que voy a añadir
    {
      
      v=FRANDOM*norma2;
     
      for(j=1; j<s; j++)
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
      
      for(j=g+1; j<s; j++)  
      {
	  P_prov[j]=P_prov[j]-fit[g] ; 
      }
      
      norma2=norma2-fit[g];
      
      //printf(" norma:   %d \n",norma2);
    }    //fin del bucle sobre los m links lanzados
  
  
  
  for(j=1;j<=m_o;j++)
    {
	g=unido[j];
	k[g]=k[g]+1;
	C[g][k[g]]=s;
	C[s][j]=g;
	M[s][j]=g;
        
        /*
	  for(i=g;i<s;i++)   
	  {
	  P[i]=P[i]+fit[g];
	  }
	*/
    }
  
  k[s]=m;
  fit[s]=(1.-eps);
  P[s]=P[s-1]+fit[s];
  norma=P[s];
 
 
}



void matriz_Conect()   //obtencion de la matriz C a partir de la M
{


    
    for(i=0;i<=N;i++)
      {
	x[i]=0;
	
	for(j=0;j<=K_MAX;j++)
	  {
	    C[i][j]=0;
	  }
      }
    
    for(i=1;i<=s;i++)     //ahora solo va de 1 a s
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
    
}


//////////////////////////////////////////
////////////////////////////////////

void jugada()  //es un paso del bucle de tiempo de "juego()" 
{

    n_coop=0;

    for(i=1;i<=s;i++)         // ahora de 1 a s!!
    {
	ben[i]=0;
    }
    
    for(i=1;i<=s;i++) // cada nodo juega con cada vecino ahora de 1 a s!!
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
    }            //fin de la partida


for (i=1;i<=s;i++)
    {
	e_aux[i]=e[i];       //guardo la estrategia del ultimo paso (para el update paralelo)
    }
    

    
  for(i=1;i<=s;i++)        //comparo estrategias    ahora de 1 a s!!
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
	      e_aux[i]=e[t];
	    
	  }
       }

     }



 for (i=1;i<=s;i++)      //actualizacion de estrategias en paralelo
      {
	e[i]=e_aux[i];
      }
    


  //calculo del fitness

  for (i=1;i<=s;i++)
    {
      fit[i]=1.0-eps+eps*ben[i];
    }
       
}

/////////////////////////////////////
///////////////////////////////////////////


void histograma_pk()            //costruccion de P(k)



{
  for(i = 0; i <= N; i++)           //inicializo
    {
	PK[i]=0.;
	PKa[j]=0.;
    }
  
  for(i =1; i <= N; i++)            //recuento
  {
      PK[k[i]]++;

      for(j=1;j<=k[i];j++)
      {
	  PKa[j]++;
      }
  }

  for(i = 1; i <= N; i++)         //normalizo
    {
      PK[i] = PK[i]/N;    
      PKa[i] = PKa[i]/N;      
     
      PK_tot[i] += PK[i];    
      PKa_tot[i] += PKa[i];  

      if(n_coop>4*m)   //recuento solo de las iteraciones validas
	{
	  PK_tot_val[i] += PK[i];    
	  PKa_tot_val[i] += PKa[i]; 
	} 

   
    }
  
 
}        

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

  

void matriz_D()        //contiene todos los pares de vecinos (segun orden de creacion del link)
{
  int i,j,h, n_links;
  
  for(i=0;i<=2*N;i++)
    {
      for(j=0;j<=m;j++)
	{
	  D[i][j]=0;
	}
    }
  
  n_links=0;
  h=1;
  for(i=1;i<=N;i++)
    {
      for(j=1;j<=m;j++)
	{
	  if(M[i][j]!=0)
	    {
	      D[h][1]=i;
	      D[h][2]=M[i][j];
	      h++;
	    }
	}
      
    }
  n_links=h;   //guardo el numero total de links


/////guardo en un fichero la red para pintarla


  sprintf(file44,"%dRED_%d-%d-b_ini%.2lf-e%.2lf_s.NET",N,tauT,tauD,b_min,eps);
  fich44=fopen(file44,"wt");
  fclose(fich44);
  
  
  fich44=fopen(file44,"at");

 fprintf(fich44,"*Vertices      %d\n",N); 

 fprintf(fich44,"*Edges"); 

  for(i=1;i<=2*N;i++)
  {
      fprintf(fich44,"\n");
      for(j=1;j<=m;j++)
      {	  
	  if(D[i][1]!=0)
	  {
	  fprintf(fich44,"%d   ",D[i][j]); 
	  }
	  
      }
  }
  fclose (fich44);
  
  /* printf("\n");
  for(i=1;i<=2*N;i++)
    {
      printf("\nD[%d][]: ",i);
      for(j=1;j<=m;j++)
	{
	  printf("%d  ",D[i][j]);
	}
    }
  // getchar();
  
  printf("fin funcion matriz_D\n"); */
}  




//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



void clustering_coef()
{
    int i,j,h,n,Nvalidos,links,nodo1,nodo2,vecino[K_MAX+1],conectividad;
    double coef, clust_coef,doubleprec, combina,aux1;
    
    for(i=1;i<=2*N;i++)
	for(j=1;j<=2;j++)
	    D_prov[i][j]=D[i][j];   //guardo la matriz de pares de links para manipularla
    
    
    
    Nvalidos=0;
    clust_coef=0.0;

    for(i=1;i<=N;i++)    //recorro la red
    {
//	printf("\nnodo %d  con conectividad %d\n",i,k[i]); 
	
	conectividad=k[i];
	
	if(conectividad>1)  
	{
	    Nvalidos++;
	    

	    links=0;
	    coef=0.0;
	    for(j=1;j<=conectividad;j++) 
	    {
		vecino[j]=C[i][j];    //guardo los vecinos del nodo i
	    }
	    
	    
	    for(j=1;j<=conectividad;j++)     //recorro los vecinos de i
	    {	      
		nodo1=vecino[j];	   
		
		
		for(n=j+1;n<=conectividad;n++)   //busco vecinos de i que lo sean entre si
		{
		    // printf("j:%d n:%d\n",j,n);  
		    nodo2=vecino[n];
		    
		    if(nodo1 != nodo2)
		    {
			
			for(h=1;h<=m*N;h++)   //recorro la matriz de pares de links
			{
			    if((D_prov[h][1] !=0)   &&  (D_prov[h][2] !=0) )
			    {
				
				if( nodo1==D_prov[h][1] && nodo2==D_prov[h][2]) 
				{
				    links ++;				    
				}
				else if(nodo1==D_prov[h][2] && nodo2==D_prov[h][1])
				{
				    links ++;				    
				}
			    }
			}
			
		    }

		   

		}	  //fin de la condicional k[i]>1 
		
	    }      //recorridos todos los vecinos de i
	    
	    
	    aux1=conectividad;
	    
	    combina=aux1*(aux1-1.0)/2.0;   //maximo de links que podria haber entre los vecinos de i
	    
	    doubleprec=links;
	    
	    coef=doubleprec/combina;
	    
	    CC[conectividad]=CC[conectividad]+coef;
	    NC[conectividad]=NC[conectividad]+1.;

	    if(e[i]==0)	    
	      {
		PC[conectividad]=PC[conectividad]+1.;
		Normandia=Normandia+1.;
	      }	      
		  	    
	    
	    clust_coef=clust_coef+coef;
//	printf("links:%d  combina:%lf  coef:%lf  clust_coef:%lf\n",links,combina,coef,clust_coef); 	
	    
	}
	
    }//recorrida toda la red
    
    
    
    
    clust_coef=clust_coef/Nvalidos;
    printf("clust_coef:%lf\n",clust_coef); 

    
    clust_coef_tot+=clust_coef;	           //voy acumulando (estadistica)

/*
/////////////// GC //////////////////////////
/////////////////////////////////////////////


   for(I=1;I<=N;I++)
    {
	warning[I]=0;
	analyse[I]=0;
    }  

   GComp=0.;
   NCLUSTERS=0;

   for(I=1;I<=N;I++)
   {
	kk=0;
	
	for(J=1;J<=N;J++)
	{
	    analyse[J]=0;
	} 

	if(k[I]>0)
	{
	    if(warning[I]==0)
	    {
		kk=1;
		warning[I]=1;
		analyse[kk]=I;

		for(J=1;J<=kk;J++)
		{ 
		    g=analyse[J]; 
		    for(i=1;i<=k[g];i++)
		    { 
			gg=C[g][i];

			if(warning[gg]==0)
			{
			    warning[gg]=1;
			    kk=kk+1;
			    analyse[kk]=gg;
			} 
		    }
		} 
	    }
	}
	else
	{warning[I]=1;}
	
	if(kk>0)
	{
	    NCLUSTERS++;
	    if(GComp<kk)
	    {
		GComp=kk;
		//Imax=analyse[1];
	    }
	}
    }

    printf("NUMERO DE CLUSTERS BI %i\n",NCLUSTERS);
    printf("GC BI %i\n",GComp);
*/
    
}






/////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


void av_path_length()
{


    int i,j,  h,cont,dist;
    int actual[N+1],di[N+1];  
      
      dist_tot=0.;
      
      for (i=1;i<=N;i++)     // para todos los nodos de la red
	{
	  dist=0;    //guardara la suma de las distancias de un nodo a todos los demas

	  for (j=1;j<=N;j++)
	    {
	      di[j]=-100;   //inicializo las banderas de los nodos como "no visitado"
	    }	               //guardara la distancia desde el nodo i a todos los demas

	  cont=1;
	  actual[cont]=i;     //guardo el nodo que estoy mirando en este momento (empiezo por mi)
	  di[i]=0;   //para que si parto de mi, no vuelva a mi mismo de nuevo

	  for(j=1;j<=cont;j++)
	    {	   

	      for(h = 1; h <= k[actual[j]]; h++)     //bucle sobre la conectividad del nodo que estoy mirando
		{
		  if(di[C[actual[j]][h]]==-100)   //si no he visitado el nodo
		    {
		      cont++;
		      actual[cont] = C[actual[j]][h]; 
		      di[C[actual[j]][h]]=di[actual[j]]+1;   //lo marco con 1+ la distancia del nodo actual
		    }
		}
	    }
	  


	  for (j=1;j<=N;j++)
	    {
	      if (di[j]>0)
		{
		  dist+=di[j];
		  //printf("dist:%d\n",dist);
		}
	    }
	  
	  dist_tot+=dist;
	  //printf("dist_tot:%lf\n",dist_tot);
	  
	}// fin del bucle a todos los nodos de la red
      


      dist_tot=dist_tot/((N-1)*N);



      D_med+=dist_tot;           //voy acumulando (estadistica)



printf("dist_tot:%lf\n",dist_tot);

   



}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



void discriminar()
{


    int i;

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
		
		S[i]=-1;
		
		Estado[i][1]=0;
	    }
	}
	
    }
    



}


///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


 
void TopologiaCP()    //Mi convenio: -1 fluct, 0 coop, +1 defect
{
     int J, I, IIII, ii,max,kk;

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

    printf("Ncc:%i  GCc:%lf\n",NCLUSTERS,GComp);
}

    

//////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////

 
void TopologiaDP()
{
     int J, I, IIII,ii,max,kk;
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
    printf("Ncd:%i  GCd:%lf\n",NCLUSTERSd,GCompd);
    //printf("GC: %lf\n NClusters: %i\n",GComp*(float)NODOS,NCLUSTERS);
    
    

   
    
    
    //printf("GCompDP:%lf  NClustersd: %i\n",GCompd,NCLUSTERSd);  
 
    
}


