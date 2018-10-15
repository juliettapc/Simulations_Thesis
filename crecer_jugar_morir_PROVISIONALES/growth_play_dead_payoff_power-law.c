//8-2-2010: Programa para simular una variacion del modelo EPA, que combinaba
// crecimiento y dinamica del DP, pero ahora ademas, añadimos muerte de los 
// nodos como funcion de su payoff, de modo que los más ricos, viven mas

//basado en el programa "growth_play_tot_Pkct.c"

//(pero) convenio de estrategias:  1==cooperador, 0==defector

//la probabilidad de morir disminuye con el payoff-->>> método Jesus de Pref. Attach.



# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

# define N 1000 //tamaño de la red
# define m_o  2      //nodos inicialmente unidos
# define m   2      //links nuevos añadidos a cada paso de tiempo

# define Niter 10   //estadistica




# define  eps 0.99

# define  tauT 10       //(si tauT=10 y tauD=1, pongo 10 nodos y juego una vez)
# define  tauD 1
 

# define  ro 0.5 
# define  R 1.0
# define  Pu 0.0
# define  Su 0.0


# define b_min    1.00
# define b_max   3.00
# define delta_b 0.1



# define jugadas_evolucion 3000      //num pasos tiempo durante los cuales juega-crece-muere, 
                                     // tras haber alcanzado el tamaño N


//////Generacion de numeros aleatorios/////////////////
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
/////////////////////////////////////////////////////////



int k[N+1],k_PA[N+1],A[N+1];  //conectividad topológica y de Preferential Attachment (PA)
double norma,norma2;
int   M[N+1][m+1], unido[m+1], Cuidado[N+1]; //  unido[] guarda a quién se ha unido, y cuidado[] quien le ha lanzadolinks
double P[N+1], P_prov[N+1];     //P_prov[] es la provisional para manipularla
int  C[N+1][N+1], x[N+1];
int  C_aux[N+1][N+1];
double PK[N+1], PK_tot[N+1], PK_coop[N+1], PK_coop_tot[N+1], PK_tot_fin[N+1], PK_coop_tot_fin[N+1];




int  e[N+1],e_aux[N+1], k_m; 
double  ben[N+1],p;
int  n_coop; //numero de coop. instantaneos
int iter,counter,n,ii;

int s,II,i,j,hueco,tpo;


double  fit[N+1],c_media,c_media_fin;       
double b;



void inicia_rand(int semilla);
void nucleo_inicial();
void jugada();
void nuevo_nodo();
void muerte_nodo();
void nuevo_nodo_hueco();
void histograma_pk();
void histograma_pk_coop();
void matriz_Conect();



char file1[256], file2[256], file3[256], file4[256];

FILE *fich1;     
FILE *fich2;   
FILE *fich3;   
FILE *fich4;   


int main()
{
  
  inicia_rand(time(0));
  







  printf("\n\nPrograma Growth & Play & Dead (payoff) (%d iter)\nN=%d  b_ini=%lf b_fin=%lf eps=%lf  TauT=%d  TauD=%d\n", Niter,N,b_min,b_max,eps,tauT,tauD);
  
  
  
  
  
  sprintf(file1,"C_vs_b_N%d_%d-%d-b_ini%.2lf-e%.2lf_tpo%d_%diter_payoff.dat",N,tauT,tauD,b_min,eps,jugadas_evolucion,Niter);  
  fich1=fopen(file1,"wt");
  fclose(fich1);
  
  
  
  
  
  
  b=b_min;    
  while(b<=b_max)     //bucle externo para barrer en b
    {
      printf("\nb=%lf\n",b);
      
      
      
      
      sprintf(file2,"Pk_Pkc_N%d_%d-%d-b%.2lf-e%.2lf_tpo%d_%diter_payoff.dat",N,tauT,tauD,b,eps,jugadas_evolucion,Niter);      
      fich2=fopen(file2,"wt");
      fclose(fich2);   
      
      

      
      sprintf(file3,"EvolCoop_N%d_%d-%d-b%.2lf-e%.2lf_%diter_payoff.dat",N,tauT,tauD,b,eps,Niter);    // id.    
      fich3=fopen(file3,"wt");
      fclose(fich3);
      
      
      
      sprintf(file4,"Conex_N%d_%d-%d-b%.2lf-e%.2lf_tpo%d_%diter_payoff.dat",N,tauT,tauD,b,eps,jugadas_evolucion,Niter);      
      fich4=fopen(file4,"wt");
      fclose(fich4);   
      
      for(i=1; i<=N; i++)       
	{
	  PK_tot[i]=0;	 
	  PK_coop_tot[i]=0;

	  PK_tot_fin[i]=0;	 
	  PK_coop_tot_fin[i]=0;
	}
      
      
      c_media=0.;
      c_media_fin=0.;
      for(iter=1;iter<=Niter;iter++)         //inicio bucle de estadistica
	{
	  printf("\niter:%d\n",iter);
 	  
	  
	  s=m_o;         //tamagno actual de la red: s=N(t)	  	 
	  
	  for(i=1;i<=N;i++)      
	    {
	      ben[i]=0;
	      e[i]=0;     //TODOS DEFECTORES      
	    }
	  
	  
	  
	  nucleo_inicial();
	  matriz_Conect(); 
	  
	  do{
	    
	    for(II=1;II<=tauD;II++) 
	      {
		jugada();
		
		//		printf("jugada  ");	      	
		n_coop=0;
		for(i=1;i<=s;i++)
		  {
		    if(e[i]==1)	    
		      n_coop++;	      
		  }
		//		    printf("n_coop:%d\n",n_coop);		    


		if(iter<=5)    //pq si no, el archivo se hace demasiado grande!!
		  {
		    fich3=fopen(file3,"at");
		    fprintf(fich3,"%d   %d\n",s,n_coop);      
		    fclose (fich3);
		  }

		
	      }


	    //Recalculo Probabilidad de attachment: 	    
	    for(i=1;i<=N;i++)
	      {P[i]=0.;}
	    
	    for(i=1;i<=s;i++)
	      {
		for(j=1;j<=i;j++)  
		  {
		    P[i]=P[i]+fit[j]; 
		  }
	      }
	    norma=P[s];
	    //printf("recalculadas prob.   ");

	      
	    //ATTACHMENT: 
	    for(II=1;II<=tauT;II++)
	      {
		if(s<N)    // parche para evitar que escriba en los vectores, más alla de s=N
		  {
		    s++;	    
		    nuevo_nodo();
		    //		    printf("s:%d  ",s);	      		   		    
		  }
	      	
		
	      }       //fin del bucle de tau_t
	    
	    
	    
	  }while(s<N);  //////////////////acaba de crecer la red	
	  
	  //intf("\nAcabo de crecer");

	  c_media=c_media+n_coop;	
	  
	  histograma_pk();
	  histograma_pk_coop();   



	  for(i=1; i<=N; i++)            //acumulo
	    {
	      PK_tot[i]+=PK[i];	
	      PK_coop_tot[i]+=PK_coop[i];   
	    }
	  
 
	  if(iter==1)    //solo una muestra
	    {
	      fich4=fopen(file4,"at");
	      for(i=1;i<=N;i++)
		{
		  fprintf(fich4,"%d:  ",i);      
		  for(j=1;j<=N;j++)
		    {
		      if(C[i][j]!=0)
			{
		      fprintf(fich4,"%d  ",C[i][j]);      
			}
		    }
		  fprintf(fich4,"\n");      
		}

	      fprintf(fich4,"\n \n");      
	      fclose (fich4);
	    }
	  
	  /////////////////////////////////////////////////////////
	  ////////////// Ahora empieza la etapa en que la red juego, gana un nodo y pierde otro

	  tpo=1;
	  do{
	    //	    printf("\n\ntpo:%d       ",tpo);

	    for(II=1;II<=tauD;II++) 
	      {
		s=N;   //por si acaso...
		jugada();
		
		//printf("jugada  "  );
		n_coop=0;
		for(i=1;i<=N;i++)
		  {
		    if(e[i]==1)	    
		      n_coop++;	      
		  }
		//printf("n_coop:%d  ",n_coop);	                         	

		
		if(iter<=5)    //solo unas muestras
		  {
		    fich3=fopen(file3,"at");
		    fprintf(fich3,"%d   %d\n",s+tpo,n_coop);      
		    fclose (fich3);
		  }	

		if(n_coop==N  || n_coop==0)
		  {
		    tpo=jugadas_evolucion;   //para que se salte el resto de los pasos temp 
		    break;

		  }		
					
	      }


	   
	    
	    //ATTACHMENT: (dentro de la funcion muerte_nodo(), calculo las prob de attach.)
	    for(II=1;II<=tauT;II++)
	      {	
		
		muerte_nodo();      //para mantener el tamaño cte, luego quito	
		
		
		//Recalculo Probabilidad de attachment: 	    
		for(i=1;i<=N;i++)
		  {P[i]=0.;}
		
		for(i=1;i<=N;i++)
		  {
		    for(j=1;j<=i;j++)  
		      {
			P[i]=P[i]+fit[j]; 
		      }
		  }
		norma=P[N];	    
		
		
		
		
		nuevo_nodo_hueco();   //primero añado, con la P attach recient calculada
		
		
		
	      }  

	    
	  
	    
	    tpo++;
	  }while(tpo<jugadas_evolucion);  /////acaban las jugadas extras 
	                                  //growth_play_dead


	  c_media_fin=c_media_fin+n_coop;	
	  
	  histograma_pk();
	  histograma_pk_coop();   



	  for(i=1; i<=N; i++)            //acumulo
	    {
	      PK_tot_fin[i]+=PK[i];	
	      PK_coop_tot_fin[i]+=PK_coop[i];   
	    }
	  
 
	  if(iter==1)    //solo una muestra
	    {
	      fich4=fopen(file4,"at");
	      for(i=1;i<=N;i++)
		{
		  fprintf(fich4,"%d:  ",i);      
		  for(j=1;j<=N;j++)
		    {
		      if(C[i][j]!=0)
			{
			  fprintf(fich4,"%d  ",C[i][j]);      
			}
		    }
		  fprintf(fich4,"\n ");      
		}
	      
	      fprintf(fich4,"\n \n");      
	      fclose (fich4);
	    }
	  
	  
	  
	}   ///////////////// fin bucle de estadistica
      
      c_media=c_media/(N*Niter);
      c_media_fin=c_media_fin/(N*Niter);
      
      printf("\n\nb:%lf   c_med:%lf   c_med_fin:%lf\n",b, c_media, c_media_fin);
      
      
      
      // Escribo el nivel de cooperacion para este valor de b
      fich1=fopen(file1,"at");
      fprintf(fich1,"%lf   %lf   %lf   %lf   %d   %d\n",b,c_media,c_media_fin,eps,tauT,tauD); 
      fclose (fich1);
      
      
      
      
      
      // Escribo Pk y Pkc, y las normalizo	
      fich2=fopen(file2,"at"); 
      for(i=1; i<N; i++)        
	{	   	 
	  if(PK_tot[i]!=0)    //si una clase de conectividad no esta presente, no escribo nada
	    {
	      fprintf(fich2,"%d   %lf   %lf   %lf   %lf\n",i, PK_tot[i]/(N*Niter),
		 PK_coop_tot[i]/PK_tot[i], PK_tot_fin[i]/(N*Niter),PK_coop_tot_fin[i]/PK_tot_fin[i]);
	    }
	}                  		
      fclose (fich2);
      
      
      
      
      
      
      
      b+=delta_b;    
      
    }   ///////////////// fin bucle de barrer en b
  
  
  
  
  
}




////////////////////////////////////////
///////////////////////////////////////


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
  
  int i, j;

  n_coop=0;

  for(i=0; i<=N; i++)
    {
      k[i]=0;     
      
      for(j=0; j<=m_o; j++)
	{
	  M[i][j] =0;
	}
      
      for(j=0;j<=N;j++)
	{
	  C[i][j]=0;
	}
    }
  
  
  for(j=1; j<=m_o ; j++) 
    k[j]=m_o-1;
  
  
  for(i=1; i<=m_o ; i++)
    {
      counter=0;
      for(j=1; j<=m_o ; j++)
	{
	  if(i!=j)
	    {
	      counter++;
	      C[i][counter]=j;
	    }
	}
    }
  

    for(i=2; i<=m_o ; i++)
    {
	for(j=1; j<=i-1 ; j++)
	{
	  M[i][j] = j;	 
	}
    }
  
    for(i=1;i<=m_o;i++)                 
    { 
	fit[i]=1.-eps;  
    }

    for(i=0; i<=N; i++)
      {
	P[i] = 0;
	P_prov[i]=0;
      }
    

  for(i=1; i<=m_o; i++)                    //prob de los nodos iniciales
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
      e[i]=1;       //cooperador
      n_coop++;
    }
  
  
}




//////////////////////////////////////
//////////////////////////////////////




void nuevo_nodo()    //el indice del nodo nuevo es s.
{
  
  double r,v;
  int i,j,q,g;
  
  r=FRANDOM;   // establezco el carcter del nuevo nodo:
  
  if(r<ro)
    {                
      e[s]=1;       //cooperador
    }


  for(i=0; i<=m; i++)
    {
      unido[i]=0;
    }
  
  for(i=1; i<=N; i++)
    {
      P_prov[i]=P[i];   //guardo la prob en un vector aux para manipularlo
    }  
  
  norma2=norma;
  P[s]=0.;       //excluyo al propio nodo recien llegado
  P_prov[s]=0.;  
  

  for(q=1; q<=m; q++)        //bucle a los m links nuevos que voy a añadir
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
      
      for(j=g+1; j<s; j++)   //recalculo todas las probabilidades
      {
	  P_prov[j]=P_prov[j]-fit[g] ; 
      }
      
      norma2=norma2-fit[g];
      
     
    }   ///////fin del bucle sobre los m links lanzados


  for(j=1;j<=m;j++)   
    {
	g=unido[j];
	k[g]=k[g]+1;
	C[g][k[g]]=s;
	C[s][j]=g;
	M[s][j]=g;
               
    }
  
  k[s]=m;
  fit[s]=(1.-eps);
  P[s]=P[s-1]+fit[s];
  norma=P[s];
 
 
}


///////////////////////////////////////////
/////////////////////////////////////////



void matriz_Conect()   //obtencion de la matriz C a partir de la M
{

  int i,j,z;
    
    for(i=0;i<=N;i++)
      {
	x[i]=0;
	
	for(j=0;j<=N;j++)
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
	
		C[i][x[i]]=z;
		x[z]++;
		C[z][x[z]]=i;	 
	      }
	  }
      }
    
}


//////////////////////////////////////////
////////////////////////////////////


void jugada()  //convenio:  1==cooperador, 0==defector
{
  double r,v;
  int w,t,i,j,y;
    
  
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
	      
	      if(e[i]==1 && e[y]==1)   //C contra C
		{
		  ben[i]+=R;
		  ben[y]+=R;
		  //printf("%d : %d  empate coop \n", i,y);
		}
	      
	      if(e[i]==1 && e[y]==0)   //C contra D
		{
		  ben[i]+=Su;
		  ben[y]+=b;
		  //printf("%d : %d  coop contra defect\n", i,y);
		}
	      
	      if(e[i]==0 && e[y]==0)   //D contra D
		{
		  ben[i]+=Pu;
		  ben[y]+=Pu;
		  
		  //printf("%d : %d  empate defect \n", i,y);
		}
	      
	      if(e[i]==0 && e[y]==1)   //D contra C
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
      e_aux[i]=e[i]; //guardo la estrategia del ultimo paso 
    }                 //  (para el update paralelo)
  
  
  
  
  for(i=1;i<=s;i++)        //comparo estrategias   
    {
      r=FRANDOM;
      r=r*k[i];
      w=(int)r+1;
      
      
      t=C[i][w];    //el vecino elegido al azar
      
      if(k[i]>k[t])
	{k_m=k[i];}
      else
	{k_m=k[t];}
      
     

      if(ben[i]<ben[t])
	{
	  p=(ben[t]-ben[i]);
	  p=p/(b*(double)k_m);
	  
	  v=FRANDOM;
	  if(v<p)
	    {
	      e_aux[i]=e[t];    //IMITACION
	      
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

  int i;

    for(i=0; i<=N; i++)           //inicializo
    {
	PK[i]=0.;	
    }
    
    for(i=1; i<=N; i++)            //recuento
    {
	PK[k[i]]++;	
    }        
    
    
}        



/////////////////////////////////////
///////////////////////////////////////////


void histograma_pk_coop()            //costruccion de P(k) de cooperadores 
{                              //      (pero en su k[i] cuentan coop y defect!)
  
  int i;
  for(i=0; i<=N; i++)           
    {
      PK_coop[i]=0.;	
    }
  
  
  for(i=0; i<=N; i++)
    {
      if(e[i]==1)	    
	PK_coop[k[i]]++;        //sin normalizar   (lo hago al imprimir)
    }   
  
}



/////////////////////////////////////
///////////////////////////////////////////

void muerte_nodo()      //elijo un nodo (más prob cuanto menor payoff) 
{               // y lo elimino junto con todas sus conexiones

  double r,normalizo,normalizo2=0., Pdie[N+1];
  int i,j,x[N+1],vecino, posicion;

  for(i=0;i<=N;i++)
    {     
      x[i]=0;
      Pdie[i]=0;      
      for(j=1;j<=N;j++)      //copio la matriz de conect, para manipularla
	{	 
	  C_aux[i][j]=C[i][j];	  
	}     
    }

  // calculo la norma de la prob.
  normalizo=0.0;
  for(i=1;i<=N;i++)
    {      
      normalizo+= ben[i];    
    } 
  
  //printf("norma:%f\n",normalizo);
  
  
  
  if(normalizo!=0)
    {
      for(i=1;i<=N;i++)
	{   	     
	  Pdie[i]=1.0-ben[i]/normalizo;     
	  //      printf("%d: ben:%f   k%d   Pdie:%f\n",i,ben[i],k[i],Pdie[i]);
	}	 
    } 
  else
    {
      for(i=1;i<=N;i++)
	{   	     
	  Pdie[i]=1.0;
	}
    }
  
  
  
  for(i=1;i<=N;i++)      //acumulada (para que el intervalo que ocupa cada nodo
    {              //sea proporc a sus benef)     
      Pdie[i]=Pdie[i]+Pdie[i-1];    
    }   
  
  normalizo2=Pdie[N];
  
  
  
  r=FRANDOM*normalizo2;     
  for(i=1;i<=N;i++)
    {      
      if(r<Pdie[i])
	{     	  
	  hueco=i;
	  //	  printf("NODO ELIMINADO: %d  (Pdie:%f), r:%f\n",hueco,Pdie[hueco],r);  
	  break;
	}
    }
  
  
  
  
  for(i=1;i<=k[hueco];i++)   //borro a v del vecindario de sus k[v] vecinos
    {     
      vecino=C[hueco][i];
     
      for(j=1;j<=N;j++)     
	{
	  if(C[vecino][j]==hueco)
	    {
	      C_aux[vecino][j]=0;	     
	      posicion=j;
	     
	      break;	     
	    }
	}
      for(j=posicion;j<N;j++)   //como queda un hueco al borrar v, 
	{         //  corro todos los vecinos una posicion hacia delante. 
	  C[vecino][j]=C_aux[vecino][j+1];	  	    
	}
      C[vecino][N]=0;    //porsiaca
      k[vecino]--;    //le resto una conexion al vecino corresp
     
    }

  k[hueco]=0;     //anulo la conenectividad de v

  for(j=1;j<=N;j++)    //borro los vecinos de v
    {
      C_aux[hueco][j]=0;
      C[hueco][j]=0;
    }
  


 
  ben[hueco]=0;
  P[hueco]=0;


  for(i=hueco+1;i<=N;i++)
    {
      P[i]=P[i]-fit[hueco];
    }
  norma=norma-fit[hueco];



}


/////////////////////////////////////////
//////////////////////////////////////
///////////////////////////////////////



void nuevo_nodo_hueco()  //  "hueco" es el indice del nodo nuevo
{
  
  double r,v;
  int i,j,q,g;
  
  r=FRANDOM;   // establezco el carcter del nuevo nodo:
  
  if(r<ro)
    {                
      e[hueco]=1;       //cooperador
    }


  for(i=0; i<=m; i++)
    {
      unido[i]=0;
    }
  
  for(i=1; i<=N; i++)
    {
      P_prov[i]=P[i];   //guardo la prob en un vector aux para manipularlo
    }  
  
  norma2=norma;

  //"hueco" ha quedado excluido al morir
 
  

  for(q=1; q<=m; q++)        //bucle a los m links nuevos que voy a añadir
    {      
      v=FRANDOM*norma2;
     
      for(j=1; j<N; j++)       //recorro la red
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
      
      for(j=g+1; j<N; j++)   //recalculo todas las probabilidades
	{
	  P_prov[j]=P_prov[j]-fit[g] ; 
      }
      
      norma2=norma2-fit[g];
      
     
    }   ///////fin del bucle sobre los m links lanzados


  for(j=1;j<=m;j++)   
    {
	g=unido[j];
	k[g]=k[g]+1;
	C[g][k[g]]=hueco;
	C[hueco][j]=g;
	M[hueco][j]=g;
               
    }
  
  k[hueco]=m;
  fit[hueco]=(1.-eps);
  P[hueco]=P[hueco-1]+fit[hueco]; //asigno un intervalo no nulo al nuevo nodo

  for(i=hueco+1;i<=N;i++)   //y corro todos los siguientes intervalos
    {
      P[i]=P[i-1]+fit[i];
    }
  norma=P[N];     
 
 
}


///////////////////////////////////////////
/////////////////////////////////////////

