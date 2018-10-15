// 17-2-10: Modelo de construccion de redes que combina el crecimiento con el juego, 
// pero ahora, para elegir a kien me uno, miro lo que ganan sus vecinos, de media



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



# define b_min    1.0
# define b_max    3.001
# define delta_b 0.1



# define jugadas_extra   5000       //despues de acabar de crecer la red  y tb momento en el que calculo la media de payoff de coop y defect por clases de conectividad




////////////Generacion de numeros aleatorios//////////
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
///////////////////////////////////////////////



int k[N+1];  //conectividad topológica 
int   M[N+1][m+1], unido[m+1], Cuidado[N+1]; //  unido[] guarda a quién se ha unido, y cuidado[] quien le ha lanzadolinks
double P[N+1], P_prov[N+1];     //P_prov[] es la provisional para manipularla
double PK[N+1], PK_tot[N+1], PK_coop[N+1], PK_coop_tot[N+1], PK_coop_tot_inf[N+1];    //  para la funcion histograma_P(k)
                                


double r,v;                     //para guardar los aleatorios
int  II,  C[N+1][N+1], steps,iter,i,j;
double norma,norma2;
int si[N+1], cont,cont2;



int  e[N+1],e_aux[N+1]; 
double  ben[N+1];
int  n_coop; //numero de coop. instantaneos
int iter,n,ii;

int s;
double fit[N+1],c_media,c_media_inf;       
double b;



void inicia_rand(int semilla);
void nucleo_inicial();
void jugada();
void nuevo_nodo();
void histograma_pk();
void histograma_pk_coop();



char file1[256],file2[256],file4[256],file7[256];

FILE *fich1;      
FILE *fich2;     
FILE *fich4;      
FILE *fich7;     


int main()
{


 printf("\n\nPrograma Growth & Play con Prob.attch-payoff_vecinos  (%d iteraciones)\nN=%d  b_ini=%.2lf b_fin=%.2lf eps=%.2lf  TauT=%d  TauD=%d  %d jugadas extra\n", Niter,N,b_min,b_max,eps,tauT,tauD,jugadas_extra);





 inicia_rand(time(0));


 sprintf(file4,"%dProp-top_%d-%d_b%.2lf-%.2lf_e%.2lf_tpo%d_%diter_neighbors.dat",N,tauT,tauD,b_min,b_max,eps,jugadas_extra,Niter);   //un solo archivo     
    fich4=fopen(file4,"wt");
    fclose(fich4);
    


 b=b_min;



  while(b<=b_max)     //bucle externo para barrer en b
    {
	printf("\nb=%lf\n",b);
	

	sprintf(file7,"%dPk%d-%d-b%.2lf-e%.2lf_%diter_neighbors.dat",N,tauT,tauD,b,eps,Niter);      
	fich7=fopen(file7,"wt");
	fclose(fich7);      


	sprintf(file2,"%dPkc%d-%d-b%.2lf-e%.2lf_tpo%d_%diter_neighbors.dat",N,tauT,tauD,b,eps,jugadas_extra,Niter);      
	fich2=fopen(file2,"wt");
	fclose(fich2);      
	



	sprintf(file1,"%dEvolCoop%d-%d-b%.2lf-e%.2lf_tpo%d_%diter_neighbors.dat",N,tauT,tauD,b,eps,jugadas_extra,Niter);  
	fich1=fopen(file1,"wt");
	fclose(fich1);


	
	c_media=0.;
	c_media_inf=0.;
	for (i=0;i<=N;i++)
	{	    
	    PK_tot[i]=0.;		 
	    PK_coop_tot[i]=0.;	
	    PK_coop_tot_inf[i]=0.;	    
	}


	
	for(iter=1;iter<=Niter;iter++)         //inicio bucle de estadistica
	{	    
	    printf("\niter:%d\n",iter);




	   
	    
	    for(i=1;i<=N;i++)      
	      {
		ben[i]=0;
		e[i]=0;     //TODOS DEFECTORES      
	      }
	    
	    
	    
	    nucleo_inicial();
	    

	  
	    s=m_o;         //tamaño actual de la red: s=N(t)	   	   

	    do{
	      
	      for(II=1;II<=tauD;II++) 
		{

		  jugada();
		 
		  n_coop=0;
		  for(i=1;i<=s;i++)
		    {
		      if(e[i]==1)	    
			n_coop++;	      
		    }
		  


		  if(iter<=5)    //pq si no, el archivo se hace demasiado grande
		    {
		      fich1=fopen(file1,"at");
		      fprintf(fich1,"%d   %d\n",s,n_coop);   		     
		      fclose (fich1);
		      //printf("%d   %d\n",s,n_coop);	
		    }
		  
		} //fin del bucle de tau_d
	      	      	     


	      //ahora no recalculo aki las P attch, pq lo hago antes de añadir cada nodo,
	      // en la funcion nuevo_nodo


	     
	      for(II=1;II<=tauT;II++)
		{

		  if(s<N)    // parche para evitar que escriba en los vectores
		    {      // más alla de s=N
		      s++;	    
		      nuevo_nodo();	      
		      //printf("\n%i  ",s);
		    }	
		  
		}       //fin del bucle de tau_t
	      
	      
	     
	      
	    }while(s<N);    ///////////////acaba de crecer la red
	    
	    c_media+=n_coop;
	    


	    
	    histograma_pk();     
	    histograma_pk_coop();
	   	    
	    for(i=1; i<=N; i++)        
	      {				 //ojo!!! no normalizo aki, sino al fprintf
		PK_tot[i] += PK[i]; 
		PK_coop_tot[i] += PK_coop[i]; 				
	      }		  	
	    



	    for(ii=1;ii<=jugadas_extra;ii++)  //bucle para continuar jugando 
	      {                             //      con la red ya crecida
				
		jugada();
		

		n_coop=0;
		for(i=1;i<=s;i++)
		  {
		    if(e[i]==1)	    
		      n_coop++;	      
		  }
		

		if(iter<=5)    //pq si no, el archivo se hace demasiado grande
		  {
		    fich1=fopen(file1,"at");
		    fprintf(fich1,"%d   %d\n",s+ii,n_coop);      
		    fclose (fich1);
		  }
				

		if(n_coop==0 || n_coop==1)   //stop
		  break;
				
	      }
	    
	    c_media_inf+=n_coop;
	    
	    
 
	   
	    histograma_pk_coop();
	    
	    for(i=1; i<=N; i++)        
	      {				 //ojo!!! no normalizo aki, sino al fprintf
		PK_coop_tot_inf[i] += PK_coop[i]; 				
	      }		  	
	    
	    


	}  ///////////////fin bucle estadistica


	c_media=c_media/(N*Niter);
	c_media_inf=c_media_inf/(N*Niter);


	printf("\nb:%lf   c_med:%lf   c_med_inf:%lf   \n\n",b, c_media,c_media_inf);



	



	fich4=fopen(file4,"at");
	fprintf(fich4,"%lf   %lf   %lf   %lf   %d   %d\n",b, c_media,c_media_inf, eps, tauT, tauD);      
	fclose (fich4);
	





	fich7=fopen(file7,"at");     // Escribo y normalizao la PK    
	for(i=1; i<N; i++)          
	  {	   
	    if(PK_tot[i]!=0)
	      {
		fprintf(fich7,"%d  %lf\n", i, PK_tot[i]/(N*Niter));
	      }
	  }                  		
	fclose (fich7);
	



	
	fich2=fopen(file2,"at");     // Escribo y normalizao la PK          
	for(i=1; i<N; i++)          
	  {	   
	    if(PK_tot[i]!=0)
	      {
		fprintf(fich2,"%d  %lf  %lf\n", i, PK_coop_tot[i]/PK_tot[i], PK_coop_tot_inf[i]/PK_tot[i]);
	      }
	  }                  		
	fclose (fich2);
	
	






	b+=delta_b;    
    } ///////////////fin bucle sobre b


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

  int counter, i,j;

    for(i=0; i<=N; i++)
    {
	k[i]=0;     
	fit[i]=0;       
      
      for(j=0; j<=m_o; j++)
      {M[i][j]=0;}
 
      for(j=0;j<=N;j++)
      {C[i][j]=0;}
    }
  
  
    for(j=1; j<=m_o ; j++) 
	k[j]=m_o-1;

   

    for(i=1; i<=m_o ;i++)
      {
	counter=0;
	for(j=1;j<=m_o;j++)
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
    
    for(i=0; i<=N;i++)
      {
	P[i]=0;
	P_prov[i]=0;
      }
    
    for(i=1; i<=m_o; i++)                    //prob de los nodos iniciales
      {                  //       (==intervalos propoc a la fitness de cada nodo)
	for(j=1; j<=i; j++)
	  {                              
	    P[i]=P[i]+fit[j];       
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



void nuevo_nodo()      // s es el indice del nuevo nodo
{
  int i,q,j,g;
  int vecino,vecino_vecino;
  double r, media,media_media;
  double ganancias[N+1]; //guarda la media de la media de lo que ganan los vecinos de i

  
  r=FRANDOM;// establezco el carcter del nuevo nodo:
  
  if(r<ro)
    {                
      e[s]=1;       //cooperador
    }


  for(i=0;i<=N;i++)
    {
      ganancias[i]=0;
      P[i]=0;
    }

  norma=0.;
  for(i=1;i<s;i++)    //recorro la red
    {
      

      media_media=0.;     // media de la media 
      for(j=1;j<=k[i];j++)  //recorro todos los vecinos del nodo i
	{	 
	  vecino=C[i][j];
	 

	  media=0.;       //media sobre los vecinos del vecino de i
	  for(g=1;g<=k[vecino];g++)
	    {	      
	      vecino_vecino=C[vecino][g];	    

	      media+=ben[vecino_vecino];
	    }
	  media=media/k[vecino];
	  media_media+=media;
	}
      media_media=media_media/k[i];
      norma+=media_media;     //tamaño del intervalo total a repartir entre las prob. attch

      ganancias[i]=media_media; 
      
    }




  //calculo del fitness: AHORA ES UNA FUNCION LINEAL DE LAS GANANCIAS DE MIS VECINOS!!!!
  
  for (i=1;i<s;i++)
    {
      fit[i]=1.0-eps+eps*ganancias[i];     
    }
  
 


  for(i=1;i<s;i++)   //establezco los intervalos de prob. attch de cada nodo
    {
      P[i]=fit[i]+ P[i-1];      
      //printf("gan(%d):%f   P(%d):%f\n",i,ganancias[i],i,P[i]);
    }
  norma=P[s-1];
 

  for(i=0; i<=m; i++)
    {unido[i]=0;}
  
  for(i=1; i<=N; i++)
    {P_prov[i]=P[i];}   //guardo la prob en un vector aux para manipularlo
  
  norma2=norma; //la norma2 es para manipularla mientras lanzo los m links


  P[s]=0.;  
  P_prov[s]=0.;       //anulo la prob de elegir al propio nodo recien llegado
  
  
  for(q=1; q<=m; q++)        //bucle a los m links nuevos que voy a añadir
    {
      
      v=FRANDOM*norma2;
     
      for(j=1; j<s; j++)
      {
	  if(v<P_prov[j])        //elijo el nodo
	  {
	      unido[q]=j;
	     
	      break;
	  }
      }
     

      g=unido[q];
      P_prov[g]=0;    //lo guardo y anulo su prob para no volver a cogerlo
      
      for(j=g+1; j<s; j++)  //desplazo todos los intervalos
	{
	  P_prov[j]=P_prov[j]-fit[g] ; 
	}
      
      norma2=norma2-fit[g];
     
      
    }   ////////fin del bucle sobre los m links lanzados
  
   
  
  for(j=1;j<=m;j++)     //actualizo C[][] y k[] de los nodos implicados
    {
      g=unido[j];
      k[g]++;
      C[g][k[g]]=s;
      C[s][j]=g;
      M[s][j]=g;
      
      
    }
  
  k[s]=m;


  fit[s]=(1.-eps);
  P[s]=P[s-1]+fit[s];     
  norma=P[s];            
  
  
}

//////////////////////////////////////////
////////////////////////////////////

void jugada()   // convenio:  0==defect, 1==coop
{

  int i,j,y,w,t,k_m;
  double r,p,v;
  
  for(i=1;i<=s;i++)        // tamagno del sist = s
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
	      if(e[i]==1 && e[y]==1)
		{
		  ben[i]+=R;
		  ben[y]+=R;		   
		}
	      
	      if(e[i]==1 && e[y]==0)
		{
		  ben[i]+=Su;
		  ben[y]+=b;		   
		}
	      
	      if(e[i]==0 && e[y]==0)
		{
		  ben[i]+=Pu;
		  ben[y]+=Pu;		    		   
		}
	      
	      if(e[i]==0 && e[y]==1)
		{
		  ben[i]+=b;
		  ben[y]+=Su;		    
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
      
      
      
      if(ben[i]<ben[t])
	{
	  p=(ben[t]-ben[i]);
	  p=p/(b*(double)k_m);
	  
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
  
  
  
}

/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////




void histograma_pk()            //costruccion de P(k)

{
  int i;
  
  for(i=0; i <=N; i++)           //inicializo
    {
      PK[i]=0.;	
    }
  
  for(i=1; i<=s; i++)            //recuento
    {
      PK[k[i]]++;	
    }         
}        




//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

  

void histograma_pk_coop()   //costruccion de P(k) de cooperadores 
{              //       (pero su k[i] es la topologica: cuentan coop y defect!)
    
  int i;
  for(i=0; i<=N; i++)           
    {
      PK_coop[i]=0.;      
    }
  
  
  for(i=0; i<=s; i++)
    {
      if(e[i]==1)	    
	PK_coop[k[i]]++;    
    }    
  
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

  
