//8-2-2010: Programa para simular una variacion del modelo EPA, que combinaba
// crecimiento y dinamica del DP, pero ahora ademas, añadimos muerte de los 
// nodos como funcion de su payoff, de modo que los más ricos, viven mas

//(pero) convenio de estrategias:  1==cooperador, 0==defector

//basado en el programa "growth_play_dead_payoff_power-law.c

//25-2-10:

// natalidad constante, pero mortalidad variable: depende de la exp de 
// la vida del nodo, asi como de su fitness, y se aplica a todos, uno por uno

// primero juego+crezco hasta s=50, y luego empiezan a morir

// probabilidad  attchment proporcional a la fitness


//ahora no hay fitness[i], sino directamente ben[i] -mas una baseline agnadida
//  a todos por igual.



# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

# define N 10000 //tamaño MAXIMO  de la red que voy a permitir
# define m_o  2      //nodos inicialmente unidos
# define m   2      //links añadidos con cada nodo nuevo



# define N_links 50000 //estimacion -por lo alto- para el tamaño de la matriz D[][]
# define K_MAX  5000   //estimacion por lo alto para para el programa en su caso


# define  beta    0.001     //si grande, los nodos viven poco, si peq, viven mucho



# define Niter 50   //estadistica


# define  tauT 1       //(si tauT=10 y tauD=1, pongo 10 nodos y juego una vez)
# define  tauD 1
 

# define  ro 0.5    // prob inicial de ser cooperador


# define  R 1.01          //notar que ahora hay una peq. fitness baseline
# define  Pu 0.01      //  para que si all-D, haya prob de attach !=0
# define  Su 0.01



# define b_min    1.0
# define b_max   3.01
# define delta_b 0.2



# define size_no_dead  50       // tamaño hasta el cual el mecanismo de muerte 
                               //  no se activa


# define jugadas_evolucion   2000    //tiempo que dejo al sistema crecer+jugar+morir
   //(pero tpo=0 es cuando empieza la simu, con  m_o). el tiempo pasa en cada JUGADA

//////Generacion de numeros aleatorios/////////////////
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
/////////////////////////////////////////////////////////



int k[N+1]; //conectividad topológica y de Preferential Attachment (PA)
double norma,norma2;
int   M[N+1][m+1], unido[m+1], Cuidado[N+1]; //  unido[] guarda a quién se ha unido, y cuidado[] quien le ha lanzadolinks
double P[N+1], P_prov[N+1];     //P_prov[] es la provisional para manipularla
int  C[N+1][N+1], x[N+1];
int  C_aux[N+1][N+1];
double PK[N+1],  PK_coop[N+1],  NK_tot[N+1],  PK_tot[N+1], PK_coop_tot[N+1];




int  e[N+1],e_aux[N+1], k_m; 
double  ben[N+1],p,r;
int  n_coop; //numero de coop. instantaneos
int iter,counter,n,ii, D[N_links+1][3], D_aux[N_links+1][3],num_links;

int s,II,i,j, tpo,size, hueco;
double tpo_vida[N+1];
double Pdie[N+1];


double  c_media_fin,size_med;       
double b, aux;
double NClusters_Med,NClusters,GCmed,GComp;


void inicia_rand(int semilla);
void nucleo_inicial();
void jugada();
void nuevo_nodo();
void muerte_nodo();
void histograma_pk();
void histograma_pk_coop();
void escribir_matriz_D();
void borrar_matriz_D(int, int);
void ver_clusters();

char file1[256], file2[256], file3[256], file4[256], file5[256], file6[256];

FILE *fich1;     
FILE *fich2;   
FILE *fich3;   
FILE *fich4;   
FILE *fich5;   
FILE *fich6;   



int main()
{
  
  inicia_rand(time(0));
  







  printf("\n\nPrograma Growth & Play & Dead (modelo payoff exp)  (%d iter)\nN_max=%d  b_ini=%lf b_fin=%lf   TauT=%d  TauD=%d\n", Niter,N,b_min,b_max,tauT,tauD);
  
  
  
  
  
  sprintf(file1,"C_y_N_vs_b_%d-%d-b_ini%.2lf_beta%.8lf_%diter.dat",tauT,tauD,b_min,beta,Niter);  
  fich1=fopen(file1,"wt");
  fclose(fich1);
  
  
  

      sprintf(file6,"Clsters_%d-%d-b%.2lf-%.2lf_beta%.10lf_%diter.dat",tauT,tauD,b_min,b_max,beta,Niter);        
      fich6=fopen(file6,"wt");
      fclose(fich6);   
  
  
  
  b=b_min;    
  while(b<=b_max)     //bucle externo para barrer en b
    {
      printf("\nb=%lf\n",b);
      
      
      
      
      sprintf(file2,"Pk_Pkc_%d-%d-b%.2lf_beta%.8lf_%diter.dat",tauT,tauD,b,beta,Niter);        
      fich2=fopen(file2,"wt");
      fclose(fich2);   
      
      

      
      sprintf(file3,"EvolCoopSize_%d-%d-b%.2lf_beta%.8lf_%diter.dat",tauT,tauD,b,beta,Niter);      
      fich3=fopen(file3,"wt");
      fclose(fich3);
      
      
      
      /*      sprintf(file4,"Conex_%d-%d-b%.2lf_beta%.8lf_%diter.dat",tauT,tauD,b,beta,Niter);        
      fich4=fopen(file4,"wt");
      fclose(fich4);   */
      


      for(i=1; i<=N; i++)       
	{	  
	  PK_tot[i]=0;
	  NK_tot[i]=0;    //esto acumulara sin normalizar,y lo usare para normalizar PK_coop
	  PK_coop_tot[i]=0;
	}
      
      
     
      size_med=0.;   
      c_media_fin=0.;

      GCmed=0.;
      NClusters_Med=0.;



      for(iter=1;iter<=Niter;iter++)         //inicio bucle de estadistica
	{
	  printf("\n\niter:%d\n",iter);
 	  	  

	  /*  if(iter<=2)   //como muestra
	    {
	      sprintf(file5,"matrizD_%d-%d-b%.2lf_beta%.8lf_iter%d.NET",tauT,tauD,b,beta,iter);
	      fich5=fopen(file5,"wt");	      
	      fclose(fich5);   
	      }*/
	 	 


	  for(i=1;i<=N;i++)
	    {
	      tpo_vida[i]=0;
	    }
	  

	  for(i=1;i<=N;i++)      
	    {
	      ben[i]=0;
	      e[i]=0;     //TODOS DEFECTORES      
	    }
	  
	  
	  
	  size=m_o;	 //tamaño real del sist (descontando huecos)
	  s=m_o;         // indice max de un nodo en la red (sin descontar huecos) 
	  



	  nucleo_inicial();        //construyo los m_o primeros nodos, todos con todos	 	  
	 


	  tpo=1;
	  do{  /////bucle de crecim. hasta los size_no_dead primeros nodos (sin muerte)
	    
	    for(II=1;II<=tauD;II++) 
	      {
		jugada();
		
		n_coop=0;
		for(i=1;i<=s;i++)
		  {
		    if(e[i]==1)	    
		      n_coop++;	      
		  }
	

		aux=n_coop;
		aux=aux/size;
		
		if(iter<=5)    //como muestra
		  {
		    fich3=fopen(file3,"at");
		    fprintf(fich3,"%d   %d   %d   %lf   %d\n",tpo,size,n_coop,aux,s);      
		    fclose (fich3);
		  }

		
		//		printf("\n\ntpo:%d, s:%d, size:%d, n_coop:%d\n",tpo,s, size, n_coop);

		for(i=1;i<=s;i++)
		  {
		    tpo_vida[i]++;  //cuando juegan envejecen, 
		  }		   //  no cuando añado un nodo nuevo		



		tpo++;
	      } /////fin del bucle de tau_d


	    
	    for(i=0;i<=N;i++)    /////////Recalculo Prob. attachment: 	    
	      {P[i]=0.;}
	    
	    for(i=1;i<=s;i++)
	      {
		P[i]=P[i-1]+ben[i];
	      }
	    norma=P[s];   
	   

	      
	   
	    for(II=1;II<=tauT;II++)
	      {
		if(s<N)    // parche 
		  {
		    s++;               // es el maximo indice posible de un nodo del sistema	    
		    size++;  //es el tamaño real de la red (descontando huecos)
		    nuevo_nodo();		   
		    //printf("s:%d  size:%d",s,size)		   
		  }
	      			
	      }     /////fin del bucle de tau_t
	    
	    
	   
	  }while(s<size_no_dead);     //// fin de la etapa de crecimiento+juego (sin muerte)
                      

           ///////////////////////////////////////////////////             
	    /////////////////////////////////////////////////
             ///////////////////////////////////////////////
	  	  	 		   	  	  
	   printf("\n\n empiezan a morir\n");

	  do{        ////// Ahora empieza la etapa en que la red juega, crece y muere
	  

	    for(II=1;II<=tauD;II++) 
	      {	
		jugada();		
	

		n_coop=0;
		for(i=1;i<=s;i++)
		  {
		    if(e[i]==1)	    
		      n_coop++;	      
		  }
	                    	
		//		printf("\n\ntpo:%d, s:%d, size:%d, n_coop:%d\n",tpo,s, size, n_coop);
		
		aux=n_coop;
		aux=aux/size;
		if(iter<=5)    //solo unas muestras
		  {
		    fich3=fopen(file3,"at");
		    fprintf(fich3,"%d   %d   %d   %lf   %d\n",tpo,size,n_coop,aux,s); 
		    fclose (fich3);
		  }
		
		
		for(i=1;i<=s;i++)
		  {
		    tpo_vida[i]++;  //cuando juegan, envejecen, 
		  }		//no cuando añado un nodod nuevo
		



		tpo++;			
	      }      /////// fin del bucle de tau_D


  
	    for(i=0;i<=N;i++)    /////////Recalculo Prob. attachment: 	 (de acuerdo a la ultima jugada)   
	      {P[i]=0.;}
	    
	    for(i=1;i<=s;i++)
	      {
		P[i]=P[i-1]+ben[i];

		if(e[i]!=2)   //pq si el nodo esta muerto, P[i]=0
		  {
		    norma=P[s];   
		  }
	      }
	  
	   



	    for(i=1;i<=N;i++)
	      {
		Pdie[i]=0.;   // Prob. morir		
	      }
	    



	    for(i=1;i<=s;i++)   //recorro todos los nodos, para ver a cuáles mato (puedo varios en una ronda!)
	      {
		
		if(e[i]!=2)   //si el nodo no esta muerto ya:
		  {
		    Pdie[i]=1.0-exp(-beta*tpo_vida[i]/ben[i]);
		   
		    //printf("s:%d, Pdie[%d]:%f\n",s,i,Pdie[i]);

		    r=FRANDOM;   
		    if(r<Pdie[i])     //lo mato
		      {
			hueco=i;    
			//	printf("    mato a:  %d,  ",i);
			muerte_nodo();      //borro al nodo y todas sus conexiones, y          
                        	           // modifico Cij y ki    
 
			size--;      
			
		      }		     
		  }
	
		if(size<=0)
		  {		   
		    printf("\n\nPoblacion extinta!");


		    if(iter<=5)    
		      {
			fich3=fopen(file3,"at");
			fprintf(fich3,"%d   %d   %d   %lf   %d\n",tpo,size,n_coop,aux,s);      
			fclose (fich3);
		      }
		    

		    tpo=jugadas_evolucion+10;
		    break;       //salgo del bucle de matar nodos		   
		  }


		for(j=hueco+1;j<=s;j++) ///////Recalculo Prob. attach. de todos los siguientes nodos al muerto
		  {
		    P[j]=P[j]-ben[hueco];
		    
		    if(e[j]!=2)   //la  norma sera la Pdie del ultimo nodo no muerto
		      {
			norma=P[j];	   
		      }
		  }		  						
		
	      }	    /////////// fin del bucle de matar nodos
	    
	    


	                  
	    for(II=1;II<=tauT;II++)
	      {	 
		if(s<N)  //parche
		  {    			 
		    s++;       //indice máx del ultimo nodo del sist
		    size++;      //tamagno real
		    
		    nuevo_nodo();       //dejo los huecos libres y voy aumentando s
		  }
		
	      }       /////// fin del bucle de tau_D

	    
	  	    	  

	    if(s>=N)   //por si el sist alcanza el tamaño max que le he dado a las matrices etc
	      {
		printf("\n\ns=N");

		
		
		if(iter<=5)    //como muestra
		  {
		    fich3=fopen(file3,"at");
		    fprintf(fich3,"%d   %d   %d   %lf   %d\n",tpo,size,n_coop,aux,s);      
		    fclose (fich3);
		  }
		
		tpo=jugadas_evolucion+10;
	      }

	  }while(tpo<=jugadas_evolucion);  /////acaban las jugadas evolucion del sist
	                                 
	  ver_clusters();


	  size_med=size_med+size;
	  c_media_fin=c_media_fin+n_coop;	
	  
	  histograma_pk();  // normalizo dentro de la funcion, por su size alcanzado correspondiente!!
	  histograma_pk_coop();    



	  for(i=1; i<=s; i++)            //acumulo
	    {
	      PK_tot[i]+=PK[i];	        //obtengo NK_tot[i] en la func. histogr()
	      PK_coop_tot[i]+=PK_coop[i];   
	    }
	  
 
	  /*	  if(iter==1)    //solo una muestra
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
	      }*/
	  	  	


	  if(iter<=2)   //como muestra
	    {
	      escribir_matriz_D();
	    }
	  


	  if(iter<=5)    // (para separar iteraciones en el archivo)
	    {
	      fich3=fopen(file3,"at");
	      fprintf(fich3,"\n"); 
	      fclose (fich3);
	    }
	  
	  
	}   ///////////////////////////// fin bucle de estadistica
      
     
      size_med=size_med/Niter;
      c_media_fin=c_media_fin/(size_med*Niter);  // ojo!! ahora no es un tamagno cte!!!
         
      GCmed=GCmed/Niter;
      NClusters_Med=NClusters_Med/Niter;



      printf("  b:%f  c_med:%f  size_med:%f   GCmed:%f  NClusters_med:%f",b,c_media_fin, size_med, GCmed, NClusters_Med);
  

      
      
      // Escribo el nivel de cooperacion para este valor de b
      fich1=fopen(file1,"at");
      fprintf(fich1,"%lf   %lf   %lf   %d   %d\n",b,c_media_fin,size_med,tauT,tauD); 
      fclose (fich1);
      
      
      fich6=fopen(file6,"at");
      fprintf(fich6,"%lf   %lf   %lf   %lf   %d   %d\n",b,size_med,GCmed,NClusters_Med,tauT,tauD); 
      fclose (fich6);
      
      
      
      
      
      // Escribo Pk y Pkc, y las normalizo	
      fich2=fopen(file2,"at"); 
      for(i=1; i<N; i++)        
	{	   	 
	  if(PK_tot[i]!=0)    //si una clase de conectividad no esta presente, no escribo nada
	    {
	      fprintf(fich2,"%d   %lf   %lf\n",i, PK_tot[i]/Niter,PK_coop_tot[i]/NK_tot[i]);
	    }
	}                  		
      fclose (fich2);
      
      
      
      
      
      
      
      b+=delta_b;    
      
    }   ///////////////// fin bucle de barrer en b
  
  
  
    exit(0); 
 
  
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

  num_links=1;    //indice para la matriz D[][]

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

	      D[num_links][1]=i;
	      D[num_links][2]=j;
	      num_links++;
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
  
    

    for(i=0; i<=N; i++)
      {
	P[i] = 0;
	P_prov[i]=0;
      }
    

  for(i=1; i<=m_o; i++)                    //prob de los nodos iniciales
  { 
    ben[i]=0.01;   //como benef. inic. les doy el baseline
      for(j=1; j<=i; j++)
      {                              
	P[i] = P[i] + ben[j];       
      }
      
      P_prov[i]=P[i];
  }
  
  norma=P[m_o];
  norma2=norma;



 // establezco como cooperadores a los  nodos iniciales:  (antes de llamar a esta func,. e[i]=0)
  for(i=1;i<=m_o;i++)                //  para todos los nodos
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

	      D[num_links][1]=s;    //almaceno tb el enlace en la matriz de pares de vecinos
	      D[num_links][2]=j;
	      num_links++;

	      break;
	  }
      }
      
      g=unido[q];
      P_prov[g]=0;    //lo guardo y anulo su prob para no volver a cogerlo
      
      for(j=g+1; j<s; j++)   //recalculo todas las probabilidades del elegido en adelante
      {
	  P_prov[j]=P_prov[j]-ben[g] ; 
      }
      
      norma2=norma2-ben[g];
      
     
    }   ///////fin del bucle sobre los m links lanzados

  // printf("nuevo:%d  vecinos:",s);
  for(j=1;j<=m;j++)   
    {
	g=unido[j];
	k[g]++;
	C[g][k[g]]=s;
	C[s][j]=g;
	M[s][j]=g;   

	//	printf("%d  ",M[s][j]);


	if(k[g]>=K_MAX)
	  {
	    printf("\n\nSobrepasado K_MAX(%d) por nodo %d\n",K_MAX,g );
	    exit(0);
	  }
            
    }
  
  // printf("\n");
  //	getchar();
  k[s]=m;
  ben[s]=0.01;    //pq en cuanto entre en el juego, ganara como minimo 0.01
  P[s]=P[s-1]+ben[s];
  norma=P[s];
 
 
}


///////////////////////////////////////////
/////////////////////////////////////////





void jugada()  //convenio:  1==cooperador, 0==defector,   2==muerto
{
  double r,v;
  int w,t,i,j,y;
    
  
  for(i=1;i<=s;i++)         
    {
      ben[i]=0;
    }
  
  for(i=1;i<=s;i++) /////////////// cada nodo juega con todos sus vecinos
    {
      if(e[i]!=2)  //si el nodo no esta muerto
	{
	  for(j=1;j<=k[i];j++)   // aunque porsiaca, SI HAY UN HUECO, SU k[] y C[][] estan a cero
	    {
	      y=C[i][j];            
	      if(y>i)        //para no repetir vecinos en una partida
		{  
		  
		  if(e[i]==1 && e[y]==1)   //C contra C
		    {
		      ben[i]+=R;
		      ben[y]+=R;		      
		    }
		  
		  if(e[i]==1 && e[y]==0)   //C contra D
		    {
		      ben[i]+=Su;
		      ben[y]+=b;		      
		    }
		  
		  if(e[i]==0 && e[y]==0)   //D contra D
		    {
		      ben[i]+=Pu;
		      ben[y]+=Pu;		      		     
		    }
		  
		  if(e[i]==0 && e[y]==1)   //D contra C
		    {
		      ben[i]+=b;
		      ben[y]+=Su;		     
		    }
		}
	      
	    }
	}    //condicional de estar vivo
    }         ////////fin de la partida (bucle a s, todos los nodos de la red)
  
  

  for (i=1;i<=s;i++)        //guardo la estrategia del ultimo paso 
    {                                         
      e_aux[i]=e[i];  //  (para el update paralelo)          
    }                 
  
  
  
  
  for(i=1;i<=s;i++)     //comparo estrategias   
    {
      if(e[i]!=2)    //si el nodo no esta muerto
	{	  
	  if(k[i]>0) //si el nodo no se ha quedado aislado
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
	}
      
      // si el nodo esta aislado, sigue a su pedal
      
    }
  
  
  for (i=1;i<=s;i++)      //actualizacion de estrategias en paralelo
    {	
      e[i]=e_aux[i];      
    }
  
  
  
  
}



/////////////////////////////////////
///////////////////////////////////////////


void histograma_pk()            //costruccion de P(k)
{

  int i;
  double aux=0;

    for(i=0; i<=N; i++)           //inicializo
    {
	PK[i]=0.;	
    }
    
    for(i=1; i<=s; i++)            //recuento
      {
	if(e[i]!=2)     //si el nodo no esta muerto
	  {
	    PK[k[i]]++;	
	  }
      }        
    

    for(i=1; i<=s; i++)            //recuento
      {
	NK_tot[i]+=PK[i];    //sin normalizar, lo usare para la norma de PK_coop
      }
    

    aux=size;
    for(i=1; i<=s; i++)            //normalizo
      {	
	PK[i]=PK[i]/aux;		
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

void muerte_nodo()       // "hueco" es el índice del nodo a eliminar (elegido dentro del main)
{                
  int i,j,x[N+1],vecino, posicion;

  for(i=0;i<=N;i++)
    {     
      x[i]=0;  


      for(j=1;j<=N;j++)      //copio la matriz de conect, para manipularla
	{	 
	  C_aux[i][j]=C[i][j];	  
	}     
    }

 

  //el "hueco" lo he elegido en la funcion principal
  
  for(i=1;i<=k[hueco];i++)   //borro a v del vecindario de sus k[v] vecinos
    {     
      vecino=C[hueco][i];
     
      for(j=1;j<=N;j++)     
	{
	  if(C[vecino][j]==hueco)
	    {
	      C_aux[vecino][j]=0;	     
	      posicion=j;
	     
	      borrar_matriz_D(vecino,hueco);

	      break;	     
	    }
	}
      for(j=posicion;j<N;j++)   //como queda un hueco al borrar v, 
	{                     //  corro todos los vecinos una posicion hacia delante. 
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
  Pdie[hueco]=0;

  e[hueco]=2;     //ya no es C ni D, sino que está muerto



}


/////////////////////////////////////////
//////////////////////////////////////
///////////////////////////////////////


 
void borrar_matriz_D(int vecino, int hueco)
{
  int i;


  for(i=1;i<=N_links;i++)
    {
      if((D[i][1]==hueco && D[i][2]==vecino) || (D[i][2]==hueco && D[i][1]==vecino)  )
	{
	  D[i][1]=0;
	  D[i][2]=0;	  
	}
    }

}



/////////////////////////////////////////
//////////////////////////////////////
///////////////////////////////////////


 
void escribir_matriz_D()
{
  int i,j,D_aux[N_links][3],relabel[N+1], marca,contador_Daux, vecino;



for(i=1;i<=N_links;i++)
    {
      D_aux[i][1]=0;     //necesito quitar los nodos muertos
      D_aux[i][2]=0; // y renombrar a los vivos (para poder pintar en pajek)
    }



/*for(i=1;i<=s;i++)
   {
     if(k[i]>0)
       {
	 printf("nodo:%d (k:%d) vecinos:%d, %d\n",i,k[i],M[i][1],M[i][2]);
       }
   }

   getchar();*/

  marca=1;
  for(i=1;i<=s;i++)
    {
      if(k[i]>0)   //si el nodo i esta vivo, tiene relabel[i]!=0 y es su nuevo
	{      // indice, descontando ya los muertos
	  relabel[i]=marca;
	  marca++;
	}
      else
	{
	  relabel[i]=0;
	}     
    }


//la matriz M[i][m]  guarda los m vecinos a los que se unio i qdo llego:  
  contador_Daux=1;
  for(i=1;i<=s;i++) //recorro todos los nodos (vivos y muertos)
    {
      if(relabel[i]!=0)
	{

	  // printf("\n%d es %d; (unido a %d, %d que son %d, %d)",i, relabel[i],M[i][1],M[i][2], relabel[M[i][1]], relabel[M[i][2]]);	
	  
	  for(j=1;j<=m;j++)   //recorro los vecinos a los que el se ha unido
	    {
	      vecino=M[i][j];


	      
	      if(relabel[vecino]!=0)     
		{		
		  
		  //printf("vecino:%d   (es %d)",vecino,relabel[vecino]);
		  //getchar();

		  D_aux[contador_Daux][1]=relabel[i];
		  D_aux[contador_Daux][2]=relabel[vecino];

		  // printf("%d-%d\n",D_aux[contador_Daux][1],D_aux[contador_Daux][2]);
		  //  getchar();

		  contador_Daux++;
		
		}
	    }
	}

    }


  /*  fich5=fopen(file5,"at");
  
  fprintf(fich5,"*Vertices      %d\n*Edges\n",size);      
  for(i=1;i<=N_links;i++)
    {
      if(D_aux[i][1] !=0)
	{
	  fprintf(fich5,"%d   %d \n",D_aux[i][1],D_aux[i][2]);      
	  
	}
    }
    fclose (fich5);*/
  

}







//////////////////////////////////////////////
////////////////////////////////////////////
/////////////////////////////////////////



void ver_clusters()
{
  int i, ii,j,label[N+1],cluster_actual[N+1],kk,gg,g;
  double max;

  
    
    for(i=1;i<=N;i++)
      {
	label[i]=0;
      }   //Indica que se ha mirado el nodo con un 1 y que no con un 0


    GComp=0.;
    NClusters=0;

   

  
    for(i=1;i<=s;i++)  //recorro la red 
      {
	
	if(e[i]!=2)  //solo nodos vivos
	  {
	   	kk=0;

		for(j=1;j<=s;j++)   //guardara el cluster que estoy mirando actualmente
   		  {     
		    cluster_actual[j]=0;
		  }
 
		if(label[i]==0)  //si no he mirado el nodo i antes
		  {

		    kk=1;
		    label[i]=1;     //lo marco como mirado
		    cluster_actual[kk]=i;   //y lo incluyo en el cluster act.

		    for(j=1;j<=kk;j++)   //bucle de indice sup. dinámico
		      {
			g=cluster_actual[j];
			for(ii=1;ii<=k[g];ii++)
			  {
			    gg=C[g][ii];
			    if(label[gg]==0)
			      {
				label[gg]=1;
				kk++;
				cluster_actual[kk]=gg;  //entra a formar parte
			      }
			  }

		      }

		    //kk++;     //ESTO QUE HACE?????????
		  }

		if(kk>1)        //si ha encontrado un  cluster
		  {
		    max=kk;		    
		    NClusters++;
		    
		    if (max>GComp)  //guardo el tamaño de la componente gigante
		      GComp=max;
	  	    
		  }	
	
		
	  }
	
      }      //acabo de recorrer la red

    
    //    printf("\nb:%.1lf  iter:%d  GComp:%lf  Nclusters:%lf\n\n",b,iter,GComp,NClusters);
    //    getchar();
    
    NClusters_Med=NClusters_Med+NClusters;    //acumulo
    GCmed=GCmed+GComp;
    

    



}
