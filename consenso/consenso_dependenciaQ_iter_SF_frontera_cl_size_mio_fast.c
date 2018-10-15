// 16-02-09: Programa para implementar Global consensus en networks
//18-2-09: añadida rutina para calcular la GC de links en consenso (Smax).
//20-2-09: añadida rutina para calcular la Sf_max  (componente gigante nodos que comparten el valor del rasgo f) ojoo! ahora el convenio es distinto (para la rutina S_max)
//25-2-09: uso de una matriz D auxiliar, de links no bloqueados, para acelerar
//   el proceso de eleccion aleatoria de links cerca del final de la simu.

//10-3-09: crer archivos con todos los datos, sin hacer medias (estas iran en un programilla de análisis)
//19-3-09: implemento la restriccion del uso del mecanismo de D_no_blocked solo a partir de que el 80 por ciento de los links esten ya bloqueados
//4-5-09: estudio la dependencia de <S> con q y con f (eliminado av_path length)
//4-5-09: añado bucle iteraciones (no las medias)
//4-5-09: cambio la forma de acabar una simu-->> qdo el número de links bloqueados sea igual al máximo posible


//20-5-09: estudio de los nodos frontera en el estado final. su pk.
//21-5-09: estudio del tamaño de los clusters al final de cada iter. histogramas.

//1-5-09: me cai de la parra: va muchisimo mas rápido si solo calculo n_links_blocked y S_max cada multiplo veces!!!!!!!!


# include <stdio.h>
# include <stdlib.h>
# include <math.h>


# define Niter 1

# define N  1000   //tamaño de la red


# define m_o  3      //nodos inicialmente unidos
# define m   3      //nodos nuevos añadidos a cada paso de tiempo (RECORDAR:  <k>=2m)

# define alfa 0.0     //Nota:  para hacer simulaciones sobre redes random, mejor usar la rutina SF.c, que sigue estrictamente el procedimiento de E-R






# define multiplo 10000     //cada cuantos pasos de tiempo calculo maginitudes y escribo archivos




# define Q  600    // numero de valores distintos que puede tomar un rasgo concreto (es igual para todos los rasgos)


# define f 10       // numero de rasgos de un nodo







///////////////////   Generacion de numeros aleatorios
 
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];

//////////////////////////////////////////





char file[256],file1[256],file2[256],file3[256],file4[256],file5[256],file6[256],file7[256];

 
int k[N+1],k_PA[N+1],A[N+1];                    //conectividad topológica y de Preferential Attachment (PA)
int   M[N+1][m+1], unido[m+1], Cuidado[N+1];         //  unido[] guarda a quién se ha unido i , y cuidado[] quien le ha lanzado ya un link
double P[N+1], P_prov[N+1];        //P_prov[] es la provisional para manipularla
double PK[N+1], PK_tot[N+1];        //  para la distribucion P(k)
double r,v;                                  //para guardar los aleatorios
int  i, j, jj,  w,g,gg, d, q, x[N+1],y,z, C[N+1][N+1],n, ValorA, steps,s,flat;
int norma, norma2, norma_aleat,tipo[m_o+1];
double dnorma_aleat;
int si[N+1], cont,cont2,iter;


            


/////// variables para el consenso

int v_i[N+1][f+1]; //matriz que guarda los valores de los f rasgos de los N nodos
double S_ij;
int D[N*m+1][m+1];	//matriz de pares de links 
int n_links_blocked,n_links_S1,tpo,S_max,Sf_max[f+1];

int NN[N+1],Phase[N+1][N+1],semilla[N+1];
int rasgo,Simil[N*m+1],n_links, contador,ii;
double GComp;    //guarda la GC de consenso (total o de un rasgo)
double Sf_max_media,aux;
int prob_D[N*m+1],chosen,cambio, imitador;  //guarda las prob de cada lin
int Nclusters,NCL[f+1],NCLUSTERS;   // numero de clusters de consenso, de cada rasgo,su media y aux para la rutina
double NCL_med;
int semilla_max_cluster;
double dist_tot,dist_tot_cons,dist[f+1],D_med, aux_links,aux_N;
int flag=0,accion,N_links_max,sizes[N+1];
double Pk_frontera_1_tot[N+1],Pk_frontera_2_tot[N+1],histogr_sizes_tot[N+1];



void inicia_rand(int semilla);
void construir_red();
void matriz_D();

int time();


void contar_blocked_links();
void contar_links_S1();
void marcar_link_consenso();
void marcar_link_consenso_rasgo();
void find_Smax();
void consenso_link();
void ver_nodos_frontera();



FILE *fich,*fich1,*fich2,*fich3,*fich4,*fich5,*fich6,*fich7;





int main()
{      
    
    inicia_rand(time(0));
 

//Guarda la Pk
    sprintf(file3,"Pk_N%d_Q%d_f%d_SF_%diter.dat",N,Q,f,Niter);
    fich3=fopen(file3,"wt");
    fclose(fich3);
    

    sprintf(file4,"Pk_frontera_N%d_Q%d_f%d_SF_%diter.dat",N,Q,f,Niter);
    fich4=fopen(file4,"wt");
    fclose(fich4);
    


    
//Guarda el tamaño de todos los cluters presentes al final de cada iter
	sprintf(file5,"Cluster_Size_N%d_Q%d_f%d_SF%diter.dat",N,Q,f,Niter);
	fich5=fopen(file5,"wt");
	fclose(fich5);


//Guarda el tamaño de todos los cluters presentes al final de cada iter
    sprintf(file7,"Histogr_Cl_Size_N%d_Q%d_f%d_SF_%diter.dat",N,Q,f,Niter);
    fich7=fopen(file7,"wt");
    fclose(fich7);
   
    
    printf("\n\n\nPrograma Global consensus on SF Networks  N=%d  Q=%d  f=%d (%d iter)\n\n",N,Q,f,Niter);  
    

    N_links_max=0;           //l= (m_o-1)+(m_o-2)+...+1  +  m(N-m_o)
    j=m_o-1;
    for(i=j;i>=1;i--)
    {
	N_links_max+=i;
    }
    N_links_max+=m*(N-m_o);

    printf("N_max_links: %d\n",N_links_max);  




	for(i = 0; i <= N ; i++)
	{
	    PK_tot[i] = 0.0;
	    Pk_frontera_1_tot[i]=0.0;
	    Pk_frontera_2_tot[i]=0.0;
	    histogr_sizes_tot[i]=0;
	}	       


    for(iter=1;iter<=Niter;iter++)
    {
	
	S_max=0;
	flag=0;
	
	printf("\niter%d\n",iter);  
//	getchar();
	
//Guarda S_max, <Sf_max> y Sf_max[rasgo] en cada paso de tiempo
	sprintf(file,"Consenso_N%d_Q%d_f%d_alfa%.1lf_SF_iter%d.dat",N,Q,f,alfa,iter);
	fich=fopen(file,"wt");
	fclose(fich);
	
//Guarda el numero del clusters de consenso y de cada rasgo  a cada paso de tiempo
	sprintf(file1,"N_Clusters_N%d_Q%d_f%d_alfa%.1lf_SF_iter%d.dat",N,Q,f,alfa,iter);
	fich1=fopen(file1,"wt");
	fclose(fich1);
    

  //Guarda los rasgos culturales todos los cluters presentes al final de cada iter
	sprintf(file6,"Cluster_Rasgos_N%d_Q%d_f%d_SF_iter%d.dat",N,Q,f,iter);
	fich6=fopen(file6,"wt");
	fclose(fich6);  




	
	
	construir_red();
	
	matriz_D();        //construyo la matriz de los links entre pares de nodos
	
	
	
	for(i=1;i<=N;i++)       //inicializo los valores de los rasgos de los nodos	
	{
	    for(j=1;j<=f;j++)
	    { 
		v=FRANDOM;
		v_i[i][j]=v*Q;    //da valores entre 0 y Q-1 a cada rasgo	   	    
	    }	
	}
	
	for(j=0;j<=N*m;j++)  //inicializo a cero las prob de los links  OJO, LAS POSICIONES DE LOS LINKS VAN DE 1 A n_links
	{
	    prob_D[j]=0;
	    
	}
	
	
	tpo=0;
	
	contador=0;
	
	accion=0;   //flag para indicar que hay que usar las rutinas de D_no_blocked
	
	
	while(1)    //bucle del tiempo de la simulacion
	{				
	    
	    
	    consenso_link();    //en esta rutina elijo un link y con prob S_ij, n1 imita a n2		
	    
	    
	    
	    
	    if(((tpo % multiplo)==0.0) || (flag==1))   //cuento y escribo en los ficheros sólo cada multiplo veces y tb la última vez cuando ya hay consenso global
	    {
		
		contar_blocked_links();						
		
		contar_links_S1();
		
		marcar_link_consenso(); // ahora marco los nodos de un link como 1 si están en de acuerdo                                    
		find_Smax();              // y como  0 si no. tambien obtengo NN[] y Phase[][]
		
		
		S_max=GComp;
		Nclusters=NCLUSTERS;
		
		
		for (rasgo=1;rasgo<=f;rasgo++)     //calculo las f Sf_max
		{
		    marcar_link_consenso_rasgo(); // ahora marco un link  como 1,2,...f   si están 
		    find_Smax();                       // de acuerdo en ese rasgo concreto, y 0 si no 
		    
		    
		    Sf_max[rasgo]=GComp;
		    NCL[rasgo]=NCLUSTERS;
		    
		}  
		
		
		
		if(n_links_blocked==N_links_max)
		{ 
		    flag=1;         //escribir los archivos la ultima vez	    	   	    	    	
		}
		
		
		
		printf("S_max(t=%d)=%d  n_links_blocked=%d  n_links_S1=%d  N_cl=%d \n",contador,S_max,n_links_blocked,n_links_S1,NCLUSTERS);
		
		aux_links=N_links_max;
		aux_N=N;
		fich=fopen(file,"at");
		fprintf(fich,"%d   %f   %f   %f",contador,n_links_blocked/aux_links,n_links_S1/aux_links,S_max/aux_N);	
		for(i=1;i<=f;i++)       
		{
		    fprintf(fich,"   %f",Sf_max[i]/aux_N);	
		}
		fprintf(fich,"\n");	
		fclose(fich);	  
		
		
		
		fich1=fopen(file1,"at");
		fprintf(fich1,"%d   %d",contador,Nclusters);	
		for(i=1;i<=f;i++)       
		{
		    fprintf(fich1,"   %d",NCL[i]);	
		}
		fprintf(fich1,"\n");	
		fclose(fich1);	 
		
		
		
		
		contador++;		
		
	    }
	    
	    
	    
	    if( flag==1)
	    { 	 
		
		ver_nodos_frontera();
		
		break;       // acabar la simulacion 
	    }	    	    	    
	    
	    
	    tpo++;    //pasos de tiempo reales de la simulacion


	}                    //fin de la simulacion
	
	
	marcar_link_consenso(); // ahora marco los nodos de un link como 1 si están en de acuerdo                                    
	find_Smax();   

	fich5=fopen(file5,"at");    //escribo el tamaño de los clusters presentes al final de cada iter
	fprintf(fich5,"iter:%d       ",iter);
	for(i=1;i<=N;i++)
	{
	    if(sizes[i]!=0)
	    {
		
		fprintf(fich5,"%d   ",sizes[i]);			
	    }
	}
	fprintf(fich5,"\n");	
	fclose(fich5);	 
	
	






//	getchar();



    }   //fin bucle iteraciones
    
    
/* Normalizo y escribo la PK_tot */
    
    for( i = 1; i <= N; i++)
    {
	PK_tot[i] = PK_tot[i]/Niter;

	Pk_frontera_1_tot[i]=Pk_frontera_1_tot[i]/Niter;
	Pk_frontera_2_tot[i]=Pk_frontera_2_tot[i]/Niter;
    }
    
    
    fich3=fopen(file3,"wt");   
    for(i=1;i<=N;i++)       
    {
	fprintf(fich3,"%d   %f\n",i,PK_tot[i]);	
    }  
    fclose(fich3);	 
    
    
    
    
    
    
    fich4=fopen(file4,"wt");   
    for(i=1;i<=N;i++)       
    {
	fprintf(fich4,"%d   %f   %f\n",i,Pk_frontera_1_tot[i],Pk_frontera_2_tot[i]);	
    }  
    fclose(fich4);	     


    
    
    printf("\n\n");
    for(i=1;i<=N;i++)
    {
	histogr_sizes_tot[i]= histogr_sizes_tot[i]/Niter;
//	printf("%i   %f\n",i,histogr_sizes_tot[i]);

    }





    fich7=fopen(file7,"wt");   
    for(i=1;i<=N;i++)       
    {
	fprintf(fich7,"%i   %f\n",i,histogr_sizes_tot[i]);	
    }  
    fclose(fich7);	 





    
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



/////////////////////////////
////////////////////////////
//////////////////////




void construir_red()
{
    double aux;  
    
    ValorA=3;
    steps=N-m_o;
    
    
    
    
    for(i = 0; i<=N; i++)
    {
	k[i] = 0;
	k_PA[i] = 0;
	A[i] = 0;
	
	for(j = 0; j <= m_o; j++)
	    M[i][j] =0;
    }
    
    
    for(j = 1; j <=m_o ; j++)    // inicializo los m_o primeros nodos conectados entre si
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
	
//	printf("\nnodo: %d \n",s+m_o);
	
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
    
    
//calculo la P(k)
    
    for(i=0; i<N; i++)
	PK[i] = 0;
    
    
    for(j=1; j<=N; j++)
        PK[k[j]]++;
    
    
    
    
    /* Normalizo PK */
    aux=N;
    for(i=1;i<=N;i++)
    {
	PK[i] = PK[i]/aux;
	PK_tot[i]+=PK[i];
    }         
    
    
 

}  






//////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////



void matriz_D()        //contiene todos los pares de vecinos (segun orden de creacion del link)
{
  int i,j,h;
  
  for(i=0;i<=m*N;i++)     //ojo!  esto es de 0 a 2N porque m=2  !!!!!
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
	    {                          //OJO!!!  VALIDO SOLO SI m=2
	      D[h][1]=i;
	      D[h][2]=M[i][j];
	      h++; 
	    }
	}
      
    }

  n_links=0;   //guardo el numero total de links   
  for(i=1;i<=m*N;i++)
  {
      if(D[i][1]!=0)
      {
	  n_links++;	  
      }
//printf("n_links:%d\n",n_links);      
  }
 
 
 
}  



////////////////////////////////////
//////////////////////////////////////////
////////////////////////////////////////






void contar_blocked_links()
{
    
    int i,j,n1,n2;
   
    
    n_links_blocked=0;
    
    
    for(i=1;i<=N_links_max;i++)    
    {
	n1=D[i][1];         //OJO, ESTO VALE SI m=2, si no, hay que modificarlo!!
	n2=D[i][2];


	S_ij=0.0;
	for(j=1;j<=f;j++)  
	{
	    if(v_i[n1][j]==v_i[n2][j])
	    {
		S_ij++;	   
	    }
	}	

	if(S_ij==0.0 || S_ij==f)   
	{
	    n_links_blocked++;
	}
    }
    /* printf("n_links_bloqued%d\n",n_links_blocked); 
       getchar();*/
}



void contar_links_S1()
{
    
    int i,j,n1,n2;
   
    
    n_links_S1=0;
    
    
    for(i=1;i<=N_links_max;i++)    
    {
	n1=D[i][1];         //OJO, ESTO VALE SI m=2, si no, hay que modificarlo!!
	n2=D[i][2];


	S_ij=0.0;
	for(j=1;j<=f;j++)  
	{
	    if(v_i[n1][j]==v_i[n2][j])
	    {
		S_ij++;	   
	    }
	}
	   
	
	if(S_ij==f)   
	{
	    n_links_S1++;
	}
    }
    
}




////////////////////////////////
///////////////////////////////
////////////////////////////////




void find_Smax()    //basado en la rutina TopologiaCP (Mi convenio era: -1 fluct, 0 coop, +1 defect) 
{                         //ahora =i es un link de consenso en el rasgo i, y =0 es no consenso

    int J, I,h,ii,kk,i,gg,todos_nodos[N+1],aislados,analyse[N+1],visitado[N+1],size[N+1];
    double histogr_sizes[N+1],aux;
    
    
    //nota: matrices NN[] y Phase inicializadas y obtenidas en la rutina:  marcar_links_consenso() y  marcar_links_consenso_rasgo()
    
    
    
    
    
    for(I=0;I<=N;I++)
    {
	visitado[I]=0;
	size[I]=0;
	histogr_sizes[I]=0;
	todos_nodos[I]=0;
    }


    ii=0;
    aislados=0;
    NCLUSTERS=0;  
    GComp=0;

    for(I=1;I<=N;I++)     //recorro la red
    {
	kk=0;
	for(J=1;J<=N;J++)
	{
	    analyse[J]=0;      //guardara los nodos que forman parte del cluster actual
	}

	if(NN[I]>0) //si no es un nodo asilado
	{
	    kk=1;
	    if(visitado[I]==0)   //si el nodo no ha sido visitado antes
	    {
		visitado[I]=1;   //lo marco como visitado
		analyse[kk]=I;   //lo marco como perteneciente al cluster actual
		todos_nodos[I]=I;



		for(jj=1;jj<=kk;jj++)
		{
		    g=analyse[jj];    //el nodo en el que estoy
		    
		    for (i=1;i<=NN[g];i++)    //miro sus vecinos
		    {
			h=Phase[g][i];
			
			if(visitado[h]==0)
			{
			    visitado[h]=1;
			    kk++;
			    analyse[kk]=h;
			    todos_nodos[h]=h;
			}
		    }
		}
	    }

	}	
	else    //si el nodo esta aislado
	{
	    visitado[I]=1; 
	    aislados++;  
	    histogr_sizes[1]++;   
	    
	    /* if(flag==1)    //solo al acabar la simu
	    {
		printf("aislado:%d\n ",I);		
		}*/
	}	
	
	///////////////////////////////////////////cuando acabo con el cluster actual:


	if(kk>GComp)      //guardo el mayor de los clusters como GComp del sistema
	{
	    GComp=kk;
	}
	

	if(kk>1)
	{
	    NCLUSTERS++;  
	}


	if(flag==1)    //solo al acabar la simu
	{	
	    if(kk>1)
	    {
		/*	printf("\ntamaño:%d       ",kk);
		for(J=1;J<=kk;J++)
		{
		    printf("%d   ",analyse[J]);
		}
		printf("\n");*/
		
		
		ii++;
		size[ii]=kk;     //guardo el tamaño del cluster que he encontrado
		histogr_sizes[kk]++;   //acumulo para el histograma de tamaños
		
		
		
		
		for(jj=1;jj<=kk;jj++)   
		{			
		    fich6=fopen(file6,"at");
		    fprintf(fich6,"nodo:%d (k:%d)          ",analyse[jj],k[analyse[jj]]);
		    
		    for(gg=1;gg<=f;gg++)       
		    {
			fprintf(fich6,"   %d",v_i[analyse[jj]][gg]);	
		    }
		    fprintf(fich6,"\n");	
		    fclose(fich6);	 		    		   
		    
		}	
		
		fich6=fopen(file6,"at");
		fprintf(fich6,"\n\n");	
		fclose(fich6);
		
	    }	
	}    
	

    }   /////////////////////////////////fin del bucle a todos los nodos de la red


 
    if(flag==1)      //solo cuando la iter ya ha acabado
    {
	aux=aislados+NCLUSTERS;    	
	for(I=1;I<=N;I++)
	{	   	    
	    histogr_sizes[I]=histogr_sizes[I]/aux;	    	   	    
	    histogr_sizes_tot[I]= histogr_sizes_tot[I] + histogr_sizes[I];

//	    printf("%d    %f   (%f)\n",I,histogr_sizes[I],aux);
	    
	}

	
	for(jj=1;jj<=N;jj++)        //guardo tb los rasgos de los nodos aislados
	{
	    if(todos_nodos[jj]==0)
	    {			
		fich6=fopen(file6,"at");
		fprintf(fich6,"nodo:%d (k:%d)          ",jj,k[jj]);
		
		for(gg=1;gg<=f;gg++)       
		{
		    fprintf(fich6,"   %d",v_i[jj][gg]);	
		}
		fprintf(fich6,"\n\n");	
		fclose(fich6);	 
				   
	    }
	}	
	printf("GComp:%f    NCLUSTERS:%d    aislados:%d",GComp,NCLUSTERS,aislados);	
    }
    
    




    
}

    

//////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////



void marcar_link_consenso ()
{
    int i,j,n1,n2;  
    
    
    for(i=1;i<=N_links_max;i++)   //inicializo el vector que guardará la similitud de un link (en total, no por rasgos)
    {
	Simil[i]=0;
    }
    
    
    for(j=1;j<=N;j++)    
    {
	NN[j]=0;              //NN[] es k[] pero solo enlaces de consenso
	for(i=1;i<=N;i++)
	{
	    Phase[j][i]=0;  //Phase[][] es C[][] pero solo entran links de consenso
	}
    }


    for(j=1;j<=N_links_max;j++)  //recorro la matriz de links
    {
	n1=D[j][1];    //OJO!!!  VALIDO SOLO SI m=2, si no, hay que añadir lineas
	n2=D[j][2];
	
	
	Simil[j]=0;    //vector que guarda la similitud en un link (el tamagno de vector es m*N)
	
	for(i=1;i<=f;i++)  //calculo S_ij
	{	       
	    if(v_i[n1][i]==v_i[n2][i])
	    {
		Simil[j]++;
	    }
	}
	

	if(Simil[j]==f)   //construyo las matrices para el subgrafo de consenso    
	{	  
	    NN[n1]++;
	    NN[n2]++;

	    Phase[n1][NN[n1]]=n2;
	    Phase[n2][NN[n2]]=n1;
	    	   
	} 	
	


    }
       

}



///////////////////////////////////////
////////////////////////////////////////
///////////////////////////////


void marcar_link_consenso_rasgo ()
{
    int i,j,n1,n2;  
    
    
    for(i=1;i<=N_links_max;i++)   //inicializo el vector que guardará la similitud de un link respecto al rasgo que este mirando
    {	
	Simil[i]=0;	      
    }
    
    
    for(j=1;j<=N;j++)   
    {
	NN[j]=0;                //NN[] es k[] pero solo enlaces de consenso
	for(i=1;i<=N;i++)
	{
	    Phase[j][i]=0;  //Phase[][] es C[][] pero solo entran nodos que comparten links de consenso
	}
    }        



//estamos mirando un rasgo concreto (rasgo =1,2,...f)

	for(j=1;j<=N_links_max;j++)  //recorro la matriz de links
	{
	    
	    n1=D[j][1];    //OJO!!!  VALIDO SOLO SI m=2, si no, hay que añadir lineas
	    n2=D[j][2];
	    
	    if(v_i[n1][rasgo]==v_i[n2][rasgo])    //si hay consenso en el link sobre el rasgo
	    {
		Simil[j]++;	      
	    }


	    if(Simil[j]==1)
	    {
		NN[n1]++;
		NN[n2]++;

		Phase[n1][NN[n1]]=n2;
		Phase[n2][NN[n2]]=n1;
	    }
	}
 
    

}

/////////////////////////////////////////////////////
//////////////////////////////////////////////////
/////////////////////////////////////////



void consenso_link()   // rutina usada en el programa consenso_Sf_max.c
{

    double aux_w, aux_r, aux;
    int i,w,r,flag;
    int nodo1,nodo2,n1,n2;
    
      
    
    aux_w=FRANDOM;     //elijo un link y guardo los nodos que une
    aux_w=aux_w*N_links_max;
    w=aux_w+1;   //por que genero numeros de 0 a  N*m -1

    nodo1=D[w][1];       
    nodo2=D[w][2];
    
    //printf("nodos: %d y %d\n",nodo1,nodo2);

    aux_r=FRANDOM;
    if(aux_r<0.5) //elijo al azar cuál imitará a cuál
    {
	n1=nodo1;
	n2=nodo2;
    }	
    else
    {
	n1=nodo2;
	n2=nodo1;	
    }

   S_ij=0.0;
    for(i=1;i<=f;i++)  //calculo la probabilidad de cambio (S_ij)
    {
	if(v_i[n1][i]==v_i[n2][i])
	{
	    S_ij++;	   
	}
    }
    aux=f;
    S_ij=S_ij/aux;   
 


    aux_w=FRANDOM; 
    if(aux_w<S_ij)    //prob de que n1 imite a n2 en un rasgo que no tengan igual
    {	
	flag=0;

	while(flag==0)
	{	    
	    if(S_ij!=0.0 && S_ij!=1.0 )
	    {    
		aux_r=FRANDOM;
		aux_r=aux_r*f; 
		r=aux_r+1;       //elijo un rasgo  (ojo!!!! la matriz se llena de 1 a f y el generador va de 0 a f-1)
	    
	
		if(v_i[n1][r] != v_i[n2][r])
		{		   
		    v_i[n1][r]=v_i[n2][r];  //n1 imita a n2		  
		    flag=1;	
		}
		
	    }
	    else              //  PARA QUE NO SE QUEDE ATASCADO EN UN BUCLE QUE ESTE BLOQUEADO
	    {			    		
		flag=1;		    		    		
	    }
	}
    }
        
    
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////




void ver_nodos_frontera()
{

    int i,j,ii,n1,n2,contador_S0,N1,N2;
    double Pk_frontera_1[N+1],Pk_frontera_2[N+1],aux;
    
    for(i=1;i<=N;i++)
    {
	Pk_frontera_1[i]=0;    //pk de nodos frontera con al menos un link de S=0
	Pk_frontera_2[i]=0;     // id. al menos 2
    }

    N1=N2=0;
    
    for(i=1;i<=N;i++)    //recorro toda la red
    {
	contador_S0=0;   //recuento del numero de links con S=0 del nodo

	n1=i;
	for(j=1;j<=k[i];j++)       //recorro todos los vecinos del nodo
	{
	    n2=C[i][j];
	    
	    S_ij=0.0;
	    for(ii=1;ii<=f;ii++)  //calculo la probabilidad de cambio (S_ij)
	    {
		if(v_i[n1][ii]==v_i[n2][ii])
		{
		    S_ij++;	   
		}
	    }
	    

	    if(S_ij==0)
	    {
		contador_S0++;
	    }
	}

	if(contador_S0>=1) //tipo 1
	{
	    Pk_frontera_1[k[i]]++;
	    N1++;              //la norma de los nodos frontera tipo 1 (=num tot nodos tipo 1)
	}
	
	if(contador_S0>=2)  //tipo 2
	{
	    Pk_frontera_2[k[i]]++;
	    N2++;         //la norma de los nodos frontera tipo 2 (=num tot nodos tipo 2)
	}

   } //fin bucle nodos de la red


    
    /* aux1=N1;             //NO ESTOY SEGURA DE LA NORMALIZACIÓN!!!!!!
    aux2=N2;
    for(i=1;i<=N;i++)
    {
	
	if(N1!=0)
	{
	    Pk_frontera_1[i]= Pk_frontera_1[i]/aux1;
	    Pk_frontera_1_tot[i]+=Pk_frontera_1[i];
	}
	
	
	if(N2!=0)
	{
	    Pk_frontera_2[i]= Pk_frontera_2[i]/aux2;	    
	    Pk_frontera_2_tot[i]+=Pk_frontera_2[i];
	}
	

	}     */



    aux=N;
    for(i=1;i<=N;i++)
    {		
	Pk_frontera_1[i]= Pk_frontera_1[i]/aux;
	Pk_frontera_1_tot[i]+=Pk_frontera_1[i];		
	
	Pk_frontera_2[i]= Pk_frontera_2[i]/aux;	    
	Pk_frontera_2_tot[i]+=Pk_frontera_2[i];		
	
    }     
    


    
    

}
















