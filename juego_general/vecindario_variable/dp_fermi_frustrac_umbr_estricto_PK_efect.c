//22-10-10: programa basado en    dp_fermi_vecindario_tpos_ks_frustrac_umbr_estricto.c
// pero sólo para calcular la P(k) efectiva, en funcion de k*. NO NECESITO JUGAR, SOLO SELECCIONAR EL VECINDARIO!!!



# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>


# define Niter   10    //estadistica

# define alfa 0.0


# define N 4000 //tamaño de la red
# define m_o  2      //nodos inicialmente unidos
# define m   2      //nodos nuevos añadidos a cada paso de tiempo

 



# define KSTAR      80      // numero de vecinos con los que juego/me comparo (entre m y N)




 

/////////////////////////Generacion de numeros aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
////////////////////////////////////////////////////////





char nombre3[256],nombre2[256],nombre1[256];

int iter; 
int k[N+1],k_PA[N+1],A[N+1];                    //conectividad topológica y de Preferential Attachment (PA)
int   M[N+1][m+1], unido[m+1], Cuidado[N+1];         //  unido[] guarda a quién se ha unido i , y cuidado[] quien le ha lanzado ya un link
double P[N+1], P_prov[N+1];        //P_prov[] es la provisional para manipularla
double PK[N+1], PK_tot[N+1];        //  para la distribucion P(k)
double r,v;                                  //para guardar los aleatorios
int  i, j, jj,  w,g, d, q, x[N+1],y,z, C[N+1][N+1],n, ValorA, steps,s,flat;
int norma, norma2, norma_aleat,tipo[m_o+1];
double dnorma_aleat;
int si[N+1], cont,cont2;
int k_star[N+1],C_star[N+1][N+1];

double Pk_eff[N+1],Pk_eff_tot[N+1];

double  kmedia_top, kmedia_eff;  //media sobre nodos

double  kmedia_top_tot, dispers_k_top_tot, kmedia_eff_tot, dispers_k_eff_tot; //media sobre iter


double resta2_ktopol, resta2_keff;



void inicia_rand(int semilla);
void construir_red();
void selecciono_vecindario();
void pk_efectiva();
void histograma_pk();
void media_dispersion_k ();

FILE *escribe3,*escribe2,*escribe1;




int main()
{
  
  inicia_rand(time(0));
  // inicia_rand(26);

    //para la Pk topologica
    sprintf(nombre1,"Pk_topologica_N%d_%diter_frustr.dat",N,Niter); 
    escribe1=fopen(nombre1,"wt");  
    fclose(escribe1);
    
    
    //para la Pk efectiva
    sprintf(nombre2,"Pk_efectiva_kstar%d_N%d_%diter_frustr.dat",KSTAR,N,Niter); 
    escribe2=fopen(nombre2,"wt");  
    fclose(escribe2);
    
  //para la <k> y <k2>
    sprintf(nombre3,"Kmed_DispersK_N%d_%diter_frustr.dat",N,Niter); 
    escribe3=fopen(nombre3,"at");  
    fclose(escribe3);
    


    for(i=1;i<=N;i++)
      {
	Pk_eff[i]=0;

      }
    
    
     kmedia_top_tot= dispers_k_top_tot= kmedia_eff_tot= dispers_k_eff_tot=0;

	for(iter=1;iter<=Niter;iter++)
	{

 
	  construir_red();
	  histograma_pk();

	    
	  selecciono_vecindario();	  
	  pk_efectiva();

	




	  media_dispersion_k ();


	}   //fin bucle estadistica
	
	
	for(i=1;i<=N;i++)
	  {
	    Pk_eff[i]=Pk_eff[i]/(Niter*N);
	    PK[i]=PK[i]/(Niter*N);
	  }
	



	  kmedia_top_tot = kmedia_top_tot/Niter;
	  dispers_k_top_tot = dispers_k_top_tot/Niter;      //la media (sobre Niter)) de las dispers sobre nodos
	  kmedia_eff_tot = kmedia_eff_tot/Niter;
	  dispers_k_eff_tot = dispers_k_eff_tot/Niter;








 printf("\n\nKSTAR:%d   kmed_top:%f   dispers_k_top:%f   kmed_eff:%f   dispers_k_eff:%f \n",KSTAR, kmedia_top_tot, dispers_k_top_tot,kmedia_eff_tot,dispers_k_eff_tot );


	//// escribo lal Pk topologica
	escribe1=fopen(nombre1,"at");
	for(i=1;i<=N;i++)
	  {
	    fprintf(escribe1,"%d  %f\n",i,PK[i]);      
	  }
	fclose(escribe1);
	



//// escribo lal Pk efectiva
	escribe2=fopen(nombre2,"at");
	for(i=1;i<=N;i++)
	  {
	    fprintf(escribe2,"%d  %f  %d\n",i,Pk_eff[i],KSTAR);      
	  }
	fclose(escribe2);
	



//// escribo lal Pk efectiva
	escribe3=fopen(nombre3,"at");
	fprintf(escribe3,"%d  %f  %f  %f  %f\n",KSTAR, kmedia_top_tot, dispers_k_top_tot,kmedia_eff_tot, dispers_k_eff_tot );      	  
	fclose(escribe3);
	





	exit(0);


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

}  



/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////





void selecciono_vecindario()
{
   
    int i,j,kmin,elegido[N+1],y,w,minimo,available[N+1],n_max_interacc,n_frustraciones,suma_vecinos_dispon; 
    double r,aux;
    int label[N+1],label_ini[N+1];    //guarda de k* a cero, el numero de interacciones libres que le quedan a un nodo


    n_max_interacc=0;
    for(i=1;i<=N;i++)  //establezco el numero de interacciones disponibles en principio para cada nodo
    {
	if(k[i]>=KSTAR)
	{
	    label_ini[i]=KSTAR; //para saber siempre cuántas interacciones posibles tenia al ppio y 
	    label[i]=KSTAR;           //  cuántas le quedan al final (si =0, las ha usado todas, si no, es que hay frustracion)
	}
	else
	{
	    label_ini[i]=k[i];
	    label[i]=k[i];
	}
	k_star[i]=0;    //guardaré el numero real de interacciones que acumula (in +out)


	for(j=1;j<=N;j++)
	{
	    C_star[i][j]=0;
	}


	elegido[i]=1;    //marco =0 los nodos a los que ya les he seleccionado sus vecinos

	n_max_interacc+=label_ini[i];

    }

 

 for(jj=1;jj<=N;jj++)//bucle para repetir la operacion N veces
    {
	for(i=1;i<=N;i++)  //establezco el numero de interacciones disponibles en principio para cada nodo
	{
	    available[i]=1;  //para evitar enlaces dobles
	}
	


	kmin=N;
	minimo=0;
	for(i=1;i<=N;i++) //busco el nodo con menor k, que será elegido para establecer sus k*(o menos) conexiones.
	{
	    if(k[i]<kmin && elegido[i]==1)
	    {	    
		kmin=k[i];
		minimo=i;		
	    }
	}
	
	elegido[minimo]=0; //lo marco para no volver a elegirlo  ¿se queda en la i que yo quiero al hacer el break????	
	
//	printf("elegido:%d  (k:%d)",minimo,k[minimo]);  
	//getchar();
	

	i=minimo;
	






///////////////lo primero miro si el nodo i va a sufrir frustracion 

	suma_vecinos_dispon=0;
	for(j=1;j<=k[i];j++)
	{
	    y=C[i][j];
	    if(label[y]>0)
	    {
		suma_vecinos_dispon++;
	    }
	}	    
	if(suma_vecinos_dispon>=label_ini[i])    // no sufrira frustracion
	{

	    //printf("  (no frustrac)");
	    while(label[i]>0 ) // elijo al azar hasta completar el max de vecinos posibles
	    {
		
		r=FRANDOM; //elijo un vecino al azar
		r=r*k[i];
		w=(int)r+1;
		y=C[i][w];
		
		if(label[y]>0 && label[i]>0 && available[y]!=0) //si es un vecino disponible (con links libres 
		{                                            //    y no elegido previamente por ese mismo nodo)
		    k_star[i]++;
		    k_star[y]++;			
		    
		    
		    C_star[i][k_star[i]]=y;     //asi garantizo que tanto i como y acumulen luego beneficios
		    C_star[y][k_star[y]]=i;
		    
		    label[y]--;        //ambos gastan una de sus posibles interacciones
		    label[i]--;

		    
		    
		    available[y]=0;  //para evitar enlaces dobles
		    
		    
		    //printf("%d-%d  ahora son amigos, con k_star respect. %d y %d   y label %d  y %d\n",i,y,k_star[i],k_star[y],label[i],label[y]);  
		    
		}	
		
	    }  //fin del while
	}
	else     //si sufrira frustracion  -->> elijo directamente por orden los vecinos disponibles con los que interact.
	{
	    //printf("  (no frustrac)");
	    for(j=1;j<=k[i];j++)
	    {
		y=C[i][j];
		if(label[y]>0 && label[i]>0 && available[y]!=0)
		{
		    k_star[i]++;
		    k_star[y]++;			
		    
		    
		    C_star[i][k_star[i]]=y;     //asi garantizo que tanto i como y acumulen luego beneficios
		    C_star[y][k_star[y]]=i;
		    
		    label[y]--;        //ambos gastan una de sus posibles interacciones
		    label[i]--;
		    
	    
		    available[y]=0;  //para evitar enlaces dobles
		    
		}
	    }
	}
	
	
	
	
    }
    

    n_frustraciones=0;
    for(i=1;i<=N;i++)
    {
	n_frustraciones+=label[i];

    }


    aux=n_frustraciones;
    aux=aux/n_max_interacc;
    //printf("%d/%d=%f \n",n_frustraciones,n_max_interacc,aux);

    //getchar();


    /*   for(i=1;i<=N;i++)
    {
      printf("k(%d):%d   kstar:%d\n",i,k[i], k_star[i]);  
	 getchar();
	 }*/

   
    

/*    if (iter==1)
    {	
	escribe1=fopen(nombre1,"wt");
	for(i=1;i<=N;i++)
	{
	    fprintf(escribe1,"%d  %d  %d\n",i,k[i],k_star[i]);      
	}
	fclose(escribe1);
    }
*/  
    
    
}



////////////////////////////////////
///////////////////////////////////////
///////////////////////////


void pk_efectiva()
{

  int i;
  

 for(i=1;i<=N;i++)
    {
      Pk_eff[k_star[i]]++;
    }


 for(i=1;i<=N;i++)
   {
     Pk_eff_tot[i]+=Pk_eff[i];
   }


}



///////////////////////
///////////////////////
///////////////////////////

void histograma_pk()
{

  int i;
  

 for(i=1;i<=N;i++)
    {
      PK[k[i]]++;
    }


 for(i=1;i<=N;i++)
   {
     PK_tot[i]+=PK[i];
   }

}





///////////////////////
///////////////////////
///////////////////////////


void media_dispersion_k ()
{

  int i;
 
  kmedia_top=0.;
  
  kmedia_eff=0.;
 

  for(i=1;i<=N;i++)
    {

      kmedia_eff=kmedia_eff + k_star[i];
    
      kmedia_top=kmedia_top + k[i];
  
    }
 

  kmedia_top=kmedia_top/N;     
  kmedia_eff=kmedia_eff/N;


  kmedia_top_tot += kmedia_top;
  kmedia_eff_tot += kmedia_eff;



 
resta2_ktopol=resta2_keff=0.;

  for(i=1;i<=N;i++)
    {
      resta2_ktopol += (k[i]-kmedia_top)*(k[i]-kmedia_top);       
      resta2_keff += (k_star[i]-kmedia_eff)*(k_star[i]-kmedia_eff);    
      
    }


  resta2_ktopol=resta2_ktopol/(N-1.0);
  resta2_keff=resta2_keff/(N-1.0);


  resta2_ktopol=sqrt(resta2_ktopol);
  resta2_keff=sqrt(resta2_keff);

  printf("KSTAR:%d   kmed_top:%f (%f)    kmed_eff:%f (%f) \n",KSTAR, kmedia_top, resta2_ktopol,kmedia_eff,resta2_keff);

  dispers_k_top_tot +=resta2_ktopol;
  dispers_k_eff_tot +=resta2_keff;


}
