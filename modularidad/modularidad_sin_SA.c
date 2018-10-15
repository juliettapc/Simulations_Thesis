/////// 13 Ag 2010:  programa para obtener la estructura de comunidades de una 
////red, mediante maximizacion de la modularidad (por Simulated Annealing). 
//la red se lee de un fichero.




# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>


/////////////////////////Generacion de numeros aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
////////////////////////////////////////////////////////




//numero de filas del archivo de entrada == numero total links de la red:
# define filas 4598    //VER  ARCHIVO!!!!!
# define N  256          //VER  ARCHIVO!!!!! 


# define T_max    20.0           //temperatura para el SA
# define T_min    0.001
# define delta_T  0.005


int D[filas+1][3];  //matriz de pares de links
int k[N+1],C[N+1][N+1];
double PK[N+1];        //  para la distribucion P(k)



int nodo[N+1],nodo_aux[N+1];   //guarda el indice de la comunidad a la que pertenece el nodo
                             // ahora y en el posible cambio, respectivamente
int nodo_elegido,NumCom,ii,jj,i, desaparecido;

int nodo_escin[N+1],nodo_aux_escin[N+1];    
  //indices de los nodos (y aux para el cambio) pero en referencia  solo a la comunidad que voy a partir

double T,M,M_old,aux,prob,exponencial;


int  NumCambiosIndiv,NumCambiosColect,flag_escin,flag_fus;     
    

int comunidad[N+1],size; //guardo los nodos que pertenecen a la comunidad a escindir (y su tamagno)

int tpo;


void inicia_rand(int semilla);
void leer_red();
void construir_C_k();
void particion_inicial(); 
double calculo_modularidad();
void cambio_indiv(int);
void escindir();
void fusionar();
void reasignar_indices();
void cambio_indiv_escin(int,int,int);
double calculo_modularidad_anidado(int);


char nombre1[256],nombre2[256];
  
FILE *archivo1;
FILE *escribe2;




int main()
{     
  
  
  
  printf("Programa para obtener la estructura de comunidades de una red \n maximizando la modularidad  -(sin SA)-.\n\n");
  
  
  
  /////////archivo de entrada
  sprintf(nombre1,"256_4_4_4_13_18_p.net");  
  
  
  
  /////////archivo de salida
  sprintf(nombre2,"Comunid_N%d_Nlinks%d_T%.2f-%.2f_sinSA.dat",N,filas,T_max,T_min); 
  escribe2=fopen(nombre2,"wt");  
  fclose(escribe2);
  
  
  
  
  inicia_rand(time(0));
  leer_red();
  
  
  
  //segun establecen en el paper: NumCambiosIndiv=N*N   y  NumCambiosColect=N
  
  
  
  particion_inicial();    //2comunidades: por la mitad, aleatoriamente
  NumCom=2;     //numero de comunidades 
  M_old=calculo_modularidad();
  
  
  
  tpo=0;
  T=T_max;
  
  
  escribe2=fopen(nombre2,"at");
  fprintf(escribe2,"%d  %f  %f  %d\n",tpo,T,M_old,NumCom);      
  fclose(escribe2);
  
  
         
       printf("\nT:%f\n",T); 
      
      for(ii=1;ii<=(N*N+N);ii++) // N*N cambios individuales y N cambios colesctivos, intercalados
	{
	  
	  prob=FRANDOM*(N*N+N);
	  
	  if(prob>N)     ////////////cambio indiv.
	    {
	      
	      aux=FRANDOM*N;
	      nodo_elegido=aux+1;      //elijo el nodo a mover
	      
	      
	      cambio_indiv(nodo_elegido);     //modifica el indice de:  nodo_aux[nodo_elegido]
	      
	      
	      
	      M=calculo_modularidad();  //de la particion nueva  
	      
	      //printf("T:%f   M:%f   NumCom:%d\n",T,M,NumCom);
	      
	      if(M>M_old)
		{
		  nodo[nodo_elegido]=nodo_aux[nodo_elegido]; //acepto el cambio 
		  //printf("  -->   acepto cambio indiv\n");


	      
		  escribe2=fopen(nombre2,"at");
		  fprintf(escribe2,"%d  %f  %f  %d\n",tpo,T,M,NumCom);      
		  fclose(escribe2);
		  
		  tpo++;
		  


		  M_old=M;	      //SOLO SI ACEPTO EL CAMBIO, CAMBIO M_OLD?????????


		}
	      /*	       else   //acepto el cambio con una peq. prob. (para escapar de max. locales)
		{
		  prob=FRANDOM;
		  exponencial=exp((M-M_old)/T);
		  if(prob<exponencial)
		      {
			nodo[nodo_elegido]=nodo_aux[nodo_elegido]; //acepto el cambio 
			//printf("  -->   acepto cambio indiv (peq.prob.)\n");

			escribe2=fopen(nombre2,"at");
			fprintf(escribe2,"%d  %f  %f  %d\n",tpo,T,M,NumCom);      
			fclose(escribe2);
			
			tpo++;
			

			M_old=M;	      //SOLO SI ACEPTO EL CAMBIO, CAMBIO M_OLD?????????

		      }
		      }	 */     	      	      
	      
	      
	     	      	      
	      
	      
	    }   //fin cambio indiv.	  	  	  	  	  
	  else      ///////////cambios colectivos 
	    {	      

	      flag_escin=flag_fus=0;  //indicadores del tipo de cambio colectivo realizado
	      
	      prob=FRANDOM;
	      
	      if(prob<0.5)
		{
		  escindir();	
		  flag_escin=1;
		}
	      else
		{		  
		  fusionar();
		  flag_fus=1;
		  
		}
	      
	      
	      
	      M=calculo_modularidad();  //de la particion nueva 
	      
	      
	      // printf("T:%f   M:%f   NumCom:%d\n",T,M,NumCom);
	      
	      if(M>M_old)  //acepto el cambio 
		{
		  if(flag_fus==1)
		    {
		      NumCom --;
		      reasignar_indices ();
		      //   printf("  -->   acepto cambio colect. (fusion)\n");     

		    }
		  else		 
		    {
		      if(flag_escin==1)    //al escindir, ya reasigno los indices de comun. a los nodos
			{
			  NumCom++;
			  //  printf("  -->   acepto cambio colect. (escision)\n");
			 
			}		      		    
		    }
		  


	      
	      
		  escribe2=fopen(nombre2,"at");
		  fprintf(escribe2,"%d  %f  %f  %d\n",tpo,T,M,NumCom);      
		  fclose(escribe2);
		  
		  tpo++;
		  

		   M_old=M;	      //SOLO SI ACEPTO EL CAMBIO, CAMBIO M_OLD?????????
		  
		  for(i=1;i<=N;i++)
		    {
		      nodo[i]=nodo_aux[i];
		    }
		}
	      /* else       //acepto el cambio con una peq. prob. (para escapar de max. locales)
		{
		  prob=FRANDOM;
		  exponencial=exp((M-M_old)/T);
		  if(prob<exponencial)
		    {
		      if(flag_fus==1)
			{
			  NumCom--;
			  reasignar_indices ();
			  //printf("  -->   acepto cambio colect. (fusion peq. prob.)\n");
			}
		      else
			{
			  if(flag_escin==1) //al escindir, ya reasigno los indices de comun. a los nodos
			    {			      
			      NumCom++;
			      //printf("  -->   acepto cambio colect. (escision peq. prob.)\n");
			    }			  
			}		     
		      

		      escribe2=fopen(nombre2,"at");
		      fprintf(escribe2,"%d  %f  %f  %d\n",tpo,T,M,NumCom);      
		      fclose(escribe2);
		      
		      tpo++;		      		      


		       M_old=M;	      //SOLO SI ACEPTO EL CAMBIO, CAMBIO M_OLD?????????


		      //printf("  -->   acepto cambio colect. (peq.prob.)\n");
		      for(i=1;i<=N;i++)
			{
			  nodo[i]=nodo_aux[i]; 
			}
		    }
		    }  */	      
	      
	      
	         
	    } /////////fin cambio colectivo
	  
	  
	  
	  
	}    //fin del for de elegir los N*N+N cambios	  
      


      printf("T:%f  (t:%d)  M_old:%f   NUmCom:%d\n\n",T,tpo,M_old,NumCom);



  
  
  
  
  
  printf(" \nRESULTADO FINAL OPTIMIZACION:\n T:%f   M:%f   NumCom:%d\n\n",T,M,NumCom);
  
  
  
  
  
  
  exit(0);
}


  
  //////////////////////////////////////////
  /////////////////////////////////////////
////////////////////////////////////////



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




//////////////////////////////////////////
/////////////////////////////////////////
////////////////////////////////////////

void leer_red()
{

    int i,c1,c2,c3;
     

    for(i=1;i<=filas;i++)
    {
	D[i][1]=0;
	D[i][2]=0;
    }
   

    archivo1=fopen(nombre1,"r");  
    for(i=1;i<=filas;i++)
    {
      fscanf(archivo1,"%d   %d   %d\n",&c1,&c2,&c3);   //la tercera columna no me interesa

	D[i][1]=c1;
	D[i][2]=c2;
    }
    fclose(archivo1);


    /*  for(i=1;i<=filas;i++)   //comprobacion
    {
	printf("%d-%d\n",D[i][1],D[i][2]);
    }
    getchar();*/
}




//////////////////////////////////////////
/////////////////////////////////////////
////////////////////////////////////////


void construir_C_k()      //parto de la matriz de pares de links D[][]
{

    int i,j,d1,d2;

    for(i=1;i<=N;i++)
    {
	k[i]=0;
	for(j=1;j<=N;j++)
	{
	    C[i][j]=0;
	}
    }



    for(i=1;i<=filas;i++)
    {
//	printf("fila%d: ",i);

	d1=D[i][1];    //nodo 1 del link
	d2=D[i][2];    //nodo 2 del link

//	printf("%d-%d\n",d1,d2);


	k[d1]++;
	k[d2]++;
	
	/*printf("k1:%d,  k2:%d\n",k[d1],k[d2]);

	getchar();*/

	C[d1][k[d1]]=d2;
	C[d2][k[d2]]=d1;

    }
    
    /* for(i=1;i<=N;i++)   //comprobacion
    {

	printf("\n%d (k:%d):",i,k[i]);  
	for(j=1;j<=k[i];j++)
	{
	    printf("%d   ",C[i][j]);    
	}
	getchar();
	}*/
}



/////////////////////////////////
/////////////////////////////////////////////
/////////////////////////////////////

void particion_inicial()    //parto la red inicial por la mitad, aleatoriamente
{

  int i;
  double p;
  
  
  for(i=1;i<=N;i++)    //inicializo a cero los nodos 
    {                  //  (no pertenecen a ninguna comunidad)
      nodo[i]=0;
      nodo_aux[i]=0;
    }
  
  
  for(i=1;i<=N;i++)
    {
      p=FRANDOM;

      if(p<0.5)
	{
	  nodo[i]=1;    //pertenece a la comunidad 1
	  nodo_aux[i]=1;
	}
      else
	{
	  nodo[i]=2;    //pertenece a la comunidad 2
	  nodo_aux[i]=2;
	}
    }


}


/////////////////////////////////
/////////////////////////////////////////////
/////////////////////////////////////


void cambio_indiv (int elegido) //cambio al nodo elegido de su comunidad a otra  
{                            //OJO: es un cambio condicional, sujeto a aprobacion en el main!!!
  int auxiliar,i;
  double p;


  for(i=1;i<=N;i++)          //guardo la config original
    {
      nodo_aux[i]=nodo[i];
    }


  if(NumCom==2)
    {
      auxiliar=nodo[elegido];
      if(auxiliar==1)
	{
	  nodo_aux[elegido]=2;
	}
      else
	{
	  if(auxiliar==2)
	    {
	      nodo_aux[elegido]=1;
	    }
	}

      //      printf("Muevo nodo %d de %d a  %d\n.",elegido,nodo[elegido],nodo_aux[elegido]);

    }
  else      //si hay mas de dos comunidades
    {
      while(1)    //busco una comunidad distinta a la suya
	{
	  p=FRANDOM*NumCom;
	  auxiliar=p+1;

	  
	  if(auxiliar!=nodo[elegido])
	    {	      	      	      
	      nodo_aux[elegido]=auxiliar;	      
	      
	      //      printf("Muevo nodo %d de %d a  %d\n.",elegido,nodo[elegido],nodo_aux[elegido]);
	      
	      break;
	    }
	}
      

    }


}

//////////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////

void fusionar()     //elijo dos comunidades y las fusiono en una sola
{

  double p,q; 
  int aux_p,aux_q, i;

  
  for(i=1;i<=N;i++)
    {
      nodo_aux[i]=nodo[i];
    }
  

  while(1)       //elijo dos cumnun. distintas al azar
    {
      p=FRANDOM*NumCom;
      q=FRANDOM*NumCom;

      aux_p=p+1;
      aux_q=q+1;

      if(p!=q)
	{
	  break;
	}
    }

  for(i=1;i<=N;i++)
    {
      if(nodo_aux[i]==q)
	{
	  nodo_aux[i]=p; 
	}
    }


  desaparecido=q;   //indice de la comunidad que ya no existe 


}


////////////////////////////
////////////////////////////////
//////////////////////////////////

void escindir()     //elijo una comunidad y la parto en dos (haciendole un SA anidado)
{

  int i, ii,NumCambios,elegido;                            
  double p;              
  double T_escin=T_max;
  int comunidad1,comunidad2;
  double Mod,Mod_old;
  int auxiliar;


  p=FRANDOM*NumCom;    //elijo la comunidad a partir
   
  comunidad1=p+1;
  comunidad2=NumCom+1;

  size=0;
  for(i=1;i<=N;i++)     //obtengo el tamaño de la comunidad en cuestion
    {            //y la parto en dos aleatoriamente

      //OJO!!  TRABAJO CON INDICES RELATIVOS:      

      nodo_escin[i]=0;  // =aux o =NumCom+1 signif. que el nodo pertenece a una comunidad. o la otra, 
                        //   =0 a ninguna de las dos

      comunidad[i]=0;

      if(nodo[i]==comunidad1)
	{
	  size++;
	  
	  comunidad[size]=i;  //guardo una lista de los size nodos que pertenecen a la comunidad

	  p=FRANDOM;
	  
	  if(p<0.5)
	    {
	      nodo_escin[i]=comunidad1;    
	      nodo_aux_escin[i]=comunidad1;
	    }
	  else
	    {
	      nodo_aux_escin[i]=comunidad2;
	    }
	}
      
    }
  
  Mod_old=calculo_modularidad_anidado(comunidad1);  //de la particion aleatoria  

 

  NumCambios=size*size;


  while(T_escin>T)    //parto otra vez de una T alta, 
    {                //y la voy bajando hasta alcanzar la actual de sistema
                  	  	  
	  
      for(ii=1;ii<=NumCambios;ii++) //en este SA anidado sólo hay cambios individuales 
	{
	  
	  p=FRANDOM*size;
	  auxiliar=p+1;
	  elegido=comunidad[auxiliar];      //elijo el nodo a mover
	  
	  


	  //cambio un nodo de una parte a la otra:

      if(nodo_escin[elegido]==comunidad1)
	{
	  nodo_aux_escin[elegido]=comunidad2;
	}
      else
	{
	  if(nodo_escin[elegido]==comunidad2)
	    {
	      nodo_aux_escin[elegido]=comunidad1;
	    }
	}
      



	  // ESTE CALCULO ES SOLO RELATIVO A LA COMUNIDAD PARTIDA, NO A LA RED ENTERA!!:
	  
	  Mod=calculo_modularidad_anidado(comunidad1);  //de la particion nueva  
	  
	  if(Mod>Mod_old)
	    {
	      nodo_escin[elegido]=nodo_aux_escin[elegido]; //acepto el cambio 

	      Mod_old=Mod;	  	 //SOLO SI ACEPTO EL CAMBIO, CAMBIO MOD_OLD  ???????
	    }
	  /* else
	     {
	       p=FRANDOM;
	       expon=exp((Mod-Mod_old)/T);
	       if(p<expon)
		 {
		   nodo_escin[elegido]=nodo_aux_escin[elegido]; //acepto el cambio     
		   
		   
		   Mod_old=Mod;	  	 //SOLO SI ACEPTO EL CAMBIO, CAMBIO MOD_OLD  ???????
		   
		 }
		 }	*/ 	  	  
	  
	  	 
	  
	} /////////fin bucle  size*size  cambios individuales
      
      
      
      
      T_escin=T_escin-delta_T;
    }           ////////////////////   fin del bucle del SA anidado
  
  

  for(i=1;i<=N;i++) //copio en la matriz auxiliar general el resultado del SA para decidir si acepto
    {
      if(nodo[i]==comunidad1)
	{
	  nodo_aux[i]=nodo_aux_escin[i];
	}
    }
  
  
}


////////////////////////////
////////////////////////////////
//////////////////////////////////

  void reasignar_indices()       //solo si ha habido fusion. y se ha aceptado el cambio
  {               
    int i;
    
    //desaparecido: indice de la comun. que ya no existe

    for(i=1;i<=N;i++)
      {
	if(nodo[i]>=desaparecido)
	  {
	    nodo[i]--;   //le resto 1 a los indices de las comunidades posteriores
	  }
      }
    
    
  }
  


////////////////////////////
////////////////////////////////
//////////////////////////////////

double calculo_modularidad()
{

  int s, f,nodo1,nodo2;
  double links_intra[N+1],auxiliar,links_inter_intra[N+1],corchete,Modularidad;

 

  for(f=1;f<=N;f++)
    {
      links_inter_intra[f]=0;
      links_intra[f]=0;
    }


  for(f=1;f<=N;f++)       //recorro todos los nodos de la red
    {
      links_inter_intra[nodo_aux[f]]+=k[f];   //calculo conectividad inter+intra sumada 
    }                                  //a los nodos de cada comunidad
  

  for(f=1;f<=filas;f++)   //calculo la suma de links inter de cada comunidad
    {
      nodo1=D[f][1];
      nodo2=D[f][2];
      
      if(nodo_aux[nodo1] == nodo_aux[nodo2])    //si ambos pertenecen a la comunidad que estoy mirando
	{
	  s=nodo_aux[nodo1];
	  links_intra[s]++;
	}
    }

  
  auxiliar=filas;     //numero de links tot.
  Modularidad=0.0;
  for(s=1;s<=NumCom;s++)
    {
      corchete=0.0;

      corchete=links_inter_intra[s]/(2.0*auxiliar);
      corchete=corchete*corchete;
      corchete=-corchete;
      corchete=corchete+(links_intra[s]/auxiliar);

      Modularidad+=corchete;
    }

  
  

  return(Modularidad);


}


////////////////////////////
////////////////////////////////
//////////////////////////////////

double calculo_modularidad_anidado(int comun1)
{

  int s, f,nodo1,nodo2,num_links_anidado;
  double links_intra[N+1],auxiliar,links_inter_intra[N+1],corchete,Modularidad;

 

  for(f=1;f<=N;f++)
    {
      links_inter_intra[f]=0;
      links_intra[f]=0;
    }


  for(f=1;f<=size;f++)       //recorro solo los nodos de la comunidad escindida
    {
      nodo1=comunidad[f];

      links_inter_intra[nodo_aux_escin[nodo1]]+=k[nodo1];   //calculo conectividad inter+intra sumada 
    }                                  //a los nodos de las dos comunidades
  

  num_links_anidado=0;    //num. tot. links de la comunidad a escindir
  for(f=1;f<=filas;f++)   //calculo la suma de links inter de cada comunidad
    {
      nodo1=D[f][1];
      nodo2=D[f][2];

      if(nodo[nodo1]==comun1 || nodo[nodo2]==comun1)
	{   // sólo si alguno de los dos pertenece a la comunidad (OJO! INDICE ORIGINAL DE LA RED ENTERA!!)
	  num_links_anidado++;

	  if(nodo_aux_escin[nodo1] == nodo_aux_escin[nodo2])    
	    {          //si ambos pertenecen a la comunidad
	      s=nodo_aux_escin[nodo1];
	      links_intra[s]++;
	    }
	}
    }
  
  
  auxiliar=num_links_anidado;     //numero de links tot.
  Modularidad=0.0;
  for(s=1;s<=NumCom;s++)
    {
      corchete=0.0;

      corchete=links_inter_intra[s]/(2.0*auxiliar);
      corchete=corchete*corchete;
      corchete=-corchete;
      corchete=corchete+(links_intra[s]/auxiliar);

      Modularidad+=corchete;
    }

  
  

  return(Modularidad);


}
