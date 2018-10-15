//   20-11-09:  programa para implementar el modelo de Axelrod de consenso cultural
//     en lattices cuadradas con condiciones de contorno periodicas.
//   tb estudiaremos si hay diferencias entre la descripcion global y por capas



# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>



# define L   20       //lado de la red  (size = LxL)


# define kmed  4      //  numero de vecinos que tiene cada nodo (TODOS IGUAL!!)
//                      OJO! SI CAMBIO kmed, HABRIA QUE REHACER LA FUNCION QUE OBITIENE LA MATRIZ C[][] y D[][]!!


# define Niter 1    //estadistica






# define multiplo 10000     //cada cuantos pasos de tiempo calculo maginitudes y escribo archivos

# define Q  5    // numero de valores distintos que puede tomar un rasgo concreto (es igual para todos los rasgos)







# define f 10       // numero de rasgos de un nodo




///////////////////   Generacion de numeros aleatorios
 
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
//numero aleatorio flotante en el intervalo (0,1]
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];

//////////////////////////////////////////





int iter,C[L*L+1][kmed+1],k[L*L+1],D[2*L*L+1][3];
int x_d[L*L+1],x_i[L*L+1],y_ar[L*L+1],y_ab[L*L+1],N_links;

int v_i[L*L+1][f+1]; //matriz que guarda los valores de los f rasgos de los N nodos
double S_ij;
int tpo,contador, flag,Simil[2*L*L+1], NN[L*L+1],Phase[L*L+1][L*L+1],analyse[L*L+1],NCLUSTERS;
double GComp,S_max, n_links_blocked, n_links_S1,Sf_max_media,Sf_max[f+1];
int rasgo;
int Nclusters,NCL[f+1],NCLUSTERS;   // numero de clusters de consenso, de cada rasgo,su media y aux para la rutina
double NCL_med;

char file1[256],file2[256];

void inicia_rand(int semilla);
void crear_lattice();
void inicializar_lattice();
void matriz_D();
void consenso_link();
void contar_blocked_links();	
void marcar_link_consenso();
void find_Smax();          
void marcar_link_consenso_rasgo();




FILE *fich1,*fich2;


int main ()
{

    inicia_rand(time(0));


    //printf("\nAxelrod global consensus on square lattice. L=%d   Q=%d   f=%d  (%d iter)\n\n\n",L,Q,f,Niter);



    crear_lattice();    //esto solo lo hago una vez
    matriz_D();      



    for(iter=1;iter<=Niter;iter++)      //bucle estadistica
    {
	inicializar_lattice();    //doy valores a los f rasgos de los nodos
	

	
//Guarda S_max, <Sf_max> y Sf_max[rasgo] en cada paso de tiempo
	sprintf(file1,"Consenso_lattice_L%d_Q%d_f%d_iter%d.dat",L,Q,f,iter);
	fich1=fopen(file1,"wt");
	fclose(fich1);
	
// Ncluster de consenso total y por rasgos en cada paso de tiempo
	sprintf(file2,"NClusters_consenso_lattice_L%d_Q%d_f%d_iter%d.dat",L,Q,f,iter);
	fich2=fopen(file2,"wt");
	fclose(fich2);


	tpo=flag=contador=0;
	while(1)      //bucle de la evolucion temporal
	{
	    
	    consenso_link();    //en esta rutina elijo un link y con prob S_ij, n1 imita a n2		


	    if(((tpo % multiplo)==0.0) || (flag==1))   //cuento y escribo en los ficheros sólo cada multiplo veces 
	    {                                         //y tb la última vez cuando ya ha llegado al estado frozen

		contar_blocked_links();	  //y tb cuento los de Sij=1

		marcar_link_consenso();


		
		rasgo=0;  //para indicarle a la siguiente rutina que estoy calculando la GC de consenso global
	                  	// (pq qdo valga 1,...,f, estará calculando la GC por capas)		
		find_Smax();
		S_max=GComp;		
		Nclusters=NCLUSTERS;	


		for(rasgo=1;rasgo<=f;rasgo++)
		{
		    marcar_link_consenso_rasgo();
		    find_Smax();		    
		    Sf_max[rasgo]=GComp;
		    NCL[rasgo]=NCLUSTERS;

		}

		printf("S_max(t=%d)=%f  n_links_bl=%f  n_links_S1=%f  N_cl=%d \n\n",contador,S_max/(L*L),n_links_blocked,n_links_S1,NCLUSTERS);



// guardo la evolucion temporal de num links, Smax tot y por capas

		fich1=fopen(file1,"at");
		fprintf(fich1,"%d   %f   %f   %f",contador,n_links_blocked/N_links,n_links_S1/N_links,S_max/(L*L));
		for(rasgo=1;rasgo<=f;rasgo++)       
		{
		    fprintf(fich1,"   %f",Sf_max[rasgo]/(L*L));	
		}
		fprintf(fich1,"\n");
		fclose(fich1);	  
		


// guardo la evolucion temporal de num clusters tot y por capas
		fich2=fopen(file2,"at");
		fprintf(fich1,"%d   %d",contador,Nclusters);	
		for(rasgo=1;rasgo<=f;rasgo++)       
		{
		    fprintf(fich2,"   %d",NCL[rasgo]);	
		}
		fprintf(fich2,"\n");	
		fclose(fich2);	 
		





		contador++;	
	    }
	

	
	  if( flag==1)
	    { 	 						
		break;       // acabar la simulacion 
	    }	    	    	    
	    
	    

	    tpo++;    //pasos de tiempo reales de la simulacion 		   
	    
	}       //fin    bucle de la evolucion temporal
	
    }    //fin  bucle estadistica
    
    
}









////////////////////////////////////////
//////////////////////////////////////////
//////////////////////////////////////////


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
//////////////////////////////////////////
/////////////////////////////////////////



void  crear_lattice()
{
  int i,ii;

    for(i=1;i<=L*L;i++) //calculo lo que tengo que sumar/restar a uno nodo con indice i (en x o en y)  (metodo = Ising)
    {
	x_d[i]=1;     
	x_i[i]=-1;

	y_ar[i]=-L;    
	y_ab[i]=L;


    }
//con las excepciones de los bordes

    for(i=1;i<=L*L;i++)
    {
	if( (i%L)==1.0 ) //nodos en la primera columna
	{
	 x_i[i]=L-1;  
	 // printf("primera col: %d  \n",i);

	}

	if( (i%L)==0.0 )    //nodos de la ultima columna 
	{
	    x_d[i]=-L+1;
	    // printf("ultima  col: %d  \n",i);
	}

	if(i<=L)     //nodos de la primera fila
	{
	    y_ar[i]=L*(L-1);
	    // printf("primera fila: %d  \n",i);
	}

	if(i>L*(L-1))     //nodos de la ultima fila
	{
	    y_ab[i]=-L*(L-1);

	    // printf("ultima  fila: %d  \n",i);
	}
    }
    
    
    
    
   


    for(i=1;i<=L*L;i++)    //caso general. lleno la matriz de conectividad C[][] por orden:  dcha, izq, arriba, abajo.
    {
	C[i][1]=i+x_d[i];    //dcha
	C[i][2]=i+x_i[i];     //izq
	C[i][3]=i+y_ar[i];      //arriba
	C[i][4]=i+y_ab[i];      //abajo
	

	k[i]=kmed;
    }
    
  

    printf("*Vertices 400\n*Edges\n");
    for(i=1;i<=L*L;i++)
      {
	for(ii=1;ii<=4;ii++)
	  {
	    printf("%d  %d\n",i,C[i][ii]);
	  }
	if(i==200)
	  {
	    getchar();
	  }	
	
      }
    
   
    
    /* for(i=1;i<=L*L;i++)
    {
	printf("%d   ",i);
	if(( i% L)==0.0)
	{
	    printf("\n");
	}
	
	}*/




}




///////////////////////////////////////////
///////////////////////////////////////////
////////////////////////////////////////

void inicializar_lattice()
{
    int i,j;
    double v;
    
    for(i=1;i<=L*L;i++)
    {
	for(j=1;j<=f;j++)
	{
	    v=FRANDOM;
	    v_i[i][j]=v*Q;    //da valores entre 0 y Q-1 a cada rasgo
	}
    }
    
}



/////////////////////////////////////////////
//////////////////////////////////////////
/////////////////////////////////////////


void matriz_D()      // construyo la matriz de pares de vecinos a partir de la de conectividad
{
    int i,j,ii,indice,vecino,C_aux[L*L+1][kmed+1];


    for(i=1;i<=L*L;i++)     //copio la C[][] en una matriz aux para ir borrando los links, y no contarlos dos veces
    {
	for(j=1;j<=kmed;j++)
	{
	    C_aux[i][j]=C[i][j];
	}
	
    }

    for(i=1;i<=2*L*L;i++)     //copio la C[][] en una matriz aux para ir borrando los links, y no contarlos dos veces
    {	
	D[i][1]=0;
	D[i][2]=0;
    }
    

    indice=1;   // recuento del numero total de links
    for(i=1;i<=L*L;i++)     //OJO!!!!! NO CONTAR 4VECES CADA LINK!!!!
    {
	for(j=1;j<=kmed;j++)
	{
	    if(C_aux[i][j]!=0)   //si no lo he contado antes ya
	    {
		D[indice][1]=i;
		D[indice][2]=C[i][j];
		indice++;
		//	printf("%d-%d   ",i,C[i][j]);
		C_aux[i][j]=0;   //esto no importa en realidad pq no volvere a pasar por el

		vecino=C[i][j];
		for(ii=1;ii<=kmed;ii++)
		{
		    if(C_aux[vecino][ii]==i)
		    {
			C_aux[vecino][ii]=0;       //anulo el link del vecino j con el nodo i para no contarlo dos veces
			//	printf("%d-%d \n",vecino,C[vecino][ii]);
		    }
		}
	    }
	    // getchar();
	}			

    }

    N_links=indice-1 ;     //pq el ultimo ++ no cuenta
     
    // printf("n_links:%d",N_links);
//     getchar();

    /* for(i=1;i<=N_links;i++)
    {
	printf("%d-%d\n",D[i][1],D[i][2]);
	//getchar();
	}*/
    
    
    
}




/////////////////////////////////////////////
//////////////////////////////////////////
/////////////////////////////////////////



void consenso_link()       //rutina copiada del programa consenso_SF_AV_p_Cl_coef.c
{

    double aux_w, aux_r, aux;
    int i,w,r,ping;
    int nodo1,nodo2,n1,n2;
          
    
    aux_w=FRANDOM;     //elijo un link y guardo los nodos que une
    aux_w=aux_w*N_links;
    w=aux_w+1;   //por que genera numeros de 0 a  N_links -1

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
 
    /* printf("%d-%d  Sij=%f   ",n1,n2,S_ij);
       getchar();*/

    aux_w=FRANDOM; 
    if(aux_w<S_ij)    //prob de que n1 imite a n2 en un rasgo que no tengan igual
    {	
	ping=0;

	while(ping==0)
	{	    
	    if(S_ij!=0.0 && S_ij!=1.0 )
	    {    
		aux_r=FRANDOM;
		aux_r=aux_r*f; 
		r=aux_r+1;       //elijo un rasgo  (ojo!!!! la matriz se llena de 1 a f y el generador va de 0 a f-1)
	    
	
		if(v_i[n1][r] != v_i[n2][r])
		{		   
		    v_i[n1][r]=v_i[n2][r];  //n1 imita a n2		  
		    ping=1;

		    // printf("%d imita a %d en el rasgo %d\n",n1, n2,r);	
		}
		
	    }
	    else              //  PARA QUE NO SE QUEDE ATASCADO EN UN BUCLE QUE ESTE BLOQUEADO
	    {			    		
		ping=1;		    		    		
	    }
	}
    }
        
    
}




////////////////////////////////////////
///////////////////////////////////////////
//////////////////////////////////////////////


void contar_blocked_links()
{

    int n1, n2,i,j,similitud;
    


    n_links_blocked=0;
    n_links_S1=0;

    for(i=1;i<=N_links;i++)
    {
	n1=D[i][1];
	n2=D[i][2];

	similitud=0;
	for(j=1;j<=f;j++)
	{
	    if(v_i[n1][j]==v_i[n2][j])
	    {
		similitud++;
	    }
	}
	/*printf("%d-%d  Sij:%d   ",n1,n2,similitud);
	  getchar();*/

	if(similitud==f)
	{
	    n_links_blocked++;
	    n_links_S1++;
	}
	else
	    if(similitud==0)
	    {
		n_links_blocked++;
	    }
    }

    if(n_links_blocked==N_links)  //para acabar la simulacion;
    {
	flag=1;
    }

    //printf("N_links_blocked:%d   N_links_S1:%d\n",n_links_blocked,n_links_S1);
    // getchar();



    


}


/////////////////////////////////////////
//////////////////////////////
////////////////////////////



void marcar_link_consenso()
{
    int i,j,n1,n2;  
    
    
    for(i=1;i<=2*L*L;i++)   //inicializo el vector que guardará la similitud de cada link (en total, no por rasgos)
    {
	Simil[i]=0;
    }
    
    
    for(j=1;j<=L*L;j++)    
    {
	NN[j]=0;              //NN[] es k[] pero solo enlaces de consenso
	for(i=1;i<=kmed;i++)
	{
	    Phase[j][i]=0;  //Phase[][] es C[][] pero solo entran links de consenso
	}
    }


    for(j=1;j<=N_links;j++)  //recorro la matriz de links
    {
	n1=D[j][1];    
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


/*    for(i=1;i<=L*L;i++)
    {
	printf("%d: ",i);
	for(j=1;j<=kmed;j++)
	{
	    printf("%d  ",Phase[i][j]);
	}
	printf("\n ");
	getchar();
	
    }*/

   
}






/////////////////////////////////////////
//////////////////////////////
////////////////////////////


void find_Smax()
{

    int I,J,kk,jj,h,i,g,aislados,visitado[L*L+1];
    
    for(I=0;I<=L*L;I++)
    {
	visitado[I]=0;
    }

    
    aislados=0;
    GComp=0;

    NCLUSTERS=0;  
    for(I=1;I<=L*L;I++)     //recorro la lattice
    {
	kk=0;
	for(J=1;J<=L*L;J++)
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
			}
		    }
		}
	    }

	}
	else    //si el nodo esta aislado
	{
	    visitado[I]=1; 
	    aislados++;  
	}

	///////////////////////////////////////////y cuando acabo con el cluster actual:


	if(kk>GComp)      //guardo el mayor de los clusters como GComp del sistema y tb los nodos que lo forman
	{	    
	    GComp=kk;	   
	}



	if(kk>1)
	{
	    NCLUSTERS++;  
	}






	
    }   /////////////////////////fin del bucle sobre los nodos de la red



}






/////////////////////////////////////////
//////////////////////////////
////////////////////////////


void marcar_link_consenso_rasgo ()
{
    int i,j,n1,n2;  
    
    
    for(i=1;i<=N_links;i++)   //inicializo el vector que guardará la similitud de un link respecto al rasgo que este mirando
    {	
	Simil[i]=0;	      
    }
    
    
    for(j=1;j<=L*L;j++)   
    {
	NN[j]=0;                //NN[] es k[] pero solo enlaces de consenso
	for(i=1;i<=L*L;i++)
	{
	    Phase[j][i]=0;  //Phase[][] es C[][] pero solo entran nodos que comparten links de consenso
	}
    }        



//estamos mirando un rasgo concreto (rasgo =1,2,...f)

	for(j=1;j<=N_links;j++)  //recorro la matriz de links
	{
	    
	    n1=D[j][1];    
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
