#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define NODOS 4000     
#define N_INICIALES 2    
#define m_links 2      
#define kmax 400     
#define A 3   

#define niter 10
#define nb 1
#define ventana 1000
#define Transient 5000
#define JUGADAS 500000   
#define JUGADAS2 10000

#define R 1.0          //Reward
#define S 0.0          //Sucker's payoff
#define P 0.0          //Punishment


#define H 0.5          //fracc. inicial de cooperantes


//Generación de números aleatorios
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];

void init_rand(int seed);
void NetworkGrowth();
void Jugada();
void TopologiaCP();
void TopologiaDP();

int C[NODOS][kmax];
int k[NODOS];               
int kPA[NODOS];             
int P_PA[NODOS];                                      
double Pk[kmax],Pk2[kmax];            
int conectado[m_links];     

double GCMED,NCLUSTERSMED,KMAX,KMAXCLUST,KMAXCLUSTi,KMAXCLUSTERc,kmaxclusterc,
       CLUSTERS,CLUSTERSd,
       NKCI[kmax],NKCT[kmax],KCI[kmax],dk,dsize,NKci[kmax],NormaLinks,
       GCMEDd,NCLUSTERSMEDd,kmaxclusti,kmaxclustid,KMAXCLUSTid,max,GComp,
       GCompd,da,Tiempos[JUGADAS2+1],NKCId[kmax],NKCTd[kmax],KCId[kmax],
       NKcid[kmax];

int IMAX,NCLUSTERS,NCLUSTERSd,NN[NODOS],Phase[NODOS][kmax],warning[NODOS],
    analyse[NODOS],semilla[NODOS],iter,aant,TC[NODOS],numvar[NODOS];
int ESTADO[NODOS][2],NonNull,NonNullCP,NonNullDP,NonNullOSC;
int i,ib,j,s,ss,a,links,aleat_int,kk,ii,IIII,J,I,g,gg,ggg,aaa,a1,a2;
double aleatorio,links_totales,Puros,DPuros,Dpuros,puros,normat,FREC[NODOS],
    FRECMED,frecmed[kmax],PK[kmax],Pvar[kmax],PVAR[kmax],PCP[kmax],Pcp[kmax];
int links_PA, links_al,Validas;
int mlinks_PA[m_links],mlinks_al[m_links];
double payoff[NODOS];    //accumulated pay-offs  
double D,kmedia,frecfinal,frecPuros,frecDPuros,frecfinal2,coo,coop1,coop2,
       dobleprec;                             
double c_instant[ventana],Pendiente,C_valor,C_valor2,c_valor,c_valor2,
       cooperantes2,cooperantes,Threshold,dc;
int c,kmayor,d,Salir,contador,c_ant,size;
int coop[kmax],conect[kmax],Marcador[NODOS];
double P_clus[kmax],P_size[NODOS],RealizK[kmax],RealizKOsc[kmax]; 
double P_size2[NODOS];  //Distrib tamano comp gigante
double Slope,Slope1,Slope2,Slope3,Slope4,Slopef,NORMAFREC,normafrec[kmax];
int nod_clus[NODOS],k_clus[NODOS];

double T=1.75,delta_T=0.1,alfa=0.0;

FILE *probk,*probk3;
FILE *k_,*fich1;
FILE *fich7,*fich8,*fich9,*fich10,*fich11;
FILE *gig;



main()
{

char nombre [256],nombre2 [256],nombre3[256],file1[40],file7[40],
     file8[40],file10[256],file11[256],gigan[40],file9[256];

init_rand(8);

 size=NODOS;

//Guarda la Distribucion de Tiempos de Cooperacion 
 sprintf(file8,"%iTimes%.2lf-bini%.3lf.dat",size,alfa,T);
 fich8=fopen(file8,"wt");
 fclose(fich8);

//Guarda la Distribucion de Frecuencias versus K
 sprintf(file9,"%iFrecs%.2lf-bini%.3lf.dat",size,alfa,T);
 fich9=fopen(file9,"wt");
 fclose(fich9);
 
//Guarda la Transicion de Cooperantes y otras cosas
 sprintf(file7,"%iTransicion%.2lf-bini%.3lf.dat",size,alfa,T);
 fich7=fopen(file7,"wt");
 fclose(fich7);

//Guarda la Topologia de Cooperantes Puros
 sprintf(file10,"%iTopCP%.2lf-bini%.3lf.dat",size,alfa,T);
 fich10=fopen(file10,"wt");
 fclose(fich10);

//Guarda la Topologia de Defectores Puros
 sprintf(file11,"%iTopDP%.2lf-bini%.3lf.dat",size,alfa,T);
 fich11=fopen(file11,"wt");
 fclose(fich11);


 Threshold=10./((double)ventana);


 
for(ib=0;ib<nb;ib++)
{
// INICIO DE VARIABLES GLOBALES

    NonNull=0;
    NonNullCP=0;
    NonNullDP=0;
    NonNullOSC=0;

    cooperantes=0.;
    cooperantes2=0.;
    Validas=0;
    Puros=0.;
    DPuros=0.;

    GCMED=0.;
    NCLUSTERSMED=0.;
    KMAXCLUSTi=0;
    KMAXCLUSTERc=0;
    GCMEDd=0.;
    NCLUSTERSMEDd=0.;
    KMAXCLUSTid=0;

    for(J=0;J<kmax;J++)
    {
	NKCT[J]=0.;
	NKCI[J]=0.;
	KCI[J]=0.;
	
	NKCTd[J]=0.;
	NKCId[J]=0.;
	KCId[J]=0.;
	RealizK[J]=0.;
	//P_clus[J]=0.;
	//P_size[J]=0.;
	//P_size2[J]=0.;
    }   

    for(i=0;i<=JUGADAS2;i++)
    {
	Tiempos[i]=0.; //Guarda la distribucion de tiempos de Coop.
    }

    normat=0.; //Norma para la distribucion de tiempos de Coop.
  
    for(i=0;i<kmax;i++)
    {
	coop[i]=0;
	conect[i]=0;
	Pk2[i]=0.0;
	PK[i]=0.;
	PVAR[i]=0.;
	PCP[i]=0.;
    }
  
    FRECMED=0.; //Frecuencia media de cambio.
    NORMAFREC=0.; //Numero de tios que han cambiado.

    for(i=0;i<kmax;i++)
    {
	frecmed[i]=0.;  //Distrib. de Frecs. en funcion de K.
	normafrec[i]=0.;  //Numero de oscilantes de conectividad K
    }

   
for(iter=0;iter<niter;iter++)
{  
    printf("%.3lf ************* %i\n",T,iter);

// INICIA VERIABLES LOCALES

    for(i=0;i<NODOS;i++)
    {
	ESTADO[i][0]=0;  //Guarda el Estado del nodo (C=1,D=0).
	ESTADO[i][1]=0;  //Guarda el tiempo en que i paso D->C.
	TC[i]=0;         //Guarda el Tiempo de Coop del nodo i.
	FREC[i]=0.;      //Guarda la Frec media del nodo i. 
	numvar[i]=0;     //Numero de veces que cambia i.
    }

    for(i=0;i<kmax;i++)
    {
	Pvar[i]=0.;
	Pcp[i]=0.;
	Pk[i]=0.;
    }

    kmaxclusti=0.;
	
    NetworkGrowth();
    

    //Inicialización de los dos tipos: cooperante y defector (C=0,D=1).
    for(i=0;i<NODOS;i++)
    {
	if(FRANDOM<H)
	    kPA[i]=0;
	else 
	    kPA[i]=1;    
    }
    
    contador=0;
    C_valor=0.;
    C_valor2=0.;

    //Variables para minimos cuadrados (PARADA):
    Slope=0.;
    Slope1=0.;
    Slope2=0.;
    Slope3=0.;
    Slope4=0.;


    //JUGAMOS (TRANSITORIO Y VENTANAS DE PARADA):
   
    for(a=0;a<JUGADAS;a++)
    {
	Jugada();
	da=a;
	
        // Si hemos pasado el transitorio	
	if(a>Transient)
	{
	    contador++;
	    c_instant[contador-1]=c-c_ant;
	    C_valor=C_valor+c;
	    dc=c;	   
	    Slope1=Slope1+dc*da;
	    Slope2=Slope2+da;
	    Slope3=Slope3+dc;
	    Slope4=Slope4+da*da;
	    
	    if(contador==ventana)
	    {
		contador=0;
		Pendiente=0;
		dobleprec=ventana;	
		Slope=dobleprec*Slope1-Slope2*Slope3;
		Slope=Slope/(dobleprec*Slope4-Slope2*Slope2);
		
		dobleprec=ventana;
		c_valor=C_valor/dobleprec;
		
		for(i=0;i<ventana;i++)
		{
		    Pendiente=Pendiente+fabs(c_instant[i]);
		}
		dobleprec=ventana;
		Pendiente=Pendiente/dobleprec;
		
		for(i=0;i<ventana;i++)
		{
		    c_instant[i]=0;
		}
		
		Salir=0;
		
		//if(fabs(Pendiente)<=Threshold)
		//{
		//printf("Jugadas:%i  %lf  %lf\n",
                //a,fabs(Pendiente),Threshold);
		//aant=a;  
		//a=JUGADAS;
		//  Salir=1;
		//}
		if(fabs(Slope)<=Threshold)
		{
		    aant=a;  
		    a=JUGADAS;
		    Salir=1;
		    Slopef=Slope;
		}		
		else
		{
		    printf("NO--- %lf\n",Slope);
		    Slopef=Slope;
		    C_valor=0.;
		    Slope=0.;
		    Slope1=0.;
		    Slope2=0.;
		    Slope3=0.;
		    Slope4=0.;
		    aant=a;
		}
	    }
	}
	
	c_ant=c;
	
    }        //Terminó la evolución del juego
    
    c_ant=c;
    
    
   
    
    if(Salir==1)
    {
	for(i=0;i<kmax;i++)
	{
	    for(c=0;c<NODOS;c++)
	    {
		if(k[c]==i)
		{
		    conect[i]++; // KTOP
		    if(kPA[c]==0)
			coop[i]++; //KTOPc
		}
	    }    
	}
    }
    
    c=0;
    s=0;
    
    printf("Iter  %i   Valido: %i - %i, Slope %.4lf <%.4lf>\n",
           iter+1,Salir,aant,Slopef,c_valor);
    
    if(c_valor>0.)
    {NonNull++;}
    
    cooperantes+=c_valor;
    
    //calculamos la frecuencia final de los cooperadores
    
    for(i=0;i<NODOS;i++)     
    {
	Marcador[i]=0;     //todos como fluct
	if(kPA[i]==0)     //si es coop
	{c++;
	Marcador[i]=1;   //es coop
	    ESTADO[i][0]=1;
	    ESTADO[i][1]=0;
	    }
	else{Marcador[i]=-1;}    //si es defect, defect

	s+=k[i];            
    }

 
    if(Salir==1)
    {
    	Validas++;
    }

    if(c_valor>0.)
    {

//JUGAMOS (SEGUNDA VENTANA GRANDE):

	for(a=0;a<JUGADAS2;a++)
	{
	    Jugada();
	
	    for(i=0;i<NODOS;i++)     
	    {
		if(kPA[i]==ESTADO[i][0]) //si ha cambiado
		{
		    if(ESTADO[i][0]==0)  //y si era defect
		    {
			ESTADO[i][0]=1;  //ahora es coop
			ESTADO[i][1]=a;
		    }
		    else       //y si era coop
		    {
			ESTADO[i][0]=0;     //ahora es defect
			numvar[i]++;
			FREC[i]=FREC[i]+(a-ESTADO[i][1]);
			ESTADO[i][1]=0;
		    }
		}
	
		if(kPA[i]==0)    //si es coop
		{
		    C_valor2=C_valor2+1.;
		    TC[i]=TC[i]+1;
		}
	    
		if(Marcador[i]==1)  //si era coop puro
		{
		    if(kPA[i]==1)    //y ahora ha sido defect
		    {
			//Pvar[k[i]]++;
			Marcador[i]=0;    //fluctuante pa siempre
		    } 
		}

		if(Marcador[i]==-1)  //si era defect puro
		{
		    if(kPA[i]==0)     //y ahora ha sido coop
		    {
			Marcador[i]=0;      //fluctuante pa siempre
			//Pvar[k[i]]++;
		    }
		}
		
	    }	
	}    //FIN VENTANA GRANDE (DISCRIMINAR)
    

	puros=0.;
	Dpuros=0.;

	for(i=0;i<NODOS;i++)     
	{
	    if(Marcador[i]==1)
	    {
		puros++;
		Pcp[k[i]]++;
	    }
	    
	    if(Marcador[i]==-1)
	    {
		Dpuros++;
	    }
	}

	if(Dpuros+puros<NODOS)
	{NonNullOSC++;}

	dobleprec=JUGADAS2;
	c_valor2=C_valor2/dobleprec;
	cooperantes2=cooperantes2+c_valor2;
	Puros=Puros+puros;
	DPuros=DPuros+Dpuros;

	for(i=0;i<NODOS;i++)     
	{

	    if(TC[i]>0)
	    {
		normat++;
		Tiempos[TC[i]]=Tiempos[TC[i]]+1.;
	    }
	}

	for(i=0;i<NODOS;i++)
	{
	    if(numvar[i]>0)
	    {
		dobleprec=numvar[i];
		FRECMED=FRECMED+FREC[i]/dobleprec;
		NORMAFREC++;
		frecmed[k[i]]=frecmed[k[i]]+FREC[i]/dobleprec;
		normafrec[k[i]]++;
		Pvar[k[i]]++;
	    }
	}

	TopologiaCP();
	TopologiaDP();



	for(i=0;i<kmax;i++)
	{
	    if(Pk[i]>0.)
	    {
		if(Dpuros+puros<NODOS)
		{
		    RealizKOsc[i]++;
		    Pvar[i]=Pvar[i]/Pk[i];
		}
		if(NCLUSTERS>0)
		{
		    RealizK[i]++;
		    Pcp[i]=Pcp[i]/Pk[i];
		}
	    }

	    else
	    {
		Pvar[i]=0;
		Pcp[i]=0;
	    }

	    PVAR[i]=PVAR[i]+Pvar[i];
	    PCP[i]=PCP[i]+Pcp[i];
	}

       

    }
    else
    {
      	c_valor2=0.;
	cooperantes2=cooperantes2+0.;
	Puros=Puros+0.;
	DPuros=DPuros+NODOS;
    }
}        //fin bucle estadística

// PROMEDIAMOS

 FRECMED=FRECMED/NORMAFREC;
 
 for(i=0;i<kmax;i++)
 {
     frecmed[i]=frecmed[i]/normafrec[i];
     dobleprec=niter;    
     PK[i]=PK[i]/dobleprec;

     dobleprec=RealizKOsc[i];

     if(RealizKOsc[i]>0)
     {
	 PVAR[i]=PVAR[i]/dobleprec;
     }
     else{PVAR[i]=-10000.0;}

     dobleprec=RealizK[i]; //OJO porque no es NonNullCP

     if(RealizK[i]>0)
     {     
	 PCP[i]=PCP[i]/dobleprec;
     }
     else{PCP[i]=-10000.;}
 }

 fich8=fopen(file8,"at");
 fich9=fopen(file9,"at");

 dobleprec=niter*NODOS;

 frecfinal=cooperantes/dobleprec;
 frecfinal2=cooperantes2/dobleprec;
 
 frecPuros=Puros/dobleprec;
 frecDPuros=DPuros/dobleprec;


 for(i=0;i<kmax;i++)     
 {
     fprintf(fich9,"%lf %i  %lf  %lf  %lf  %lf  %i  %i\n",
             T,i,frecmed[i],PVAR[i],PCP[i],PK[i],Validas,NonNull);
 }
 fprintf(fich9,"\n");
 fclose(fich9);

 for(i=1;i<=JUGADAS2;i++)     
 {
     Tiempos[i]=Tiempos[i]/normat;
     fprintf(fich8,"%lf   %i   %lf\n",T,i,Tiempos[i]); 
 }
 fprintf(fich8,"\n");
 fclose(fich8);	           


///////// Cosas de Topologia

 dobleprec=NonNullCP;
 
 if(dobleprec>0)
 {
     GCMED=GCMED/dobleprec;
     NCLUSTERSMED=NCLUSTERSMED/dobleprec;
     KMAXCLUSTi=KMAXCLUSTi/dobleprec;
 }
 else
 {
     GCMED=0.;
     NCLUSTERSMED=0.;
     KMAXCLUSTi=0.;
 }

 dobleprec=NonNullDP;

 if(dobleprec>0)
 {
     GCMEDd=GCMEDd/dobleprec;
     NCLUSTERSMEDd=NCLUSTERSMEDd/dobleprec;
     KMAXCLUSTid=KMAXCLUSTid/dobleprec;
 }
 else
 {
     GCMEDd=0.;
     NCLUSTERSMEDd=0.;
     KMAXCLUSTid=0.;
 }

 fich7=fopen(file7,"at");

 fprintf(fich7,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %i  %i  %i  %i  %i  %lf\n",
	 T,frecfinal,frecfinal2,frecPuros,frecDPuros,
	 FRECMED,GCMED,GCMEDd,NCLUSTERSMED,NCLUSTERSMEDd,KMAXCLUSTi,
	 KMAXCLUSTid,Validas,NonNull,NonNullCP,NonNullDP,NonNullOSC,NORMAFREC);


printf("b:%lf  alfa:%lf   GCMED:%lf  GCMEDd:%lf  NCLUSTERSMED:%lf  NCLUSTERSMEDd:%lf    NonNullCP:%i  NonNullDP:%i  \n",
	 
	 T,alfa,GCMED,GCMEDd,NCLUSTERSMED,NCLUSTERSMEDd,NonNullCP,NonNullDP);



 fclose(fich7);

 fich10=fopen(file10,"at");
 dobleprec=NonNullCP;

 for(J=0;J<kmax;J++)
 { 
     if(NKCT[J]>0){KCI[J]=KCI[J]/NKCT[J];}
     else{KCI[J]=0.;}

     if(NonNullCP>0){NKCI[J]=NKCI[J]/dobleprec;} 
     else{NKCI[J]=0.;}
     
     fprintf(fich10,"%lf  %i  %lf  %lf\n",T,J,KCI[J],NKCI[J]); 
 }
 fprintf(fich10,"\n");
 fclose(fich10);


 fich11=fopen(file11,"at");
 dobleprec=NonNullDP;

 for(J=0;J<kmax;J++)
 { 
     if(NKCTd[J]>0){KCId[J]=KCId[J]/NKCTd[J];}
     else{KCId[J]=0.;}

     if(NonNullDP>0){NKCId[J]=NKCId[J]/dobleprec;} 
     else{NKCId[J]=0.;}
     
     fprintf(fich11,"%lf  %i  %lf  %lf\n",T,J,KCId[J],NKCId[J]); 
 }
 fprintf(fich11,"\n"); 
 fclose(fich11);

 T+=delta_T;

}

}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////



void init_rand(int seed)
{
int ial;
int dummy;

srand((unsigned)seed);
for(ial=0;ial<111;ial++)
  rand();
ip=128;
ip1=ip-24;
ip2=ip-55;
ip3=ip-61;
for(ial=0;ial<256;ial++)
  ira[ial]=(unsigned)rand()+(unsigned)rand();
for(ial=0;ial<1111;ial++)
  dummy=RANDOM;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

void NetworkGrowth()
{
   
    do
    {
      
	gg=1;
links=0;
for(i=0;i<N_INICIALES;i++)
      links=links+i;

links_totales=links+(NODOS-N_INICIALES)*m_links;



//Inicializamos

k[0]=N_INICIALES-1;
kPA[0]=N_INICIALES-1;
P_PA[0]=kPA[0]+A;                          

//nodos del núcleo inicial
 for(i=1;i<N_INICIALES;i++)          
 {
     k[i]=N_INICIALES-1;      
     kPA[i]=N_INICIALES-1;
     P_PA[i]=P_PA[i-1]+kPA[i]+A;     
 }
 
 for(i=N_INICIALES;i<NODOS;i++) 
 {
     k[i]=0;
     kPA[i]=0;      
     P_PA[i]=P_PA[N_INICIALES-1];
 }

 for(i=0;i<NODOS;i++)
     for(j=0;j<kmax;j++)
	 C[i][j]=-2;
 
 for(i=0;i<N_INICIALES;i++)
 {                                               
     for(j=0;j<i,j<k[i];j++)
	 C[i][j]=j;
     for(j=i+1;(j-1)<k[i];j++)
	 C[i][j-1]=j;
 }

//Comienza el crecimiento de la red

 for(i=N_INICIALES;i<NODOS;i++) 
 {                               
     links_PA=0;
     links_al=0; 
                                
     for(j=0;j<m_links;j++)    
     {         
         if(FRANDOM<alfa)       
	 {
	     aleatorio=FRANDOM*NODOS;
	     aleat_int=(int)aleatorio;                        
	     conectado[j]=aleat_int;
	     
	     for(aaa=0;aaa<kmax;aaa++)
	     {
		 if(C[i][aaa]==-2)                           
		 {
		     for(a1=0;a1<aaa;a1++)
			 if (conectado[j]==C[i][a1] || conectado[j]==i)
			 {
			     aleatorio=FRANDOM*NODOS;
			     aleat_int=(int)aleatorio;      
			     conectado[j]=aleat_int;
			 }
		     aaa=kmax; 
		 }                        
	     }
	     for(aaa=0;aaa<j;aaa++)
	     {
		 if (conectado[j]==conectado[aaa] || conectado[j]==i)
		 {
		     aleatorio=FRANDOM*NODOS;
		     aleat_int=(int)aleatorio;                  
		     conectado[j]=aleat_int;
		 }
	     }
	     
	     for(aaa=0;aaa<j;aaa++)
	     {
		 if (conectado[j]==conectado[aaa] || conectado[j]==i)
		 {                                        
		     conectado[j]=i-1;
		     for(a1=0;a1<j;a1++)
		     {
			 if (conectado[a1]==conectado[j] || conectado[j]==i)
			     conectado[j]=i-2;
		     }
		 }
	     }                               
	     
	     if(kPA[conectado[j]]==0)
	     {                     
		 mlinks_al[links_al]=conectado[j];
		 links_al++;                     
	     }
	 }
         else                    //bucle (3b) en el artículo
	 {                                 
	     aleatorio=FRANDOM*P_PA[NODOS-1];
	     for(s=0;s<NODOS;s++)
	     {
		 if(P_PA[s]>aleatorio)
		 {
		     conectado[j]=s;
		     a2=s;
		     s=NODOS;
		 }
	     }
	     
	     for(aaa=0;aaa<kmax;aaa++)
	     {
		 if(C[i][aaa]==-2)                           
		 {
		     for(a1=0;a1<aaa;a1++)
			 if (conectado[j]==C[i][a1] || conectado[j]==i)
			 {
			     aleatorio=FRANDOM*P_PA[NODOS-1];
			     for(s=0;s<NODOS;s++)
			     {
				 if(P_PA[s]>aleatorio)
				 {
				     conectado[j]=s;
				     a2=s;
				     s=NODOS;
				 }
			     }
			 }
		     aaa=kmax; 
		 }                        
	     }
	     for(aaa=0;aaa<j;aaa++)
	     {
		 if (a2==conectado[aaa])
		 {
		     aleatorio=FRANDOM*P_PA[NODOS-1];
		     for(s=0;s<NODOS;s++)
		     {
			 if(P_PA[s]>aleatorio)
			 {
			     conectado[j]=s;
			     a2=s;
			     s=NODOS;
			 }
		     }
		 }
	     }     
	     
	     for(s=conectado[j];s<NODOS;s++)
	     {                      
		 P_PA[s]=P_PA[s]-kPA[conectado[j]]-A;                      
	     }
	     mlinks_PA[links_PA]=conectado[j];
	     links_PA++;
                 
	 }
     }
     
     for(j=0;j<links_PA;j++)
     {
         ss=mlinks_PA[j];                    
         kPA[ss]++;         
         for(s=ss;s<NODOS;s++)
	 {
	     P_PA[s]=P_PA[s]+kPA[ss]+A;
	 }
     }
     
     for(j=0;j<links_al;j++)
     {
         ss=mlinks_al[j];
         for(s=ss;s<NODOS;s++)
	 {
	     P_PA[s]=P_PA[s]+A;
	 }
     }
      
     if(k[i]==0)
     {
         for(s=i;s<NODOS;s++)
	 {
	     P_PA[s]=P_PA[s]+A;
	 }   
     }
     
     for(j=0;j<kmax;j++)//Nuevos links a la matriz C
     {
         if(C[i][j]==-2)
	 {
	     for(s=0;s<m_links;s++)
	     {
		 C[i][j+s]=conectado[s]; 
	     }
	     j=kmax;
	 }
     }
       
     for(s=0;s<m_links;s++)
     {
         ss=conectado[s];
         C[i][k[i]+s]=ss;
     }
      
     for(s=0;s<m_links;s++)
     {
         ss=conectado[s];
         C[ss][k[ss]]=i;
         k[ss]++;
     }
     
     k[i]+=m_links;
          
 }
 
 for(i=0;i<NODOS;i++)
     for(j=0;j<k[i];j++)
     {
	 if(C[i][j]==i || C[i][j]<0)
	     gg=0;
     }

    }while (gg==0);


     
//Ha acabado la formación de la red

for(i=0;i<kmax;i++)     //Construimos el vector de P(k)
 {
     Pk[i]=0.0;
     for(j=0;j<NODOS;j++)
     {
	 if(i==k[j])
	 {
	     Pk[i]++;
	 }
     }
     Pk[i]=Pk[i]/((double)NODOS);
     PK[i]=PK[i]+Pk[i]; 
 }
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

void Jugada()
{
   ss=0;
   for(i=0;i<NODOS;i++)
         {
         payoff[i]=0.0;
         if(kPA[i]==0)
             ss++;
         }
   for(i=0;i<NODOS;i++)//Calculamos el "accumulated payoff" 
   {
       if(kPA[i]==0) //Si el nodo i es cooperador
       {
	   for(j=0;j<k[i];j++)
	   {
	       if(kPA[C[i][j]]==0)
	       {
                   payoff[i]+=R;                   
	       }
	       else
	       {
                   payoff[i]+=S;                   
	       }
	   }
       }
       else    //Si el nodo i es defector
       {
	   for(j=0;j<k[i];j++)
	   {
	       if(kPA[C[i][j]]==0)
	       {
                   payoff[i]+=T;                   
	       }
	       else
	       {
                   payoff[i]+=P;                   
	       }
	   }
       } 
   }

   for(i=0;i<NODOS;i++) //Cambios cooperador<->defector 
   {
       aleatorio=FRANDOM*k[i];
       c=(int)aleatorio;
       c=C[i][c];
       if(c==-2)
       {printf("ERROR-----------------\n");}

       if(payoff[i]<payoff[c])
       {
           D=T-S;
           if(k[i]>k[c])  
	       kmayor=k[i];
           else
	       kmayor=k[c];
           if(FRANDOM<((payoff[c]-payoff[i])/(D*kmayor)))
	   {
	       kPA[i]=kPA[c];             
	   }           
       }         
   }

   c=0;
   for(i=0;i<NODOS;i++)    
   {
       if(kPA[i]==0)
	   c++;
   }

}


void TopologiaCP()
{
    
    for(J=0;J<NODOS;J++)
    {
	NN[J]=0;  //NN[] es k[] pero solo enlaces entre cooperantes 
	for(IIII=0;IIII<kmax;IIII++)
	{
	    Phase[J][IIII]=0;  //Phase[][] es C[][] pero sólo entre cooperantes
	}
    }

    for(J=0;J<NODOS;J++)
    {
	if(Marcador[J]==1) //Si es Coop. Puro
	{	
	    for(IIII=0;IIII<k[J];IIII++)
	    {
		if(C[J][IIII]>J)
		{
		    if(Marcador[C[J][IIII]]==1)
		    {
			NN[J]++;
			NN[C[J][IIII]]++;
			Phase[J][NN[J]-1]=C[J][IIII];
			Phase[C[J][IIII]][NN[C[J][IIII]]-1]=J;
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

    IMAX=0;
    
    for(I=0;I<NODOS;I++)
    {warning[I]=0;}   //Indica que se ha mirado si el nodo está en un cluster 
                      //de cooperadores con un 1 y que no se ha mirado con un 0
    GComp=0.;
    NCLUSTERS=0;

    //////

    aaa=0;
    for(I=0;I<NODOS;I++)
    {
	kk=0;

	for(J=0;J<NODOS;J++)
	{analyse[J]=0;}              //Indica los nodos que forman parte del 
	                             //cluster que se mira en ese momento

	if(NN[I]>0)
	{
	    if(warning[I]==0)
	    {
		semilla[NCLUSTERS]=I; //Indica el nodo menor de cada cluster
		kk=0;
		warning[I]=1;
		analyse[kk]=I;
		for(J=0;J<=kk;J++)
		{ 
		    g=analyse[J]; 
		    for(ii=0;ii<NN[g];ii++)
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

	if(kk!=0)
	{
	    max=kk;
	    NCLUSTERS++;
	    CLUSTERS++;
	    max=max;
	    if(GComp<max)
	    {
		GComp=max;
		aaa=kk;
		IMAX=NCLUSTERS-1; //Sera indice del vector semilla para la GC
	    }
	    //P_size[kk]+=1.0;
	}	
    }
   
    //P_clus[NCLUSTERS]+=1.0;   
    //P_size2[aaa]+=1.0;
	    
    
    //printf("GC: %lf\n NClusters: %i\n",GComp*(float)NODOS,NCLUSTERS);
           
    if(NCLUSTERS>0)
    {
	NonNullCP++;
	NCLUSTERSMED=NCLUSTERSMED+NCLUSTERS;

// ANALISIS DE LA COMPONENTE GIGANTE
    
    //kmaxclust=0;

    kmaxclusti=10000; // Para la Kmin del cluster
    kmaxclusterc=0;

    for(I=0;I<NODOS;I++)
    {warning[I]=0;}

    for(J=0;J<NODOS;J++)
	{
	    analyse[J]=0;
	}

    I=semilla[IMAX];
    kk=0;
    warning[I]=1;
    analyse[kk]=I;
    
    NormaLinks=0.;
    for(J=0;J<kmax;J++)
    { 
	NKci[J]=0.;
    }

    for(J=0;J<=kk;J++)
    { 

	g=analyse[J];
	ggg=k[g];
	NKCT[ggg]++;  	//Guarda número de nodos del cluster en función de 
                        //la conectividad k[g]

	ggg=k[g];
	KCI[ggg]=KCI[ggg]+NN[g]; //Guarda la conectividad entre cooperadores 
                                 //de los nodos del cluster 
				 //en función de la conectividad k[g]
	
	ggg=NN[g];
	NKci[ggg]++;	//Guarda el número de nodos del cluster en función de 
                        //la conectividad entre cooperadores NN[g] 
	NormaLinks++;

	for(ii=0;ii<NN[g];ii++)
	{ 	    
	    gg=Phase[g][ii];	    
	    if(warning[gg]==0)
	    {
		warning[gg]=1;
		kk=kk+1;
		analyse[kk]=gg;
	    }			
	}

	if(NN[g]>1)
	{
	    if(NN[g]<kmaxclusti)
	    {
		kmaxclusti=NN[g];
	    }	
	}
    }
    
    GCMED=GCMED+GComp;
    KMAXCLUSTi=KMAXCLUSTi+kmaxclusti;
	

    for(J=0;J<kmax;J++)
    { 
	NKCI[J]=NKCI[I]+NKci[J]/NormaLinks;
    }

}
// TRAS LAS ITERS!!!!!!!!!!!!!!!!

// PROMEDIAMOS
}


/////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


void TopologiaDP()
{
    
    for(J=0;J<NODOS;J++)
    {
	NN[J]=0;  //NN[] es k[] pero solo enlaces entre defectores 
	for(IIII=0;IIII<kmax;IIII++)
	{
	    Phase[J][IIII]=0;  //Phase[][] es C[][] pero sólo entre defectores
	}
    }

    for(J=0;J<NODOS;J++)
    {
	if(Marcador[J]==-1) //Si es Def. Puro
	{	
	    for(IIII=0;IIII<k[J];IIII++)
	    {
		if(C[J][IIII]>J)
		{
		    if(Marcador[C[J][IIII]]==-1)
		    {
			NN[J]++;
			NN[C[J][IIII]]++;
			Phase[J][NN[J]-1]=C[J][IIII];
			Phase[C[J][IIII]][NN[C[J][IIII]]-1]=J;
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

    IMAX=0;
    
    for(I=0;I<NODOS;I++)
    {warning[I]=0;}   //Indica que se ha mirado si el nodo está en un cluster 
                      //de defectores con un 1 y que no se ha mirado con un 0
    GCompd=0.;
    NCLUSTERSd=0;

    //////

    aaa=0;

    for(I=0;I<NODOS;I++)
    {
	kk=0;
	for(J=0;J<NODOS;J++)
	{analyse[J]=0;}              //Indica los nodos que forman parte del 
	                             //cluster que se mira en ese momento

	if(NN[I]>0)
	{
	    if(warning[I]==0)
	    {
		semilla[NCLUSTERS]=I; //Indica el nodo menor de cada cluster
		kk=0;
		warning[I]=1;
		analyse[kk]=I;
		for(J=0;J<=kk;J++)
		{ 
		    g=analyse[J]; 
		    for(ii=0;ii<NN[g];ii++)
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

	if(kk!=0)
	{
	    max=kk;
	    NCLUSTERSd++;
	    CLUSTERSd++;
	    max=max;
	    if(GCompd<max)
	    {
		GCompd=max;
		aaa=kk;
		IMAX=NCLUSTERSd-1; //Sera indice del vector semilla para la GC
	    }
	    //P_size[kk]+=1.0;
	}	
    }
   
    //P_clus[NCLUSTERS]+=1.0;   
    //P_size2[aaa]+=1.0;
	    
    
    //printf("GC: %lf\n NClusters: %i\n",GComp*(float)NODOS,NCLUSTERS);
           
    if(NCLUSTERSd>0)
    {
	NonNullDP++;
	
	NCLUSTERSMEDd=NCLUSTERSMEDd+NCLUSTERSd;
    

// ANALISIS DE LA COMPONENTE GIGANTE
    
	//kmaxclust=0;
	
	kmaxclustid=10000; // Para la Kmin del cluster
//	kmaxclusterc=0;

	for(I=0;I<NODOS;I++)
	{warning[I]=0;}

	for(J=0;J<NODOS;J++)
	    {
		analyse[J]=0;
	    }

	I=semilla[IMAX];
	kk=0;
	warning[I]=1;
	analyse[kk]=I;
    
	NormaLinks=0.;
	for(J=0;J<kmax;J++)
	{ 
	    NKcid[J]=0.;
	}

	for(J=0;J<=kk;J++)
	{ 
	    
	    g=analyse[J];
	    ggg=k[g];
	    NKCTd[ggg]++;  	//Guarda número de nodos del cluster en función de 
	                        //la conectividad k[g]

	    ggg=k[g];
	    KCId[ggg]=KCId[ggg]+NN[g]; //Guarda la conectividad entre defectores 
	                               //de los nodos del cluster 
				       //en función de la conectividad k[g]
	
	    ggg=NN[g];
	    NKcid[ggg]++;	//Guarda el número de nodos del cluster en función de 
                                //la conectividad entre defectores NN[g] 
	    NormaLinks++;

	    for(ii=0;ii<NN[g];ii++)
	    { 	    
		gg=Phase[g][ii];	    
		if(warning[gg]==0)
		{
		    warning[gg]=1;
		    kk=kk+1;
		    analyse[kk]=gg;
		}			
	    }

	    if(NN[g]<kmaxclustid)
	    {
		kmaxclustid=NN[g];
	    }	
	}
    
	GCMEDd=GCMEDd+GCompd;
	KMAXCLUSTid=KMAXCLUSTid+kmaxclustid;
	

	for(J=0;J<kmax;J++)
	{ 
	    NKCId[J]=NKCId[I]+NKcid[J]/NormaLinks;
	}

    }
// TRAS LAS ITERS!!!!!!!!!!!!!!!!

// PROMEDIAMOS
}
