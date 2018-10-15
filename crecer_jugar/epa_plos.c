/////////////////////////////////////////////////////////////////////////////////////////
//     Source code for generating complex via Evolutionary Preferential Attachment     //
/////////////////////////////////////////////////////////////////////////////////////////

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# define N 4000 // Number of nodes of the network
# define m_o  2 // Number of nodes contained in the initial (fully-connected) core
# define m   2  // Number of links launched by each node (m<=m_o)
# define K_MAX 4000 //Overstimated number of the maximum degree a node can reach

///////////////////////////////////////////////////////////
// Variables for Random numbers generation: 
//
#define FNORM (2.3283063671E-10F)
#define RANDOM ((ira[ip++]=ira[ip1++]+ira[ip2++])^ira[ip3++])
#define FRANDOM (FNORM*RANDOM)
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
double r,v; 
//
///////////////////////////////////////////////////////////

// Some auxiliary variables and counters:
int  II,i,j,jj,g,gg,y,q,steps,t,w,counter; 

///////////////////////////////////////////////////////////
// Variables for Network topology: 
//
int k[N];  //Degree of nodes      
int  C[N][K_MAX]; //For each node, it stores the nodes it is attached to
int unido[m];  //List of nodes attached to by the last node incorporated
double P[N], P_prov[N];  //Probabilities of attachment
double norma,norma2; //Norms for the probabilities of attachment
int s; //Intant size of the network
double fit[N]; //Fitness of nodes (entering in the attachment probability)
double eps; //Weight of the node's fitness into the attachment probability (i.e. the degree of selection).
//
///////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////
//Variables for Evolutionary Dynamics:
//
double b,R,Pu,Su; //Variables of the Dilemma: b (Temptation), R (Reward), Pu (Punishment), Su (Suckers payoff)
int  e[N]; //State (strategy) of each node (we use "0" for cooperators and "1" for defectors)
int  e_a[N]; //Same as e[N], used for updating the strategy in parallel
int k_m;  //It will allocate the maximum degree of the neighbors of a node
double  ben[N]; //Payoff of the nodes after playing a round robin
double p; //It will allocate the probability that a node changes it  strategy 
double ro; // Probability that a new node plays as cooperator in its first game
int  n_coop; //Number of cooperators in the network
//
///////////////////////////////////////////////////////////

int tauT; //Number of consecutives node aditions   	
int tauD; //Number of consecutives game round robins


void inicia_rand(int semilla);
void nucleo_inicial();
void jugada(); 
void nuevo_nodo(); 


char file[256];  
FILE *fich; 


int main()
{

//////////////////////////
// SET PARAMETERS:

    eps=0.99;
    
    tauT=10;
    tauD=1;

    ro=0.5;
 
    b=1.5;    
    R=1.;
    Pu=0.;
    Su=0.;  
//
//////////////////////////

 //OUTPUT FILE:
    sprintf(file,"%d-Net-%i-%i-%.2lf-%.2lf.dat",N,tauT,tauD,b,eps);        
    fich=fopen(file,"w");

 
//////////////////////////
// INITIALIZE VARIABLES:

    steps=N-m_o;   //Number of nodes that will be added to the network initial core
    s=m_o;         //Initial size of the network
    n_coop= 0;  
    
    for(i=0;i<N;i++)      
    {
	ben[i]=0;
	e[i]=1;          
    }
  
    inicia_rand(time(0));
  
    nucleo_inicial(); //SETS THE INITIAL CORE OF THE NETWORK (m_o nodes connected all-2-all)
//
///////////////////////////

    
///////////////////////////
// NETWORK EVOLUTION LOOP:

    do{
    
        //1.-GAME
	
	for(II=1;II<=tauD;II++) // Loop of Evolutionary Dynamics
	{
	    jugada();  // Game round-robin and strategy exchanges 
	
            //Count Cooperators:

	    n_coop=0;

	    for(i=0;i<s;i++)
	    {
		if(e[i]==0)	    
		    n_coop++;	      
	    }
      
	}

        //2.-Recalculate Attachment Probabilities (that contain the benefit of the nodes in the LAST round-robin)

	for(i=0;i<N;i++)
	{P[i]=0.;}

	P[0]=fit[0];

	for(i=1;i<s;i++)
	    P[i]=P[i-1]+fit[i]; 
	
	norma=P[s-1];
	
        //3.-ATTACHMENT:
	
	
	
	for(II=1;II<=tauT;II++)  //Loop of network growth via EPA
	{
	    if(s<N)     //parche para no pasarme de s=N
	    {		
	    	s++;	    
		nuevo_nodo();  //Incorporation of a new node launching m links 	      							
		printf("s=%d\n",s);
	    }
	    	    
	}
	
}while(s<N);

// END NETWORK EVOLUTION LOOP
///////////////////////////////////////////////
	printf("fin del while\n");
///////////////////////////////////////////////
// WRITE THE RESULTING NETWORK

    for(i=0;i<N;i++)
    {
	fprintf(fich,"%i  %i  %i  ",i,e[i],k[i]);
	
	for(II=0;II<k[i];II++)
	    fprintf(fich,"%i  ",C[i][II]);
	
	fprintf(fich,"\n");
    }
    
    fclose(fich);

//
///////////////////////////////////////////////
  
} //END main()



/////////////////////////////////////////////////
/// INITIALIZE RANDOM NUMBER GENERATOR 

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

// END inicia_rand()
//////////////////////////////////////////////////


//////////////////////////////////////////////////  
// INITIALIZE VARIABLES OF THE CORE OF THE NETWORK

void nucleo_inicial() 
{

// 1.- INITIALIZE TOPOLOGICAL VARIABLES

    for(i=0;i<N;i++)
    {
	k[i] = 0;     
 
	for(j=0;j<K_MAX;j++)
	    C[i][j]=0;
    }
  
  
    for(j=0;j<m_o;j++) 
	k[j]=m_o-1;


    for(i=0;i<m_o;i++)
    {
	counter=0;

	for(j=0;j<m_o;j++)
	{
	    if(i!=j)
	    {
		C[i][counter]=j;
		counter++;
	    }
	}
    }
  
//2.- INITIALIZE FITNESS  

    for(i=0;i<m_o;i++)                 
    { 
	fit[i]=1.-eps;  
    }


//3.- INITIALIZE ATTACHMENT PROBABILITIES

    for(i=0;i<N;i++)
    {
	P[i] = 0;
	P_prov[i]=0;
    }
  
    for(i=0;i<m_o;i++)                    
    {       
	for(j=0; j<=i; j++)
	{                              
	    P[i]=P[i]+fit[j]; 
	}
      
	P_prov[i]=P[i];
    }
  
    norma=P[m_o-1];
    norma2=norma;

//4.- INITIALIZE STRATEGIES (ALL COOPERATORS)
  
    for(i=0;i<m_o;i++)
    {     
	e[i]=0;       
	n_coop++;
    }
  
}

// END nucleo_inicial()
////////////////////////////////////////////////////



////////////////////////////////////////////////////
// Addition of a new node to the network via EPA

void nuevo_nodo() 
{
  
//Initialize variables:

  for(i=0; i<m; i++)
    {unido[i]=0;}
  
  for(i=0; i<N; i++)
    {P_prov[i]=P[i];} 
  
  norma2=norma;
  P[s-1]=0.;  
  P_prov[s-1]=0.;  
  
//Set initial Strategy of the new node:
  
  r=FRANDOM;
  
  if(r<ro)
    {                
      e[s-1]=0; //New node is borned as Cooperator
    }


//Set the m links launched:
  
  for(q=0;q<m;q++)        
  {
      
     //Choose  the node target
      
	v=FRANDOM*norma2;
     
	for(j=0; j<s-1; j++)
	{
	    if(v<P_prov[j])  
	    {
		unido[q]=j;  //j is selected for attachment with th new node
		break;
	    }
	}
      
      //Set the probability of the selected node to 0 (avoid multiple links) 
	
	g=unido[q];
	P_prov[g]=0; 
      
	for(j=g+1; j<s-1; j++)  
	{
	    P_prov[j]=P_prov[j]-fit[g] ; 
	}
	
	norma2=norma2-fit[g];
      
  }    
  
//Update the network topology with the new m links 

  for(j=0;j<m;j++)
  {
      g=unido[j];
      k[g]=k[g]+1;
      C[g][k[g]-1]=s-1;
      C[s-1][j]=g;
  }
  
//Set variables of the new added node
 
  k[s-1]=m;
  fit[s-1]=(1.-eps);
  P[s-1]=P[s-2]+fit[s-1];
  norma=P[s-1];

}

// END nuevo_nodo()
/////////////////////////////////////////////



/////////////////////////////////////////////
// Game Round-robin

void jugada()  
{
    //1.-INITIALIZE:
   
    n_coop=0;

    for(i=0;i<s;i++)      
    {
	ben[i]=0;
    }
    
    //2.-EVERY COUPLE OF LINKED NODES PLAY
    //AND THE NODES ACCUMULATE THE PAYOFF 
    //OBTAINED PLAYING WITH ALL THEIR NEIGHBORS:
   
    for(i=0;i<s;i++)
    {
	for(j=0;j<k[i];j++)
	{
	    y=C[i][j];
            
	    if(y>i)  
	    {  
		if(e[i]==0 && e[y]==0)  //TWO COOPERATORS
		{
		    ben[i]+=R;
		    ben[y]+=R;
		}
		
		if(e[i]==0 && e[y]==1)  //COOPERATOR AND DEFECTOR
		{
		    ben[i]+=Su;
		    ben[y]+=b;
	     	}

		if(e[i]==1 && e[y]==1)  //TWO DEFECTORS
		{
		    ben[i]+=Pu;
		    ben[y]+=Pu;
		}
		
		if(e[i]==1 && e[y]==0)  //DEFECTOR AND COOPERATOR
		{
		    ben[i]+=b;
		    ben[y]+=Su;
		}
	    }
	  
	}
    } //END OF ROUND ROBIN


    for(i=0;i<s;i++)
	e_a[i]=e[i];

    //3.-STRATEGY UPDATES:

    for(i=0;i<s;i++)
    {
        //CHOOSE A NEIGHBOR TO COMPARE
	
	r=FRANDOM;
	r=r*k[i];
	w=(int)r;
    
	t=C[i][w];    
      
        //LOOK FOR THE MAX CONNECTIVITY 
	
	if(k[i]>k[t])
	{k_m=k[i];}
	else
	{k_m=k[t];}


       //IF THE PAYOFF OF t IS LARGER THAN THAT OF i, i MAY CHANGE ITS STRATEGY 
     
	if(ben[i]<ben[t])
	{
	    p=(ben[t]-ben[i]);
	    p=p/(b*(double)k_m);

	    v=FRANDOM;
	    if(v<p)
	    {
		e_a[i]=e[t];
	    
	    }
	}
	else
	{
	    e_a[i]=e[i];
	}

    }
    
    for(i=0;i<s;i++)
	e[i]=e_a[i];

    //4.-COMPUTATION OF NODES FITNESS (THAT WILL ENTER IN THE ATTACHMENT PROBABILITIES)

    for(i=0;i<s;i++)
	fit[i]=1.-eps+eps*ben[i];
       
}

// END jugada()
/////////////////////////////////////////////
