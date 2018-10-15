//////////////////////////////////////////////////
//////////////////////////////////////////////////
///  OJO!! ESTE PROGRAMA SOLO SE USA SI T>=1   ///
///  SI T<1 HAY QUE CAMBIAR LOS SUMANDOS EN LA ///
///  SUBRUTINA "ECUACIONES"                    ///
//////////////////////////////////////////////////
//////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#define Nlatt 1000
#define e 2.718281828

#define NormRANu ((float) (0.5*4.656612595521636e-10))
#define NormRANu2 ((float) (0.5*4.656612595521636e-10))  

unsigned long int irr[256];
unsigned long int ir1;
unsigned char ind_ran,ig1,ig2,ig3;

int I,J,II,III,JJ,marca,NT,IT,t;

double cIo[Nlatt+1],c[Nlatt+1],R,dc[Nlatt+1],dI,Lc,Norma,
       K[Nlatt+1],K1[Nlatt+1],K2[Nlatt+1],gama,
       K3[Nlatt+1],K4[Nlatt+1],cIoo[Nlatt+1],dJJ,dJ,
       T,pi,wi[Nlatt+1],delta_T,Nint,tiempo,P[Nlatt+1],
       Npart,delta_t,dNlatt,Acoplo[Nlatt+1],Acoplo2[Nlatt+1],
       Acoplo3[Nlatt+1],Cmed,C,dNmedia,r,kmed,Time,KMED,fact;


void ini_ran(int SEMILLA)
{
  int INI,FACT,SUM,i;
 
  srand(SEMILLA);

  INI=SEMILLA;
  FACT=67397;
  SUM=7364893;

  for(i=0;i<256;i++)
    {
      INI=(INI*FACT+SUM);
      irr[i]=INI;
    }
  ind_ran=ig1=ig2=ig3=0;

}

/* Esta funcion produce un numero aleatorio entre 0 y 1 */
void Random()
{
  //float r;

  ig1=ind_ran-24;
  ig2=ind_ran-55;
  ig3=ind_ran-61;
  irr[ind_ran]=irr[ig1]+irr[ig2];
  ir1=(irr[ind_ran]^irr[ig3]);
  ind_ran++;
  r=ir1*NormRANu;

  //  return r;
}

void ecuacion();

int main()
{
  FILE *fich1;
  FILE *fich2;
  FILE *fich3;
  char file1[40],file3[40];

  fich2=fopen("data_out/mapCooperacion.dat","w");

  ini_ran(1);
  pi=acos(-1.);



  Nint=10000000;
    
  dNlatt=Nlatt;

  T=1.5;
  delta_T=-0.05;
  NT=3;

  gama=-2.9;
  Norma=0.;

// CONDICIONES INICIALES Y P(k)

  for(J=1;J<=Nlatt;J++)
  {
      cIoo[J]=0.;
      dJ=J;

// TARGETTED:
//      if(J>10)               //distribucion inicial de cooperadores sobre los de k > k*
//      {cIoo[J]=1.;}

// RANDOM:
      Random();        //distribucion inicial de cooperadores al azar
      cIoo[J]=r;

//SCALE-FREE:
      P[J]=pow(dJ,gama);      //creo la distribucion de conectividad

//ERDOS-RENYI
/*      fact=1;
      for(I=1;I<=J;I++)
      {
	  dJ=I;
	  fact=fact*dJ;      
      }
      dJ=J;
      KMED=(gama-1.)/(gama-2.);
      P[J]=pow(e,-KMED)*pow(KMED,dJ)/fact;
*/
      Norma=Norma+P[J];
  }
  
  kmed=0.;

  for(J=1;J<=Nlatt;J++)
  {
      dJ=J;
      P[J]=P[J]/Norma;    //normalizo la p(k) y calculo <k> de una SF
      kmed=kmed+dJ*P[J];
  }


  for(IT=1;IT<=NT;IT++)    //bucle para barrer sobre distintos valores de T
  {

      printf("%lf\n",T);
     
      marca=1;


      sprintf(file3,"data_out/mapEvolution%.4lf.dat",T);
      fich3=fopen(file3,"wt");
      fclose(fich3); 

      sprintf(file1,"data_out/mapAll-Evolution%.4lf.dat",T);
      fich1=fopen(file1,"wt");
      fclose(fich1);    

      tiempo=0.;
      
      Cmed=0.;
      
      for(J=1;J<=Nlatt;J++)
      {
	  cIo[J]=cIoo[J];
      }

      for(II=1;II<=Nint;II++)   // bucle iteraciones del mapa
      {  

	
	  for(J=1;J<=Nlatt;J++)
	  {
	      c[J]=cIo[J];
	  }
	  
	  Lc=0.;

	  for(J=1;J<=Nlatt;J++)
	  {
	      dJ=J;         //dJ es double,  J es entero
	      Lc=Lc+(dJ*P[J]*c[J]);
	  }

	  Lc=Lc/kmed;       //prob de que un vecino sea cooperador


	  for(I=1;I<=Nlatt;I++)
	  {
	      ecuacion();   //   calculo c_k(t+1)    ---->  ==actua sobre  cIo[I]
	  }
	
	 
	  C=0.;
	  
	  for(J=1;J<=Nlatt;J++)
	  {
	      C=C+cIo[J]*P[J];     //calculo el tot de coop (ya indep de k)
	  }
	  
	  fich3=fopen(file3,"at");        
	  fprintf(fich3,"%i  %lf\n",II,C);
	  fclose(fich3);

	  if(II>(Nint-2000))
	  {
	      Cmed=Cmed+C;
	  }

	  if(II==marca)
	  {
	      
	      fich1=fopen(file1,"at");        
	      for(J=1;J<=Nlatt;J++)
	      {
		  fprintf(fich1,"%i  %i  %lf\n",II,J,cIo[J]);
	      }
	      fprintf(fich1,"\n");
 
	      fclose(fich1);

	      marca=marca+50;
	      printf("Coop=%lf\n",C);
	  }
	
      }
    
      Cmed=Cmed/2000.;    // hago la media de cooperacion cada 2000 iteraciones del mapa
    
      // fclose(fich1);

      fprintf(fich2,"%lf  %lf\n",T,Cmed);      

      printf("%lf  %lf\n\n",T,Cmed);

      T=T+delta_T;
    
  }

  
  fclose(fich2);
  
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


void ecuacion()          
  {

      int tope;
      double dtope;
      
//PRIMERO:

      dI=I;
      dtope=dI*T;
      tope=dtope;
      tope=tope+1;

      Acoplo[I]=0.;

      for(JJ=1;JJ<=Nlatt;JJ++)
      {
	  if(JJ>=dtope)
	  {	  
	      dJJ=JJ;
	      Acoplo[I]=Acoplo[I]+P[JJ]*Lc*((dJJ-T*dI)*(1.-c[I])*c[JJ] 
					-(T*dJJ-dI)*(1.-c[JJ])*c[I]);
	  }
      }

      Acoplo[I]=Acoplo[I]/(T*kmed);

//SEGUNDO:

      dtope=dI/T;
      tope=dtope;
      tope=tope+1;

      Acoplo2[I]=0.;

      for(JJ=1;JJ<=I;JJ++)
      {
	  if(JJ>=dtope)
	  {dJJ=JJ;
	  Acoplo2[I]=Acoplo2[I]+dJJ*P[JJ]*Lc*(T*dJJ-dI)*(1.-c[JJ])*c[I];
	  }
      }

      Acoplo2[I]=Acoplo2[I]/(T*dI*kmed);

//TERCERO:

      dtope=dI*T;
      tope=dtope;

      if(tope>Nlatt)
      {tope=Nlatt;}

      Acoplo3[I]=0.;
      for(JJ=I+1;JJ<=Nlatt;JJ++)
      {
	  if(JJ<=dtope)
	  {
	      dJJ=JJ;
	      Acoplo3[I]=Acoplo3[I]+P[JJ]*Lc*(T*dJJ-dI)*(1.-c[JJ])*c[I];
	  }
	  else
	  {
	      break;
	  }
      }

      Acoplo3[I]=Acoplo3[I]/(T*kmed);

      cIo[I]=c[I]+Acoplo[I]-Acoplo2[I]-Acoplo3[I];
  }


  
 
