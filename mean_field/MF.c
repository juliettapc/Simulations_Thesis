#include <stdio.h>
#include <math.h>
#define Nlatt 1000


#define NormRANu ((float) (0.5*4.656612595521636e-10))
#define NormRANu2 ((float) (0.5*4.656612595521636e-10))  

unsigned long int irr[256];
unsigned long int ir1;
unsigned char ind_ran,ig1,ig2,ig3;


int I,J,II,III,JJ,marca,NT,IT,t;

double cIo[Nlatt+1],c[Nlatt+1],R,dc[Nlatt+1],dI,Lc,Norma,
    K[Nlatt+1],K1[Nlatt+1],K2[Nlatt+1],gama,                  //se utiliza el Runge-kutta para la resuluc numerica de la ec. difer.
       K3[Nlatt+1],K4[Nlatt+1],cIoo[Nlatt+1],dJJ,dJ,
       T,pi,wi[Nlatt+1],delta_T,Nint,tiempo,P[Nlatt+1],
       Npart,delta_t,dNlatt,Acoplo[Nlatt+1],Acoplo2[Nlatt+1],
       Acoplo3[Nlatt+1],Cmed,C,dNmedia,r,kmed;



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

  fich2=fopen("data_out/Cooperacion.dat","w");

  ini_ran(1);
  pi=acos(-1.);


  delta_t=0.05; 
  T=2000.;
  Npart=T/delta_t;
  Nint=Npart;
    
  dNlatt=Nlatt;

  T=0.95;
  delta_T=0.02;
  NT=10;

  gama=-3.;
  Norma=0.;

  for(J=1;J<=Nlatt;J++)         //establezco condiciones iniciales de cooperacion (targetted o al azar), así como la P(k) SF o ER
  {
      cIoo[J]=0.;
      dJ=J;

// TARGETTED:
//      if(J>10)
//      {cIoo[J]=1.;}

// RANDOM:
      Random();
      cIoo[J]=r;

//SCALE-FREE:
      P[J]=pow(dJ,gama);

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

  for(J=1;J<=Nlatt;J++)            //normalizo p(k) y calculo <k>
  {
      dJ=J;
      P[J]=P[J]/Norma;
      kmed=kmed+dJ*P[J];
  }


  for(IT=1;IT<=NT;IT++)    //bucle para barrer en T
  {

      printf("%lf\n",T);
     
      marca=1;


      sprintf(file3,"data_out/Evolution%.4lf.dat",T);
      fich3=fopen(file3,"wt");
      fclose(fich3); 

      sprintf(file1,"data_out/All-Evolution%.4lf.dat",T);
      fich1=fopen(file1,"wt");
      fclose(fich1);    

      tiempo=0.;
      for(J=1;J<=Nlatt;J++)
      {
	  cIo[J]=0.;
	  c[J]=0.;
	  dc[J]=0.;
      }
      
      Cmed=0.;
      
      for(J=1;J<=Nlatt;J++)
      {
	  cIo[J]=cIoo[J];
      }

      for(II=1;II<=Nint;II++)        //pasos del runge-kutta
      {  

	
	  for(J=1;J<=Nlatt;J++)      // condiciones iniciales para el algoritmo
	  {
	      c[J]=cIo[J];
	  }
	  
	  Lc=0.;

	  for(J=1;J<=Nlatt;J++)    //calculo lc (probab. de que un vecino sea cooperador)
	  {
	      dJ=J;
	      Lc=Lc+(dJ*P[J]*c[J]);
	  }

	  Lc=Lc/kmed;


	  for(I=1;I<=Nlatt;I++)
	  {
	      ecuacion();      //dentro de esta funcion se calcula, ademas de dc[], K[]=delta_t*dc[]
	      K1[I]=K[I];
	      c[I]=cIo[I]+K1[I]/2.;	  
	  }
	
	  Lc=0.;

	  for(J=1;J<=Nlatt;J++)
	  {
	      dJ=J;
	      Lc=Lc+(dJ*P[J]*c[J]);
	  }

	  Lc=Lc/kmed;


	  tiempo=tiempo+delta_t/2.;
	  
	  for(I=1;I<=Nlatt;I++)
	  {
	      ecuacion();
	      K2[I]=K[I];
	      c[I]=cIo[I]+K2[I]/2.;
	  }	
       
	  Lc=0.;

	  for(J=1;J<=Nlatt;J++)
	  {
	      dJ=J;
	      Lc=Lc+(dJ*P[J]*c[J]);
	  }

	  Lc=Lc/kmed;

	  for(I=1;I<=Nlatt;I++)
	  {
	      ecuacion();
	      K3[I]=K[I];
	      c[I]=cIo[I]+K3[I];
	  }
	
	  Lc=0.;

	  for(J=1;J<=Nlatt;J++)
	  {
	      dJ=J;
	      Lc=Lc+(dJ*P[J]*cIo[J]);
	  }

	  Lc=Lc/kmed;

	  for(I=1;I<=Nlatt;I++)
	  {
	      ecuacion();
	      K4[I]=K[I];
	      cIo[I]=cIo[I]+K1[I]/6. +K2[I]/3. + K3[I]/3. + K4[I]/6.;
	  }
	
	  C=0.;
	  
	  for(J=1;J<=Nlatt;J++)
	  {
	      C=C+cIo[J]*P[J];
	  }
	  
	  fich3=fopen(file3,"at");        
	  fprintf(fich3,"%lf  %lf\n",tiempo,C);
	  fclose(fich3);

	  if(II>(Nint-20000))
	  {
	      Cmed=Cmed+C;
	  }

	  if(II==marca)
	  {
	      
	      fich1=fopen(file1,"at");        
	      for(J=1;J<=Nlatt;J++)
	      {
		  fprintf(fich1,"%lf  %i  %lf\n",tiempo,J,cIo[J]);
	      }
	      fprintf(fich1,"\n");
 
	      fclose(fich1);

	      marca=marca+50;
	      printf("Coop=%lf\n",C);
	  }
	
      }
    
      Cmed=Cmed/20000.;
    
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

//TOTAL:
      
      dc[I]=Acoplo[I]-Acoplo2[I]-Acoplo3[I];
     
/*
      if(c[I]>=1.000)
      {
	  if(dc[I]>0.000)
	  {
	      dc[I]=0.;
	      c[I]=1.000;
	  }
      }    

      if(c[I]<=0.000)
      {
	  if(dc[I]>0.000)
	  {
	      dc[I]=0.;
	      c[I]=0.000;
	  }
      }      
*/

      K[I]=delta_t*dc[I];
  }


  
 
