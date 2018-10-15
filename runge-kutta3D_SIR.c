//programa para implementar el metodo de Runge-Kutta de resolucion
			//de ecuaciones diferenciales (en tres dimensiones)
//16-3-2010: aplicado a la resolucion del modelo SIR de epidemias



# include <time.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>





# define N 3    //numero de dimensiones (num. ec. dif. del sist.)



# define lambda   0.20  //tasa de infeccion
# define mu   1.0      //tasa de recuperacion
# define k    6       // numero medio de contactos por individuo y unid. tpo



double x_old[N+1],  k_1[N+1], k_2[N+1], k_3[N+1], k_4[N+1];
double   t, t_final,t_o, delta_t;
int i,j;



double funcion1(double, double );
double funcion2(double, double );
double funcion3(double);

char file1[256];

FILE *escribe1;


int main()
{
  


  
  for(i=1;i<=N;i++)    //k_*[i=1]  corresponde con las m1,m2,m3,m4
    {                  //k_*[i=2]  corresponde con las k1,k2,k3,k4
      k_1[i]=0;        //k_*[i=3]  corresponde con las l1,l2,l3,l4
      k_2[i]=0;
      k_3[i]=0;
      k_4[i]=0;

      x_old[i]=0.0;
    }
  

  
  x_old[3]=0.0;  //c.i. reestablecidos
  x_old[2]=0.0001;  //c.i.  infectados
  x_old[1]=1.0-x_old[2]-x_old[3];  //c.i. susceptibles
  


  sprintf(file1,"SIR_lambda%.2lf_mu%.2lf_k%d_rho_ini%lf.dat",lambda,mu,k,x_old[2]);   
  escribe1=fopen(file1,"wt");
  fclose(escribe1);
  


  

  t=0;

  t_final=1000;
  
  delta_t=0.001;
  
 
  
  

  escribe1=fopen(file1,"wt");
  while(t<=t_final)
    {
            
      fprintf(escribe1,"%f  ",t);
      
      
      
      k_1[1]=funcion1(x_old[1],x_old[2])*delta_t;    //las m1, m2, ...
      
      k_2[1]=funcion1(x_old[1]+0.5*k_1[1],x_old[2]+0.5*k_1[2])*delta_t;	  
      
      k_3[1]=funcion1(x_old[1]+0.5*k_2[1],x_old[2]+0.5*k_2[2])*delta_t;	  
      
      k_4[1]=funcion1(x_old[1]+k_3[1],x_old[2]+k_3[2])*delta_t;
      
      
	  

      
      k_1[2]=funcion2(x_old[1],x_old[2])*delta_t;     //las k1, k2, ...
      
      k_2[2]=funcion2(x_old[1]+0.5*k_1[1],x_old[2]+0.5*k_1[2])*delta_t;	  
      
      k_3[2]=funcion2(x_old[1]+0.5*k_2[1],x_old[2]+0.5*k_2[2])*delta_t;	  
      
      k_4[2]=funcion2(x_old[1]+k_3[1],x_old[2]+k_3[2])*delta_t;
      
      
	  
      
      
      k_1[3]=funcion3(x_old[2])*delta_t;    //las l1, l2, ...
      
      k_2[3]=funcion3(x_old[2]+0.5*k_1[2])*delta_t;	  
      
      k_3[3]=funcion3(x_old[2]+0.5*k_2[2])*delta_t;	 	  
      
      k_4[3]=funcion3(x_old[2]+k_3[2])*delta_t;
      
      
      
      
      
      for(i=1;i<=N;i++)   
	{	  	 
	  x_old[i]=x_old[i] + (1.0/6.0)*(k_1[i]+2*k_2[i]+2*k_3[i]+k_4[i]);
	  //ojo con la fraccion, pq 1/6=0!!!!!
	  
	  fprintf(escribe1,"%lf  ",x_old[i]); 
	}

         fprintf(escribe1,"%lf  ",x_old[1]+x_old[2]+x_old[3]);   
	 fprintf(escribe1,"\n"); 
      
      if ( x_old[1]<=0.0 || x_old[2]<=0.0 || x_old[3]>=1.0 )   //si todos estan ya recuperados
	{                                     // o no quedan susceptibles o infectados
	  break;
	}      
      

      printf("%f+%f+%f=%f\n",x_old[1],x_old[2],x_old[3],x_old[1]+x_old[2]+x_old[3]);
      //getchar();

      t+=delta_t;                
    }            //fin del algoritmo
  
  fclose(escribe1);
  
  
  
  exit(0);
  
}



//////////////////////////////
//////////////////////////////////////
/////////////////////////////////



double funcion1(double x,double y)  //ec. dif. del numero de susceptibles      
{
  double resultado;
  
  
  resultado=-lambda*k*y*x;
  return (resultado); 
    
}




//////////////////////////////
//////////////////////////////////////
/////////////////////////////////



double funcion2(double x,double y)  //ec. dif. del numero de infectados      
{
  double resultado;
  
  
  resultado=-mu*y+lambda*k*y*x;
  return (resultado); 
    
}




//////////////////////////////
//////////////////////////////////////
/////////////////////////////////



double funcion3(double y)  //ec. dif. del numero de recuperados      
{
  double resultado;
  
  
  resultado=mu*y;
  return (resultado); 
    
}
