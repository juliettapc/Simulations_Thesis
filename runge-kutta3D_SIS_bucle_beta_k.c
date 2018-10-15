//programa para implementar el metodo de Runge-Kutta de resolucion
			//de ecuaciones diferenciales (en tres dimensiones)
//16-3-2010: aplicado a la resolucion del modelo SIR de epidemias



# include <time.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>





# define N 2    //numero de dimensiones (num. ec. dif. del sist.)



# define beta_min   0.0  //tasa de infeccion (como es una prob. va de 0 a 1)
# define beta_max   1.0  
# define delta_beta  0.01    


# define gamma  1.0     //tasa de recuperacion

# define k    6     //numero de contactos


double x_old[N+1],  k_1[N+1], k_2[N+1], k_3[N+1], k_4[N+1];
double   t, t_final,t_o, delta_t, beta;
int i,j;



double funcion1(double, double );
double funcion2(double, double );


char file1[256], file2[256];

FILE *escribe1;
FILE *escribe2;

int main()
{
  

  sprintf(file2,"SIS_vs_beta%.2lf-%.2lf_gamma%.2lf_k%d.dat",beta_min,beta_max,gamma,k);   
  escribe2=fopen(file2,"wt");
  fclose(escribe2);


  beta=beta_min;
  
  while(beta<=beta_max)
    {
      
      
      for(i=1;i<=N;i++)    //k_*[i=1]  corresponde con las m1,m2,m3,m4
	{                  //k_*[i=2]  corresponde con las k1,k2,k3,k4
	  k_1[i]=0;        
	  k_2[i]=0;
	  k_3[i]=0;
	  k_4[i]=0;
	  
	  x_old[i]=0.0;
	}
      
      
      
      
      x_old[2]=0.001;  //c.i.  infectados
      x_old[1]=1.0-x_old[2];  //c.i. susceptibles
      
      
      
      sprintf(file1,"Evol_temp_SIS_beta%.2lf_gamma%.2lf_i_ini%lf_k%d.dat",beta,gamma,x_old[2],k);   
      escribe1=fopen(file1,"wt");
      fclose(escribe1);
      
      
      
      
      
      t=0;
      
      t_final=50;
      
      delta_t=0.01;
      
      
      
      
      
      escribe1=fopen(file1,"wt");    
      fprintf(escribe1,"%lf  %lf  %lf  %lf\n",t,x_old[1],x_old[2],x_old[1]+x_old[2]); 	     
      printf("%f+%f=%f\n",x_old[1],x_old[2],x_old[1]+x_old[2]);
      //getchar();
      
      
      while(t<=t_final)
	{                 
	  
	  
	  k_1[1]=funcion1(x_old[1],x_old[2])*delta_t;    //las m1, m2, ...
	  
	  k_2[1]=funcion1(x_old[1]+0.5*k_1[1],x_old[2]+0.5*k_1[2])*delta_t;	  
	  
	  k_3[1]=funcion1(x_old[1]+0.5*k_2[1],x_old[2]+0.5*k_2[2])*delta_t;	  
	  
	  k_4[1]=funcion1(x_old[1]+k_3[1],x_old[2]+k_3[2])*delta_t;
	  
	  
	  
	  
	  
	  k_1[2]=funcion2(x_old[1],x_old[2])*delta_t;     //las k1, k2, ...
	  
	  k_2[2]=funcion2(x_old[1]+0.5*k_1[1],x_old[2]+0.5*k_1[2])*delta_t;	  
	  
	  k_3[2]=funcion2(x_old[1]+0.5*k_2[1],x_old[2]+0.5*k_2[2])*delta_t;	  
	  
	  k_4[2]=funcion2(x_old[1]+k_3[1],x_old[2]+k_3[2])*delta_t;
	  
	  
	  
	  
	  
	  
	  
	  
	  x_old[1]=x_old[1] + (1.0/6.0)*(k_1[1]+2*k_2[1]+2*k_3[1]+k_4[1]);	 
	  x_old[2]=x_old[2] + (1.0/6.0)*(k_1[2]+2*k_2[2]+2*k_3[2]+k_4[2]);
	  //ojo con la fraccion, pq 1/6=0!!!!!
	  
	  
	  
	  fprintf(escribe1,"%lf  %lf  %lf  %lf  %lf  %d\n",t,x_old[1],x_old[2],x_old[1]+x_old[2],gamma,k); 
	  
	  
	  
	  if ( x_old[2]<=0.0 ) //si no quedan infectados  
	    {                         
	      break;                   
	    }      
	  
	  
	  
	  //printf("%f+%f=%f\n",x_old[1],x_old[2],x_old[1]+x_old[2]);
	  //getchar();
	  
	  t+=delta_t;                
	}            //fin del algoritmo      
      fclose(escribe1);


      printf("%f: s%f  i%f\n",beta,x_old[1],x_old[2]);

      escribe2=fopen(file2,"at");    
      fprintf(escribe2,"%lf  %lf  %lf  %lf  %d\n",beta,x_old[1],x_old[2],gamma,k);
      fclose(escribe2);
      
      beta+=delta_beta;
    } 
  exit(0);
  
}



//////////////////////////////
//////////////////////////////////////
/////////////////////////////////



double funcion1(double x,double y)  //ec. dif. del numero de susceptibles      
{
  double resultado;
  
  
  resultado=-beta*y*x*k+gamma*y;
  return (resultado); 
    
}




//////////////////////////////
//////////////////////////////////////
/////////////////////////////////



double funcion2(double x,double y)  //ec. dif. del numero de infectados      
{
  double resultado;
  
  
  resultado=beta*y*x*k-gamma*y;
  return (resultado); 
    
}



