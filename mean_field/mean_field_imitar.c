//     Mean-Field:   Implentación de las ecuaciones dif. de la dinámica replicador
//     sobre redes complejas, en aproximacion de campo medio
//  Compartimentalizacion de los cooperadores y defectores en clases de conectividad
// Regla de imitación incondicional, no probabilistica (elimina el factor topológico en la dinámica)


# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# define  K_max  1000  // ---->>> ponerlo a mano.
# define  k_star  50     //para el targetted

# define t_fin 100000
# define gamma 3.0      //  exponente de la power law  

# define rho   0.9   // concentracion inicial independiente de k (para cuando no targetted)

# define T   1.2
# define R   1.0
# define P   0.0
# define S   0.0

//# define  ER

# define SF


# define targetted


double Pk[K_max+1];
long double k_factorial;   //para la distribucion de conectividad  (SF o exp)
int k,k_o, t,entero,k1_inf,k2_inf;
double P_dc[K_max+1],P_cd[K_max+1]; // Probabilidades de cambio de estrat para cada clase de conectividad
double PayC[K_max+1],PayD[K_max+1],dif_pay_cd,dif_pay_dc;  //Payoff de un coop. y defect. con conectiv. k
double lc;   // Prob de que un vecino sea coop.  y exponente de la distrib.
double Ck[K_max+1],Ck_new[K_max+1], Ck_dt[K_max+1],C_med;   //fraccion de coop con conectiv.  k  y su derivada
double aux, aux1, tita_cd,tita_dc,beta,norma;
 int i, jj,kk;   ///contadores
double k_media,max;     /////////// //luego tomaré la parte entera????????

char nombre1[256], nombre2[256];
FILE *fich1,*fich2; 

int main()
{
    
  
    
    k_o=1;   // conectividad mínima
    k_media=0.0;      ///ojo! no va metida a mano!

    

#ifdef ER
    sprintf(nombre1,"Mean-field_ER_Kmax%d_T%.2lf_t%d.dat",K_max,T,t_fin);    // fichero para guardar la <c>(t)
    fich1=fopen(nombre1,"wt");
    fclose(fich1);
    
    sprintf(nombre2,"Pk_Mean-field_ER_Kmax%d_T%.2lf_t%d.dat",K_max,T,t_fin);    // fichero para guardar la P(k)
    fich2=fopen(nombre2,"wt");
    fclose(fich2);
#endif
    
    
 
#ifdef SF
    sprintf(nombre1,"Mean-field_SF_Kmax%d_T%.2lf_t%d.dat",K_max,T,t_fin);    // fichero para guardar la <c>(t)
    fich1=fopen(nombre1,"wt");
    fclose(fich1);
    
    sprintf(nombre2,"Pk_Mean-field_SF_Kmax%d_T%.2lf_t%d.dat",K_max,T,t_fin);    // fichero para guardar la P(k)
    fich2=fopen(nombre2,"wt");
    fclose(fich2);
#endif   








# ifdef targetted      //sólo para el caso de SF
 
#ifdef SF
    sprintf(nombre1,"Mean-field_SF_Kmax%d_targetted_k%d_T%.2lf_t%d.dat", K_max,k_star,T,t_fin);    // fichero para guardar la <c>(t)
    fich1=fopen(nombre1,"wt");
    fclose(fich1);
    
    sprintf(nombre2,"Pk_Mean-field_SF_Kmax%d_targetted_k%d_T%.2lf_t%d.dat", K_max,k_star,T,t_fin);    // fichero para guardar la P(k)
    fich2=fopen(nombre2,"wt");
    fclose(fich2);
#endif   

# endif





    
    
    for(i=k_o;i<=K_max;i++)      //establezco concentraciones iniciales
    {
	Ck[i]=rho;
       
	//printf("Ck[%d]=%lf\n",i,Ck[i]);	
	
	Ck_dt[i]=0.0;     //////////////////  Velocidad inicial =0    ???
	Ck_new[i]=0.0;
    }
 
  
    
# ifdef targetted
    
    for (i=k_o;i<=K_max;i++)
    {
	Ck[i]=0.0;
	if(i>=k_star)
	    Ck[i]=1.0;	
    }
    
# endif
    

     
    
    norma=0.0;
    for(i=k_o;i<=K_max;i++)      //defino la distribucion de conectividad
    {

#ifdef SF	 
	aux=pow(k_o,(gamma-1.0));     
	aux1=pow(i,(-gamma));
	Pk[i]=(gamma-1.0)*aux1/aux;
	norma+=Pk[i];       
#endif



#ifdef ER	 
	k_factorial=1;
	for(jj=1;jj<=i;jj++)
	{
	    aux=jj;
	    k_factorial=k_factorial*aux;
	
	    //printf("factorial de %d=%lf\n",i,k_factorial);	
	}

	k_media=(gamma-1.)*k_o/(gamma-2.);  ///para la ER sí la pongo yo a mano (expresion copiada de Jesus)

	Pk[i]=exp(-k_media);
//printf("P[%d]=%lf\n",i,Pk[i]);

	aux=pow(k_media,i);             
	Pk[i]=Pk[i]*aux/k_factorial;

	//printf("P[%d]=%lf\n",i,Pk[i]);
	norma+=Pk[i];  
# endif


	
    }    
    



    
    fich2=fopen(nombre2,"at");                   
    
    for(i=k_o;i<=K_max;i++)      //normalizo  la distribucion de conectividad
    {
	
	Pk[i]=Pk[i]/norma;
	
	fprintf(fich2,"%d   %lf\n",i,Pk[i]);     //  imprimo la  P(k)
	
    }
    fclose (fich2);	
    
    //printf("Hecho el histograma P(k)\n");



#ifdef SF	 
//  Calculo la k-media   (solo para la SF!!!)
    
       k_media=0.0;
    for(i=k_o;i<=K_max;i++)      //normalizo  la distribucion de conectividad
    {
	k_media+=i*Pk[i];
    }
#endif


    printf("k_media=%lf\n",k_media);

//getchar();



 t=1;
    
 do      //////////////////////////////////// //  bucle pasos de tiempo
    {   	
	//printf("t=%d\n",t);


//   Calculo la <C>(t)  y la imprimo

	C_med=0.0;
	for(i=k_o;i<=K_max;i++)       
	{      
	    //printf("Ck[%d]=%lf\n",i,Ck[i]);
	    C_med+=Ck[i]*Pk[i];               ///////////  falta normalizar  ???????!!!!!
	}


	fich1=fopen(nombre1,"at");
	fprintf(fich1,"%d   %lf\n",t,C_med);
	fclose (fich1);	
	
	printf("<C>(%d)=%lf\n",t,C_med);

	//getchar();

// Calculo la prob de tener un vecino coop  para cada clase de conectividad

	for(i=k_o;i<=K_max;i++)
	{	
	    
	    PayC[i]=0.0;	
	    PayD[i]=0.0;
	}

	for(i=k_o;i<=K_max;i++)         //Calculo los payoffs de cada clase de conectividad
	{
	    lc=0;
	    for (kk=k_o;kk<=K_max;kk++)    //sumatorio en k'
	    {
		lc+=kk*Pk[kk]*Ck[kk]/k_media;
	    }
	    

	    PayC[i]=i*lc;	
	    PayD[i]=T*i*lc;	    //        ojo! válido sólo para P=S=0; R=1   (PD)
	    //   printf("PayC[%d](%d)=%lf\n",i,t,PayC[i]); 
	    //   printf("PayD[%d](%d)=%lf\n",i,t,PayC[i]); 
	    
	}
	//	getchar();	
	
	
// Calculo las probabilidades de cambio de estrat para cada clase de conectividad:  imitacion incondicional
	
	
	for(i=k_o;i<=K_max;i++)
	{

	    P_dc[i]=0.0;
	    P_cd[i]=0.0;
	    
     k1_inf=T*i;
     k2_inf=i/T;

	    for (kk=k1_inf;kk<=K_max;kk++)    //sumatorio en k'
	    {
		P_dc[i]=P_dc[i]+ kk*Pk[kk]*Ck[kk]/k_media;

        //printf("P_dc[%d](%d)=%lf\n",kk,t,P_dc[i]);
	    }	    

         for (kk=k2_inf;kk<=K_max;kk++)    //sumatorio en k'
	    {
          P_cd[i]=P_cd[i]+ kk*Pk[kk]*(1.0-Ck[kk])/k_media;
          //printf("P_cd[%d](%d)=%lf\n",kk,t,P_cd[i]);
         }

	     // getchar();
	}
	
	
	
//  Las ecuaciones diferenciales  -->>  mapa discreto:   Ck(t+1)=Ck(t)+Ck_dt(t)   (delta_t =1)
	
	for(i=k_o;i<=K_max;i++)
	{
	    Ck_dt[i]=(1.0-Ck[i])*P_dc[i]-Ck[i]*P_cd[i];   //calculo las velocidades
	    // printf("C_dt[%d](%d)=%lf\n",i,t,Ck_dt[i]); 
	    
	}
	//getchar();	


	for(i=k_o;i<=K_max;i++)
	{
	    Ck_new[i]=Ck[i]+Ck_dt[i];
	    // printf("Ck[%d](%d)=%lf\n",i,t,Ck[i]); 
	}
	
	//getchar();



	for(i=k_o;i<=K_max;i++)       //    reasigno
	{
	    Ck[i]=Ck_new[i];
	}


	t++;
    }
    while (t<=t_fin);       
    
       
}
