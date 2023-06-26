#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
//#include <omp.h>

// Metropolis Multihit Algorithmus

const int L = 128;
void rb();
void sweep(double b);
int Ham(int i, int j);
int lat[L+2][L+2] = {0};

double max = (double)RAND_MAX;
const int mh = 5;
double rand_num;
double Hami1;


int main()
{
    // Zeitmessung
    time_t t1;
    time_t t2;
    time_t t3;
    t1 = time(NULL);
    srandom(time(0));

    // hot start
    for (int i=1;i<L+1;++i)
    {
        for (int j=1;j<L+1;++j)
        {
            lat[i][j]= -1+2*(random()%2);
            

        }
    }
    
    // Temperatur
    double beta[35] ={0};
    for (int i =0; i<30; ++i)
    beta[i+2]= 0.15 +0.6/30*(double)i;
    for (int i =0;i<2;++i)
    beta[i] = 0.1+0.04*(double)i;

    for (int i=0;i<3;++i)
    beta[i+32]= 0.8+ 0.1*(double)i;

    // Datei aller Messdaten
    FILE *datei;
    datei = fopen("result3a.txt","w");

    // Datei mit Mittelwerten
    FILE *datei2;
    datei2 = fopen("short_results3a.txt","w");

    //#pragma omp parallel for
    for (int n=0; n<35;++n)
    {
        // thermalisation
        if (beta[n]<=0.4)
        {
        //termalization
        for (int t=0;t<6000;++t)
        {
            sweep(beta[n]);
        }
        }
        else if (beta[n]>0.4 & beta[n]<0.47)
        {
        //termalization
        for (int t=0;t<10000;++t)
        {
            sweep(beta[n]);
        }
        }
        else if (beta[n]>=0.47 & beta[n]<=0.7)
        {
        //termalization
        for (int t=0;t<7000;++t)
        {
            sweep(beta[n]);
        }
        }
        else
        {
        //termalization
        for (int t=0;t<20000;++t)
        {
            sweep(beta[n]);
        }
        }
        
        
    
        // Arrays in die die Messungen geschrieben werden Anzahl Messungen= iter
        int Hamilton_store[300] ={0};
        int mag_store[300] = {0};
        int mag_sq_store[300] = {0};
        int iter = 300;
        
        for ( int N =0; N<iter; ++N)
        {
            
            for ( int i=1; i<L+1; ++i)
            {
                for (int j =1; j<L+1; ++j)
                {
                    // Gesamthamiltonian
                    Hamilton_store[N] += Ham(i,j);
                    // Gesamtmagnetisierung
                    mag_store[N] += lat[i][j];
                }
            }
            
            // Betrag der Magnetisierung
            if (mag_store[N]<0)
            mag_store[N] *= -1;
            // mag_sq_store[N] = mag_store[N]*mag_store[N];

        // Dekorrelierung
        if (beta[n]<=0.4)
        {
        //termalization
        for (int t=0;t<6000;++t)
        {
            sweep(beta[n]);
        }
        }
        else if (beta[n]>0.4 & beta[n]<0.47)
        {
        //termalization
        for (int t=0;t<10000;++t)
        {
            sweep(beta[n]);
        }
        }
        else if (beta[n]>=0.47 & beta[n]<=0.7)
        {
        //termalization
        for (int t=0;t<7000;++t)
        {
            sweep(beta[n]);
        }
        }
        else
        {
        //termalization
        for (int t=0;t<20000;++t)
        {
            sweep(beta[n]);
        }
        }
              
    
            fprintf(datei,"%i\t",Hamilton_store[N]);
            fprintf(datei,"%i\n",mag_store[N]);


            

        }
       

    
        
        
        double energy_mean= 0;
        double mag_mean =0 ;
        // double mag_sq_mean = 0;
        for (int N =0; N<iter;++N)
        {

            energy_mean += (double)Hamilton_store[N]/iter;
            mag_mean += (double)mag_store[N]/iter;
            // mag_sq_mean += (double)mag_sq_store[N]/(double)iter;

        }
        
       
        fprintf(datei2,"%.6e\t",energy_mean);
        fprintf(datei2,"%.6e\n",mag_mean);

       

        

        // Zwischenergebnis
        if ((n+1)%5==0)
        {
            printf("%i\n",n);
            t2 = time(NULL);
            printf("%ld\n",t2-t1);
        }

        


    }
    

    fclose(datei);   
    fclose(datei2);
    t2=time(NULL);
    printf("%ld",t2-t1);


    
   

    

    return 0;
}


void rb()
{
    //#pragma omp parallel for
    for (int i=0;i<L+2;++i)
    {
        lat[i][L+1] = lat[i][1];
        lat[L+1][i] = lat[1][i];
        lat[i][0] = lat[i][L];
        lat[0][i] = lat[L][i];
    }
    

    
}

int Ham(int i, int j)
{
    int H;
    H = lat[i][j]*(lat[i-1][j]+lat[i][j+1]+lat[i][j-1]+lat[i+1][j]);
    return -H;


}

void sweep(double b)
{
    //#pragma omp parallel for
    for (int i=1;i<L+1;++i)
    {
        for (int j=1;j<L+1;++j)
        {
            Hami1 = (double) Ham(i,j);
            
            for (int mhits=0;mhits<mh;++mhits)
            {
                // ich halte mich mit dieser Zeile an den Pseudo Code aus dem Skript, bei dem auch mit einer 50%igen Wahrscheinlichkeit nichts passiert.
                if (random()%2==0)
                {
                    // wenn der (Teil-)Hamiltonian positiv ist ist auch die differenz von H_i und H_j negativ
                    if (Hami1>=0)
                    {
                        lat[i][j] *= -1;
                        Hami1 *= -1;
                        
                    }
                    else
                    {
                        // die differenz ist 2 * Hami
                        rand_num = ((double)random())/(max);
                        if (exp(b*2*Hami1)>=rand_num)
                        {
                            lat[i][j] *= -1;
                            Hami1 *= -1;
                        }
                    }
                }
            }
        }
    }           
    rb();
}

