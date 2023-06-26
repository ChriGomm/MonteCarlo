#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
//#include <omp.h>

const int L = 128;
void rb();

void sweep(double b,double h);
int Ham(int i, int j);
int lat[L+2][L+2] = {0};

double max = (double)RAND_MAX;
double rand_num;
double prob;
double Hami1;



int main()
{

    clock_t t1;
   clock_t t2;
    time_t t3;
    t1 = clock();
    srandom(time(0));

    

    
    // initializing
    for (int i =1; i<L+1; ++i)
    {
        for ( int j = 1; j<L+1; ++j)
        {
            lat[i][j] = -1+2*(random()%2);
            
        }
    }
    rb();
 
    // Der Algorithmus läuft vollkommen gleich ab, wie die anderen für mehr Kommentierung siehe dort.
    
    // Anzahl an positiven h-Werten
    const int N_h = 7;

    double beta = 0.6;

    // externes Magnetfeld
    double h_f[2*N_h] ={-0.85,-0.3,-0.08,-0.07,-0.06,-0.05,-0.04,0.04,0.05,0.06,0.07,0.08,0.3,0.85};
    
    FILE *datei;
    datei = fopen("result4b_2.txt","w");

    FILE *datei2;
    datei2 = fopen("short_results4b_2.txt","w");

    int n =0;
    //#pragma omp parallel for
    for (int m=0; m<(4*N_h-1);++m)
    {
        
        if (m>2*N_h-1)
        n = 4*N_h-2-m;
        else
        n=m;


        // if (h_f[n]>=0.4 & h_f[n]<0.53)
        // {
        // //termalization
        // for (int t=0;t<13000;++t)
        // {
        //     sweep(beta);
        // }
        // }
        // else if (h_f[n]>=0.53 )
        // {
        // //termalization
        // for (int t=0;t<6000;++t)
        // {
        //     sweep(beta);
        // }
        // }
        // else
        // {
        // //termalization
        // for (int t=0;t<4000;++t)
        // {
        //     sweep(beta);
        // }
        // }
        
        for ( int t=0;t<5000;++t)
        {
            sweep(beta,h_f[n]);
        }
    
        // drawing phase
        double Hamilton_store[100] ={0};
        int mag_store[100] = {0};
        int iter = 100;
        int pure_ham = 0;

        for ( int N =0; N<iter; ++N)
        {
              
            pure_ham = 0;
            
            for ( int i=1; i<L+1; ++i)
            {
                for (int j =1; j<L+1; ++j)
                {
                    mag_store[N] += lat[i][j];
                    pure_ham += Ham(i,j);
                }
            }
            
            Hamilton_store[N] = (double)pure_ham/2-(double)mag_store[N]*h_f[n];
            
            
            if (mag_store[N]<0)
            mag_store[N] *= -1;

            fprintf(datei,"%.4e\t",Hamilton_store[N]);
            fprintf(datei,"%i\n",mag_store[N]);

            // if (h_f[n]>=0.4 & h_f[n]<0.53)
            // {
            // for (int t=0; t<13000;++t)
            // {
            //     sweep(beta);
            // }
            // }
            // else if (h_f[n]>=0.53)
            // {
            // for (int t=0; t<6000;++t)
            // {
            //     sweep(beta);
            // }
            // }
            // else
            // {
            for (int t=0; t<5000;++t)
            {
                sweep(beta,h_f[n]);
            }
            // }
              
    
            


            

        }
       

    
        
        
        double energy_mean= 0;
        double mag_mean =0 ;
        for (int N =0; N<iter;++N)
        {

            energy_mean += (double)Hamilton_store[N]/(double)iter;
            mag_mean += (double)mag_store[N]/(double)iter;

        }
        
       
        fprintf(datei2,"%.6e\t",energy_mean);
        fprintf(datei2,"%.6e\n",mag_mean);

        


        if ((n+1)%5==0)
        {
            printf("%i\n",n);
            t2 = clock();
            printf("%ld\n",(t2-t1)/CLOCKS_PER_SEC);
        }

        


    }
    

    fclose(datei);   
    fclose(datei2);
   


    
   

    

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
void sweep(double b,double h)
{
    //#pragma omp parallel for
    for (int i=1;i<L+1;++i)
    {
        for (int j=1;j<L+1;++j)
        {
            Hami1 = (double)Ham(i,j)-(double)lat[i][j]*h;
            
            rand_num = ((double)random())/(max);
            prob = exp((b*Hami1))/(exp(b*Hami1)+exp(-b*Hami1));
            if (rand_num<prob)
                    {
                        lat[i][j] *= -1;
                       
                    }
                    
            
                
            
                    
                
        }

        }

    rb();
    
}

