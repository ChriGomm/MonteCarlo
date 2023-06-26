#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
// #include <omp.h>

// Thermalisierungsuntersuchung für Abhängigkeit von externem Magnetfeld für den Heat Bath Algorithmus.

const int L = 128;
void rb();

void sweep(double b, double h);
int Ham(int i, int j);
int lat[L+2][L+2] = {0};

double max = (double)RAND_MAX;




int main()
{

    clock_t t1;
    clock_t t2;
    time_t t3;
    t1 = clock();
    srandom(time(0));

    
    
    double beta[3];
    beta[0] = 0.5;
    beta[1] = 0.6;
    beta[2] = 0.8;
    // beta[3] = 0.53;
    // beta[4] = 0.6;
    // beta[5] = 0.7;
    // beta[6] = 0.8;
    double h_f[7] ={0.04,0.05,0.06,0.07,0.08,0.3,0.85};
    // for (int i =0; i<30; ++i)
    // h_f[i+2]= 0.05 +0.01/25*(double)i;
    // for (int i =0;i<2;++i)
    // h_f[i] = 0.01+0.02*(double)i;

    // for (int i=0;i<3;++i)
    // h_f[i+32]= 0.44+ 0.1*(double)i;

    // for ( int i=0;i<7;++i)
    // h_f[i]=0.01/7*(double)(i+1);



    FILE *datei;
    datei = fopen("therm_test_hb_h.txt","w");

    
    //#pragma omp parallel for
    // Temperatur Schleife
    for (int it=0;it<3;++it)
    {
        // Feldschleife
    for (int n=0; n<7;++n)
    {
        
        
        for (int i=1;i<L+1;++i)
    {
        for (int j=1;j<L+1;++j)
        {
            lat[i][j]= -1+2*(random()%2);
            
        }
    }
    rb();   
    
        // drawing phase
        double Hamilton_store=0;
        int mag_store=0;
        int mag_partial;
        int iter = 100000;
        
        int pure_ham;

        for ( int N =0; N<iter; ++N)
        {
             
            
            
    
            pure_ham = 0;
            Hamilton_store = 0;
            mag_store = 0;
            
            for ( int i=0; i<L; ++i)
            {
                for (int j =0; j<L; ++j)
                {
                    mag_store += lat[i][j];
                    pure_ham += Ham(i,j);
                }
            }
            
            Hamilton_store = (double)pure_ham-(double)mag_store*h_f[n];
            
            

            fprintf(datei,"%.2f\t",Hamilton_store);
            fprintf(datei,"%i\n",mag_store);
        
        
            sweep(beta[it],h_f[n]);
        
        
        
        
        

        
        
              
    
            


            

        }
       

    
        
        
        
        

        


    }}
    

           fclose(datei); 
            t2 = clock();
            printf("%ld\n",(t2-t1)/CLOCKS_PER_SEC);
        
      
   
    

    
   

    

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
    double Hami1;
    double rand_num;
    double prob;
    
    for (int i=1;i<L+1;++i)
    {
        for (int j=1;j<L+1;++j)
        {
            Hami1 = (double)Ham(i,j)-(double)lat[i][j]*h;

            
            rand_num = ((double)random())/((double)max);
            prob = exp((b*Hami1))/(exp(b*Hami1)+exp(-b*Hami1));
            if (rand_num<prob)
                    {
                        lat[i][j] *= -1;
                       
                    }
                    
            
        }
    }
    rb();
}

