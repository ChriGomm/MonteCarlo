#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
// #include <omp.h>

// Thermalisierungsuntersuchung für den Heat Bath Algorithmus ohne externes Magnetfeld
// genauere Kommentieung ist in anderen Porgrammen zu finden.


const int L = 128;
void rb();

void sweep(double b, double h);
int Ham(int i, int j);
int lat[L+2][L+2] = {0};

double max = (double)RAND_MAX;
double Hami1;
double rand_num;
double prob;



int main()
{

    clock_t t1;
    clock_t t2;
    time_t t3;
    t1 = clock();
    srandom(time(0));

    
    for (int i=1;i<L+1;++i)
    {
        for (int j=1;j<L+1;++j)
        {
            lat[i][j]= -1+2*(random()%2);
            

        }
    }
    
    double beta[7];
    beta[0] = 0.55;
    beta[1] = 0.57;
    beta[2] = 0.49;
    beta[3] = 0.53;
    beta[4] = 0.6;
    beta[5] = 0.7;
    beta[6] = 0.8;
    double h_f[7] ={0};
    

    



    FILE *datei;
    datei = fopen("therm_test_hb_woh.txt","w");

    
    //#pragma omp parallel for
    for (int n=0; n<7;++n)
    {
        
        
        for (int i=1;i<L+1;++i)
    {
        for (int j=1;j<L+1;++j)
        {
            lat[i][j]= -1+2*(random()%2);
            
        }
    }
        
    
        // drawing phase
        
        int mag_store=0;
        
        int iter = 100000;
        
        int pure_ham;

        for ( int N =0; N<iter; ++N)
        {
            
            
            
    
            pure_ham = 0;
            
            mag_store = 0;
            
            for ( int i=0; i<L; ++i)
            {
                for (int j =0; j<L; ++j)
                {
                    mag_store += lat[i][j];
                    pure_ham += Ham(i,j);
                }
            }
            
            
            fprintf(datei,"%i\t",pure_ham);
            fprintf(datei,"%i\n",mag_store);
        
            sweep(beta[n],0);
        
        
              
    
            


            

        }
       

    
        
         
        
        

        


    }
    

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

