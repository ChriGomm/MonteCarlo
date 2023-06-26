#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
//#include <omp.h>

const int L = 128;
void rb();

void sweep(double b, int mh);
int Ham(int i, int j);
int lat[L+2][L+2];
double max = (double)RAND_MAX;


int main()
{

    time_t t1;
    time_t t2;
    time_t t3;
    t1 = time(NULL);
    srandom(time(0));

    
    // ausführlichere Komentierung ist in anderen Programmen zu finden.

            
    double beta[7];
    beta[0] = 0.3;
    beta[1] = 0.46;
    beta[2] = 0.49;
    beta[3] = 0.53;
    beta[4] = 0.6;
    beta[5] = 0.7;
    beta[6] = 0.8;

    FILE *datei;
    datei = fopen("therm_test_met_1l_3.txt","w");

    // multihit loop
    for (int mhits=3;mhits<7;++mhits)
    {
    for (int n=0;n<7;++n)
    {
        for (int i=1;i<L+1;++i)
    {
        for (int j=1;j<L+1;++j)
        {
            lat[i][j]= -1+2*(random()%2);
            

        }
    }

        int Hamilton_store[100000] ={0};
        int mag_store[100000] = {0};
        int iter = 100000;
        
        for ( int N =0; N<iter; ++N)
        {
            
            for ( int i=1; i<L+1; ++i)
            {
                for (int j =1; j<L+1; ++j)
                {
                    Hamilton_store[N] += Ham(i,j);
                    mag_store[N] += lat[i][j];
                }
            }
            
            



            sweep(beta[n],mhits);
            
              
    
            fprintf(datei,"%i\t",Hamilton_store[N]);
            fprintf(datei,"%i\n",mag_store[N]);

            

       



       
        

        }
    }
    }
    
    t2 = time(NULL);
            printf("%ld\n",t2-t1);

    fclose(datei);   
    
   


    
   

    

    return 0;
}


void rb()
{
    //#pragma omp parallel for
    for (int i=1;i<L+1;++i)
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
    
    H = lat[i][j]*(lat[i+1][j]+lat[i][j+1]+lat[i-1][j]+lat[i][j-1]);
    return -H;


}
void sweep(double b, int mh)
{
    //#pragma omp parallel for
    for (int i=1;i<L+1;++i)
    {
        for (int j=1;j<L+1;++j)
        {
            double Hami1 = (double) Ham(i,j);
            
            for (int mhits=0;mhits<mh;++mhits)
            {
                if (random()%2==0)
                {
                    
                    if (Hami1>=0)
                    {
                        lat[i][j] *= -1;
                        Hami1 *= -1;
                        // Hamiltonian += -Hami;
                    }
                    else
                    {
                        double rand_num = ((double)random())/((double)max);
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

