#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
//#include <omp.h>

const int L = 128;
// boundary conditions implementation
void rb();

void sweep(double b);
int Ham(int i, int j);
// da omp für mich nicht funktioniert hat, habe ich das Programm wieder auf ein Gitter umgeschrieben um möglicherweise Zeit zu sparen.
int lat[L+2][L+2] = {0};

// maximalwert des random number generators in C
double max = (double)RAND_MAX;
double rand_num;
double prob;
double Hami1;



int main()
{
    // Zeitmessung 
    clock_t t1;
   clock_t t2;
    time_t t3;
    t1 = clock();
    srandom(time(0));

    
    
    // Initialisierung der Temperatur
    double beta[35] ={0};
    for (int i =0; i<30; ++i)
    beta[i+2]= 0.15 +0.6/30*(double)i;
    for (int i =0;i<2;++i)
    beta[i] = 0.1+0.04*(double)i;

    for (int i=0;i<3;++i)
    beta[i+32]= 0.8+ 0.1*(double)i;

    // Datei in die alle Messwerte geschrieben werden.
    FILE *datei;
    datei = fopen("result4a_4.txt","w");

    // Datei mit den relavanten Mittelwerten von eps und |m|
    FILE *datei2;
    datei2 = fopen("short_results4a_4.txt","w");

    // hot start, außerhalb der Temperaturschleife ergibt glattere Resutlate
    for (int i=1;i<L+1;++i)
    {
        for (int j=1;j<L+1;++j)
        {
            lat[i][j]= -1+2*(random()%2);
            

        }
    }
    rb();
    //#pragma omp parallel for
    for (int n=0; n<35;++n)
    {
    
    // for (int i=0;i<L;++i)
    // {
    //     for (int j=0;j<L;++j)
    //     {
    //         lat[i][j]= -1+2*(random()%2);
            

    //     }
    // }
        // 
        // Termalization
        if (beta[n]>=0.4 & beta[n]<0.53)
        {
        //termalization
        for (int t=0;t<13000;++t)
        {
            sweep(beta[n]);
        }
        }
        else if (0.6>=beta[n]>=0.53 )
        {
        //termalization
        for (int t=0;t<13000;++t)
        {
            sweep(beta[n]);
        }
        }
        else
        {
        //termalization
        for (int t=0;t<10000;++t)
        {
            sweep(beta[n]);
        }
        }
        
        
    
        // drawing phase
        
        int Hamilton_store[100] ={0};
        int mag_store[100] = {0};
        int iter = 100;
        
        for ( int N =0; N<iter; ++N)
        {
            
            for ( int i=1; i<L+1; ++i)
            {
                for (int j =1; j<L+1; ++j)
                {
                    // Gesamthamiltonian wird berechnet
                    Hamilton_store[N] += Ham(i,j);
                    // Gesamtmagnetisierung wird berechnet
                    mag_store[N] += lat[i][j];
                }
            }
            
            // betrag der Magnetisierung wird berechnet
            if (mag_store[N]<0)
            mag_store[N] *= -1;

            // gemessene Werte werden in Datei geschrieben
            fprintf(datei,"%i\t",Hamilton_store[N]);
            fprintf(datei,"%i\n",mag_store[N]);

            // Dekorrelierung
            if (beta[n]>=0.4 & beta[n]<0.53)
            {
            for (int t=0; t<13000;++t)
            {
                sweep(beta[n]);
            }
            }
            else if (beta[n]>=0.53)
            {
            for (int t=0; t<13000;++t)
            {
                sweep(beta[n]);
            }
            }
            else
            {
            for (int t=0; t<10000;++t)
            {
                sweep(beta[n]);
            }
            }
              
    
            


            

        }
       

    
        
        // Berechnung der Mittelwerte
        double energy_mean= 0;
        double mag_mean =0 ;
        for (int N =0; N<iter;++N)
        {

            energy_mean += (double)Hamilton_store[N]/2/(double)iter;
            mag_mean += (double)mag_store[N]/(double)iter;

        }
        
    //    Mittelwerte werden in datei geschrieben
        fprintf(datei2,"%.6e\t",energy_mean);
        fprintf(datei2,"%.6e\n",mag_mean);

        

    // Zwischenstand wird ausgegeben 
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
void sweep(double b)
{
    //#pragma omp parallel for
    for (int i=1;i<L+1;++i)
    {
        for (int j=1;j<L+1;++j)
        {
            Hami1 = (double)Ham(i,j);
            
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

