#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
//#include <omp.h>

// Metropolis Multihit Algorithmus

const int L = 128;
void rb();
void sweep(double b, double expon[]);
int Ham(int i, int j);
int lat[130][130] = {0};

double max = (double)RAND_MAX; 
const int mh = 1;
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
    FILE *datei_ham;
    datei_ham = fopen("Ham_store.txt","w");

    // Datei mit Mittelwerten
    FILE *datei_mag;
    datei_mag = fopen("Mag_store.txt","w");

    //#pragma omp parallel for
    for (int n=0; n<35;++n)
    {
        if (n== 18 || n==19 || n== 24 || n==33){
        double expo[5] = {exp(beta[n]*(-8)),exp(beta[n]*(-4)),0,exp(beta[n]*(4)),exp(beta[n]*(8))};
        // printf("%f\n",beta[n]);
        // for (int it=0;it<5;it++)
        // printf("%f\t",expo[it]);
        // printf("\n");       
        if (beta[n]<0.44){

        
        
        for (int i=1;i<L+1;++i)
    {
        for (int j=1;j<L+1;++j)
        {
            lat[i][j]= 1;
            

        }
    }}
    else {
        for (int i=1;i<L+1;++i)
    {
        for (int j=1;j<L+1;++j)
        {
            lat[i][j]= -1+2*(random()%2);
            

        }
    }
    
    }
    
    
        // Arrays in die die Messungen geschrieben werden Anzahl Messungen= iter
        int Hamilton_store =0 ;
        int mag_store = 0;
        int mag_sq_store = 0;
        int iter = 300000;
        
        for ( int N =0; N<iter; ++N)
        {
           Hamilton_store = 0;
           mag_store = 0; 
            for ( int i=1; i<L+1; ++i)
            {
                for (int j =1; j<L+1; ++j)
                {
                    // Gesamthamiltonian
                    Hamilton_store += Ham(i,j);
                    // Gesamtmagnetisierung
                    mag_store += lat[i][j];
                }
            }
            
            
            if (mag_store<0)
            mag_store *= -1;
         
              
            sweep(beta[n],expo);
            fprintf(datei_ham,"%i\t",Hamilton_store);
            fprintf(datei_mag,"%i\t",mag_store);


            

        }
       

            fprintf(datei_ham,"\n");
            fprintf(datei_mag,"\n");
        
        }

       

        

        // Zwischenergebnis
        if ((n+1)%5==0)
        {
            printf("%i\n",n);
            t2 = time(NULL);
            printf("%ld\n",t2-t1);
        }

        


    }
    

    fclose(datei_ham);   
    fclose(datei_mag);
    t2=time(NULL);
    printf("%ld\n",t2-t1);


    
   

    

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

void sweep(double b, double expon[])
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
                        if (expon[(int)(2+Hami1/2)]>=rand_num)
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

