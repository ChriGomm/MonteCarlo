#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
//#include <omp.h>

const int L = 32;
// boundary conditions implementation
void rb();
void spinScreener(int visited[32][32],int cluster[32][32], double exponential, int start_x, int start_y, int ref);
void sweep(double exponential);
int Ham(int i, int j);
// da omp für mich nicht funktioniert hat, habe ich das Programm wieder auf ein Gitter umgeschrieben um möglicherweise Zeit zu sparen.
int lat[32][32] = {0};

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
    datei = fopen("mag_wolff.txt","w");

    // Datei mit den relavanten Mittelwerten von eps und |m|
    FILE *datei2;
    datei2 = fopen("ham_wolff.txt","w");

    // FILE *cluster;
    // cluster = fopen("suveillance.txt","w");
    // hot start, außerhalb der Temperaturschleife ergibt glattere Resutlate
    
    // rb();
    //#pragma omp parallel f
    for (int n=0; n<35;++n)
    {
        if (n>=13){
    double expo = exp(-2*beta[n]);
   
   for (int i=0;i<L;++i)
    { 
        for (int j=0;j<L;++j)
        {
            lat[i][j]= -1+2*(random()%2);
            

        }
    }
        int Hamilton_store = 0;
        int mag_store = 0;
        int iter = 400;
        
        for ( int N =0; N<iter; ++N)
        {
            
            Hamilton_store =0;
            mag_store = 0;
            for ( int i=0; i<L; ++i)
            {
                for (int j =0; j<L; ++j)
                {

                    // Gesamthamiltonian wird berechnet
                    Hamilton_store += Ham(i,j);
                    // Gesamtmagnetisierung wird berechnet
                    mag_store += lat[i][j];
                }
            }
            sweep(expo);
            // betrag der Magnetisierung wird berechnet
            // if (mag_store<0)
            // mag_store *= -1;

            // gemessene Werte werden in Datei geschrieben
            fprintf(datei2,"%i\t",Hamilton_store);
            fprintf(datei,"%i\t",mag_store);

        }
       

    fprintf(datei2,"\n");
    fprintf(datei,"\n");
        
       
    // Zwischenstand wird ausgegeben 
        if ((n+1)%5==0)
        {
            printf("%i\n",n);
            t2 = clock();
            printf("%ld\n",(t2-t1)/CLOCKS_PER_SEC);
        }

        
        }

    }
    

    fclose(datei);   
    fclose(datei2);
    // fclose(cluster);
    return 0;
}


// void rb()
// {
//     //#pragma omp parallel for
//    for (int i=0;i<L+2;++i)
//     {
//         lat[i][L+1] = lat[i][1];
//         lat[L+1][i] = lat[1][i];
//         lat[i][0] = lat[i][L];
//         lat[0][i] = lat[L][i];
//     }
    
    

// }
int Ham(int i, int j)
{
    int H;
    H = lat[i][j]*(lat[(L+i-1)%L][j]+lat[i][(j+1)%L]+lat[i][(L+j-1)%L]+lat[(i+1)%L][j]);
    return -H;


}

void spinScreener(int visited[32][32],int cluster[32][32], double exponential,int start_x,int start_y,int ref){
    if (visited[start_x][start_y]==0){
        visited[start_x][start_y] = 1;
        if (lat[start_x][start_y]==-1){
            cluster[start_x][start_y] = 1;
            double rand_num = (double)rand();
            if (rand_num<=(1-(exponential))*max){
                lat[start_x][start_y] *= 1;
            }
            spinScreener(visited,cluster, exponential,(start_x+1)%(L),start_y,ref);
            spinScreener(visited,cluster, exponential,(L+start_x-1)%(L),start_y,ref);
            spinScreener(visited,cluster, exponential,start_x,(start_y+1)%(L),ref);
            spinScreener(visited,cluster, exponential,start_x,(L+start_y-1)%(L),ref);
    }}
    }


void sweep(double exponential)
{
    int visited[32][32] = {0};
    // for (int i=0;i<L;i++){
    //     printf("%i\t",visited[10][i]);
    // }    
    // printf("\n");
    int rand_x = (int)(((double)random())/(max)*L);
    int rand_y = (int)(((double)random())/(max)*L);
    int cluster[32][32];
    int ref = lat[rand_x][rand_y];
    spinScreener(visited,cluster, exponential,rand_x,rand_y,ref);
    for (int i =0;i<L;i++){
        for (int j=0;j<L;j++){
            if (cluster[i][j]==0){
                printf("\033[0;37m");
                if (lat[i][j]==1){
                    printf("+1 ");
                }
                else{
                printf("%i ",lat[i][j]);
                }
                               
            }
            else if (i==rand_x && j==rand_y){
                printf("\033[0;31m");
                if (lat[i][j]==1){
                    printf("+1 ");
                }
                else{
                printf("%i ",lat[i][j]);
                }
            }
            else{
                printf("\033[0;32m");
                if (lat[i][j]==1){
                    printf("+1 ");
                }
                else{
                printf("%i ",lat[i][j]);
                }
            }
        }
        printf("\n");
    }
    printf("\n");
}

