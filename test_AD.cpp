#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float ***buildGrid(int numRows, int numCols, int numLevels);
void acceleration_deposition( int N, float ***a_grid, double *x, double *y, double *z, double gs, int GN, int mode, double *a);

int main(void){
    double gs = 1.0;
    int GN = 5;
    int N = 1;
    double M[N], x[N], y[N], z[N];
    M[1] = 10.0;
    x[1] = 1.3;
    y[1] = 2.2;
    z[1] = 0.8;
    int mode = 3;
    float ***a_grid = buildGrid(GN,GN,GN);
    // Note we need to allocate a place to pass the 3-d function.
    double a[N];
    
     for(int i = 0; i<GN; i++){
         for(int j = 0; j<GN; j++){
             for(int k = 0; k<GN; k++){
                 a_grid[i][j][k] = 1.0;// provides the values a where we should give the value of real a when implementation.
             }
         }
     }

    acceleration_deposition(N, a_grid, x, y, z, gs, GN, mode, a);
    
    for(int n = 0; n<N; n++){
        printf( "a[%2d] = %5.3f \n",n, a[n] );
    }
    return 0;
}

float ***buildGrid(int numRows, int numCols, int numLevels)
{
    float ***levels;
    levels = (float * * *)malloc(numLevels *sizeof(float *)); //Contains all levels

    int rowIndex, levelIndex;

    for (levelIndex = 0; levelIndex < numLevels; levelIndex++)
    {
        float **level = (float * *)malloc(numRows * sizeof(float *)); //Contains all rows

        for(rowIndex = 0; rowIndex < numRows; rowIndex++)
        {
            level[rowIndex] = (float *)malloc(numCols * sizeof(float)); //Contains all columns
        }

        levels[levelIndex] = level;
    }

    return levels;
}

void acceleration_deposition( int N, float ***a_grid, double *x, double *y, double *z, double gs, int GN, int mode, double *a)
{
/*
 mode 1: NGP
 mode 2: CIC
 mode 3: TSC
 N is total number of particles.
 gs is grid size.
 GN is total grid number. Note that I assume we have a equilateral box.
 a_grid is the input acceleration matrix.
 a is output acceleration of one component for every particle. (a zero matrix pointer)
*/
    for(int n = 0; n<N; n++){
        a[n] = 0.0;
    }
    double dx, dy, dz;
    if(mode == 1)
    {
        for(int n = 0; n<N; n++)
        {
            for(int i = 0; i<GN; i++)
            {
                dx = fabs(x[n]-i*gs);
                for(int j = 0; j<GN; j++)
                {
                    dy = fabs(y[n]-j*gs);
                    for(int k = 0; k<GN; k++)
                    {
                        dz = fabs(z[n]-k*gs);
                        if(dx<=0.5*gs && dy<=0.5*gs && dz <=0.5*gs)
                        {
                           a[n] += a_grid[i][j][k];
                        }
                    }
                }
            }
        }
    }
    if(mode == 2)
    {
        for(int n = 0; n<N; n++)
        {
            for(int i = 0; i<GN; i++)
            {
                dx = fabs(x[n]-i*gs);
                for(int j = 0; j<GN; j++)
                {
                    dy = fabs(y[n]-j*gs);
                    for(int k = 0; k<GN; k++)
                    {
                        dz = fabs(z[n]-k*gs);
                        if(dx<=gs && dy<=gs && dz<=gs)
                        {
                            a[n] += (1.0-dx/gs)*(1.0-dy/gs)*(1.0-dz/gs)*a_grid[i][j][k];
                        }
                    }
                }
            }
        }
    }
    if(mode == 3)
    {
        double wx, wy, wz; //weighted function
        double ai; //
        for(int n = 0; n<N; n++)
        {
            for(int i = 0; i<GN; i++)
            {
                dx = fabs(x[n]-i*gs);
                if(dx<=0.5*gs)
                {
                    wx = 0.75-dx*dx/gs/gs;
                }
                else if(dx>=0.5*gs && dx <= 1.5*gs)
                {
                    wx = 0.5*(1.5-dx/gs)*(1.5-dx/gs);
                }
                else wx = 0.0;
                
                for(int j = 0; j<GN; j++)
                {
                    dy = fabs(y[n]-j*gs);
                    if(dy<=0.5*gs)
                    {
                        wy = 0.75-dy*dy/gs/gs;
                    }
                    else if(dy>=0.5*gs && dy <= 1.5*gs)
                    {
                        wy = 0.5*(1.5-dy/gs)*(1.5-dy/gs);
                    }
                    else wy = 0.0;
                    
                    for(int k = 0; k<=GN; k++)
                    {
                        dz = fabs(z[n]-k*gs);
                        if(dz<=0.5*gs)
                        {
                            wz = 0.75-dz*dz/gs/gs;
                        }
                        else if(dz>=0.5*gs && dz <= 1.5*gs)
                        {
                            wz = 0.5*(1.5-dz/gs)*(1.5-dz/gs);
                        }
                        else wz = 0.0;
                        
                        a[n] += wx*wy*wz*a_grid[i][j][k];
                        
                    }
                }
            }
        }
    }
}

