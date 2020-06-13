#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <mpi.h>

using namespace std;

int main( int argc, char *argv[] )
{
const int N = 2;                 //number of bodies simulating
double position[N][3] = {
                        {10, 0, 0},
                        {-10, 0, 0}
                        }; //positions of bodies in 3-dimension, [N][0] for x, [N][1] for y, [N][2] for z
double velocity[N][3] = {
                        {-0.1, 0, 0},
                        {0.1, 0, 0}
                        }; //velocity of bodies in 3-dimension, [N][0] for x, [N][1] for y, [N][2] for z
double mass[N] = {1e10, 1e10};        //mass of bodies in 3-dimension
printf("position: x1=%lf, y1=%lf, z1=%lf, x2=%lf, y2=%lf, z2=%lf\n", position[0][0], position[0][1], position[0][2], position[1][0], position[1][1], position[1][2]);
double acceleration[N][3];
float G = 6.67e-11;
float r[3];
float sp = 0.01;//softening parameter
   for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
         for (int k=0; k<3; k++){
            r[k] = fabs(position[i][k]-position[j][k]);
            acceleration[i][k] += G*mass[j]*r[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5));
                                }
                             }
                          }

//initial acceleration derivative of time calculate
   double da[N][3];
   float rv[3];//relative speed between bodies
   for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
         for (int k=0; k<3; k++){
            r[k] = fabs(position[i][k] - position[j][k]);
            rv[k] = fabs(velocity[i][k] - velocity[j][k]);
            da[i][k] += G*mass[j]*(rv[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5)) + 3*(rv[k]*r[k])*r[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5)));
                                }
                             }
                          }

//first timestep. Use the minima timestep as the global timestep
float t = 0;//time passed in simulation
float dt[N];//timestep based on all particles' properties respectively.
float tmin;//save the minimum timestep in this variable. All particles follow this timestep.
float timescale = 50;//the total time of simulation we want(to be given)
for (int i=0; i<N; i++){
      dt[i] = 0.01*( fabs(pow((pow(acceleration[i][0], 2) + pow(acceleration[i][1], 2) + pow(acceleration[i][2], 2)), 0.5))/fabs(pow((pow(da[i][0], 2) + pow(da[i][1], 2) + pow(da[i][2], 2)), 0.5)) );
//      t[i] += dt[i];
                       }
tmin = dt[0];
for (int i=0; i<N; i++){
   if ( dt[i] < tmin ) tmin = dt[i];
                       }
   printf("%f\n", tmin);
   t += tmin;
double af[N][3];
double daf[N][3];
double a2[N][3];//second order derivative of acceleration
double a3[N][3]; //third order derivative of acceleration

while ( t<timescale )
{
//predict position and velocity (to be corrected)
   for (int i=0; i<N; i++){
   for (int k=0; k<3; k++){
      position[i][k] += ( tmin*velocity[i][k] + pow(tmin, 2)*acceleration[i][k]/2 + pow(tmin, 3)*da[i][k]/6 );
      velocity[i][k] += ( tmin*acceleration[i][k] +  pow(tmin, 2)*da[i][k]/2 );
                          }
                          }

//update acceleration and its derivative
   for (int i=0; i<N; i++){
   for (int j=0; j<N; j++){
   for (int k=0; k<3; k++){
      r[k] = fabs(position[i][k]-position[j][k]);
      af[i][k] = G*mass[j]*r[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5));
      rv[k] = fabs(velocity[i][k] - velocity[j][k]);
      daf[i][k] = G*mass[j]*(rv[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5)) + 3*(rv[k]*r[k])*r[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5)));
                          }
                          }
                          }

//high order correction
   for (int i=0; i<N; i++){
   for (int k=0; k<3; k++){
      a2[i][k] = ( -6*(acceleration[i][k] - af[i][k]) - 2*tmin*(2*da[i][k] + daf[i][k]) )/pow(tmin, 2);
      a3[i][k] = ( 12*(acceleration[i][k] - af[i][k]) - 6*tmin*(da[i][k] + daf[i][k]) )/pow(tmin, 3);
                          }
                          }

//final acceleration, position and velocity
   for (int i=0; i<N; i++){
   for (int k=0; k<3; k++){
      acceleration[i][k] = af[i][k];
      da[i][k] = daf[i][k];
      position[i][k] += ( pow(tmin, 4)*a2[i][k]/24 + pow(tmin, 5)*a3[i][k]/120 );
      velocity[i][k] += ( pow(tmin, 3)*a2[i][k]/6 + pow(tmin, 4)*a3[i][k]/24 );
                          }
                          }
//check the result
 printf("position: x1=%lf, y1=%lf, z1=%lf, x2=%lf, y2=%lf, z2=%lf\n", position[0][0], position[0][1], position[0][2], position[1][0], position[1][1], position[1][2]);

//new timestep to next step
double aaf[N];
double aaf2[N];
double adaf[N];
double aa3[N];
   for (int i=0; i<N; i++){
      aaf[i]  = pow((pow(af[i][0], 2) + pow(af[i][1], 2) + pow(af[i][2], 2)), 0.5);
      adaf[i] = pow((pow(daf[i][0], 2) + pow(daf[i][1], 2) + pow(daf[i][2], 2)), 0.5);
      aa3[i]  = pow((pow(a3[i][0], 2) + pow(a3[i][1], 2) + pow(a3[i][2], 2)), 0.5);
      aaf2[i] = pow((pow((a2[i][0] + dt[i]*a3[i][0]), 2) + pow((a2[i][1] + dt[i]*a3[i][1]), 2) + pow((a2[i][2] + dt[i]*a3[i][2]), 2)), 0.5);
      dt[i]   = pow(((aaf[i]*aaf2[i] + pow(adaf[i], 2))/(adaf[i]*aa3[i] + pow(aaf2[i], 2))), 0.5);
			  }
   tmin = dt[0];
   for (int i=0; i<N; i++){
      if ( dt[i] < tmin ) tmin = dt[i];
                          }
   printf("%f\n", tmin);
   t += tmin;

}
   return EXIT_SUCCESS;
}
