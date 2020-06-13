#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <mpi.h>

using namespace std;

int N;                 //number of bodies simulating
double position[N][3]; //positions of bodies in 3-dimension, [N][0] for x, [N][1] for y, [N][2] for z
double velocity[N][3]; //velocity of bodies in 3-dimension, [N][0] for x, [N][1] for y, [N][2] for z
double mass[N];        //mass of bodies in 3-dimension

void hermite(double position, double velocity, double mass)
{
//initial acceleration calculate
   double acceleration[N][3];
   float G = 6.67e-11;
   float r[3];
   float sp = 0.01;//softening parameter (unknown what to set yet)
   for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
	 for (int k=0; k<3; k++){
	    r[k] = abs(position[i][k]-position[j][k]);
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
            r[k] = abs(position[i][k] - position[j][k]);
	    rv[k] = abs(velocity[i][k] - velocity[j][k]);
            da[i][k] += G*mass[j]*(rv[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5)) + 3*(rv[k]*r[k])*r[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5)));
                                }
                             }
                          }

//first timestep. Use the minima timestep as the global timestep
   float t = 0;//time passed in simulation
   float dt[N];//timestep based on all particles' properties respectively.
   float tmin;//save the minimum timestep in this variable. All particles follow this timestep.
   float timescale;//the total time of simulation we want(to be given)
   for (int i=0; i<N; i++){
      dt[i] = 0.01*( abs(pow((pow(acceleration[i][0], 2) + pow(acceleration[i][1], 2) + pow(acceleration[i][2], 2)), 0.5))/abs(pow((pow(da[i][0], 2) + pow(da[i][1], 2) + pow(da[i][2], 2)), 0.5)) );
//      t[i] += dt[i];
			  }
   for (int i=0; i<N; i++){
      tmin = dt[0];
      if ( dt[i] < tmin ) tmin = dt[i];
			  }
   t += tmin;

//from now on, we want the simulation to go on by itself until it reach the time scale for simulation to complete.
while ( t<timescale ){
//predict position and velocity (to be corrected)
   for (int i=0; i<N; i++){
      for (int k=0; k<3; k++){
	 position[i][k] += ( tmin*velocity[i][k] + pow(tmin, 2)*acceleration[i][k]/2 + pow(tmin, 3)*da[i][k]/6 );
	 velocity[i][k] += ( tmin*acceleration[i][k] +  pow(tmin, 2)*da[i][k]/2 )
			     }
			  }

//update acceleration and its derivative (to be corrected)
   double af[N][3];
   double daf[N][3];
   for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
         for (int k=0; k<3; k++){
	    r[k] = abs(position[i][k]-position[j][k]);
            af[i][k] = G*mass[j]*r[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5));
	    rv[k] = abs(velocity[i][k] - velocity[j][k]);
            daf[i][k] = G*mass[j]*(rv[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5)) + 3*(rv[k]*r[k])*r[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5)));
                                }
                             }
                          }
   
//high order correction
   double a2[N][3];//second order derivative of acceleration
   double a3[N][3]; //third order derivative of acceleration
   for (int i=0; i<N; i++){
      for (int k=0; k<3; k++){
	 a2[i][k] = ( -6*(acceleration[i][k] - af[i][k]) - 2*tmin(2*da[i][k] + daf[i][k]) )/pow(tmin, 2);
	 a3[i][k] = ( 12*(acceleration[i][k] - af[i][k]) - 6*tmin(da[i][k] + daf[i][k]) )/pow(tmin, 3);
                             }
                          }                                }

//final acceleration, position and velocity
   for (int i=0; i<N; i++){
      for (int k=0; k<3; k++){
//	 acceleration[i][k] += ( tmin*da[i][k] + pow(dt[i], 2)*a2[i][k]/2 + pow(tmin, 3)*a3[i][k]/6 );
	 acceleration[i][k] = af[i][k];
	 da[i][k] = daf[i][k];
	 position[i][k] += ( pow(tmin, 4)*a2[i][k]/24 + pow(tmin, 5)*a3[i][k]/120 );
	 velocity[i][k] += ( pow(tmin, 3)*a2[i][k]/6 + pow(tmin, 4)*a3[i][k]/24 );
                             }
                          }   

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
//      t[i]   += dt[i];
			  }
   for (int i=0; i<N; i++){
      tmin = dt[0];
      if ( dt[i] < tmin ) tmin = dt[i];
                          }
   t += tmin;
//if t < given time scale, return to beginning of while loop, start another iteration.
	}
}
