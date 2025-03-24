#include <stdio.h>
#include "RKF78.h"
#include <math.h>

/* Es compila com:
 * gcc -Wall -O3 -o Test_RKF78Sys.exe Test_RKF78Sys.c RKF78.c -lm
 */

/* Es potposar el contingut que es vulgui i el nom que es vulgui per exemple ODE_Sys_Parameters --> Manuela */
typedef struct {
    double lambda;
    double Theta;
    double drift;
} ODE_Sys_Parameters;

/* The equation describing the motion of a pendulum can be described by
 * the single dependent variable $\theta$ representing the angle the
 * pendulum makes with the vertical.
 * The coefficients of the equation depend on the length of the pendulum,
 * the mass of the bob, and the gravitational constant.
 * Assuming a coefficient value of k, the equation is
                $\theta'' + k\theta = 0$
 * and one possible set of initial conditions is
        $\theta(t_{0}) = 1,$
        $\theta'(t_{0}) = 0.$
 * This second order equation can be written as a system
        y'_1 = y_2,
        y'_2 = -k * y_1

        y_1(0) = 1, y_2(0) = 0
 * (Recall that this transformation is accomplished by the change of variables $y_1=\theta$ and $y_2=\theta'$.)
 *
 * The solution for these initial conditions is
        \theta(t) = y_1(t) = cos(sqrt(k) * t)
*/

void Example60D(double t, double y[], unsigned n, double f[], void *Params){
    f[0] = y[1];   f[1] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[0];
    f[2] = y[3];   f[3] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[2];
    f[4] = y[5];   f[5] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[4];
    f[6] = y[7];   f[7] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[6];
    f[8] = y[9];   f[9] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[8];
    f[10] = y[11]; f[11] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[10];
    f[12] = y[13]; f[13] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[12];
    f[14] = y[15]; f[15] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[14];
    f[16] = y[17]; f[17] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[16];
    f[18] = y[19]; f[19] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[18];
    f[20] = y[21]; f[21] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[20];
    f[22] = y[23]; f[23] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[22];
    f[24] = y[25]; f[25] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[24];
    f[26] = y[27]; f[27] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[26];
    f[28] = y[29]; f[29] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[28];
    f[30] = y[31]; f[31] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[30];
    f[32] = y[33]; f[33] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[32];
    f[34] = y[35]; f[35] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[34];
    f[36] = y[37]; f[37] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[36];
    f[38] = y[39]; f[39] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[38];
    f[40] = y[41]; f[41] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[40];
    f[42] = y[43]; f[43] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[42];
    f[44] = y[45]; f[45] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[44];
    f[46] = y[47]; f[47] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[46];
    f[48] = y[49]; f[49] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[48];
    f[50] = y[51]; f[51] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[50];
    f[52] = y[53]; f[53] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[52];
    f[54] = y[55]; f[55] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[54];
    f[56] = y[57]; f[57] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[56];
    f[58] = y[59]; f[59] = -(((ODE_Sys_Parameters *) Params)->lambda)*y[58];
}

int main ( ) {
    int tfin;
    double hmin = 1.e-8, hmax = 1.e-3, h = 1.e-6, err;
    double t=0.0, y[60] = { 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                            1.0, 0.0, 1.0, 0.0 };

    ODE_Sys_Parameters Example60D_Parameters = {2.5, 1.e12, 7.0};

    printf ( "#### Test the RKF78 function (C version)\n" );
    printf ( "####   Solve 30 times the vector ODE:\n" );
    printf ( "####      x'_1 =  x_2\n" );
    printf ( "####      x'_2 = -kx_1\n" );
    printf ( "####   This system is equivalent to the second order system z\" = - kz.\n" );


    for(tfin=1; tfin < 101; tfin++) {
		int retval;     
		double kk = sqrt(Example60D_Parameters.lambda), sol;
		register int i;

		while(t+h < tfin) { retval = RKF78Sys(&t, y, 60U, &h, &err, hmin, hmax, 1.e-15, &Example60D_Parameters, Example60D);
			sol = cos(kk*t);
			for(i=0; i < 60; i += 2) if(fabs(y[i] - sol) > 2.e-10) printf("ALARM: %20.15lf %d %20.15lf %20.15lf %.15g (%d) %.15g\n", t, i, y[i], y[i+1], err, retval, fabs(y[i]-sol));
//		printf("%20.15lf %20.15lf %20.15lf %.15g (%d) %.15g\n", t, y[0], y[1], err, retval, fabs(y[0]-sol));
		}
		h = tfin - t;
		retval = RKF78Sys(&t, y, 60U, &h, &err, hmin, hmax, 1.e-15, &Example60D_Parameters, Example60D);
		sol = cos(kk*t);
		for(i=0; i < 60; i += 2) if(fabs(y[i] - sol) > 2.e-10) printf("ALARM: %20.15lf %d %20.15lf %20.15lf %.15g (%d) %.15g\n", t, i, y[i], y[i+1], err, retval, fabs(y[i]-sol));
		printf("%20.15lf %20.15lf %20.15lf %.15g (%d) %.15g\n", t, y[0], y[1], err, retval, fabs(y[0]-sol));
    }
    return 0;
}
