#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <fftw3.h>
#include "Propagator.h"
#include "SCFT.h"

#define PARAMETER_SET

//Simulates single chain of homopolymer.
const double N = 1.0; // Chain length
const double q_0 = 1.0;

#ifdef PARAMETER_SET
//Fredrickson.
const int N_s = 512; // Number of integration steps along the chain
const int M_x = 512; // Grid points in x direction
const double R_g =  15 ; // chain gyration radius (unperturbed radius of gyration)
const double b = sqrt(6) * R_g / sqrt(N);//Kuhn segment length
const double L = 10.0 * R_g; // box length
#else
//Chen,Kim and Alexander-Kaiz
const int N_s = 2001; // Number of integration steps along the chain
const int M_x = 256; // Grid points in x direction
const double R_g =  1 ; // chain gyration radius (unperturbed radius of gyration)
const double b = sqrt(6) * R_g / sqrt(N);//Kuhn segment length
const double L = 3.2 * R_g; // box length
#endif
const double delta_s = N/(N_s); // chain integration step size
const double delta_x = L / M_x;
/*
 * Mathematical identity
 */

inline double sech(double x) {
    return 1.0 / cosh(x);
}

/*
 * Potential function w(x)
 */
/*
inline double w(double x) {
#ifdef PARAMETER_SET
    return (1.0  - 2.0 * pow(sech(3.0 * (x - L / 2.0)/ (2 * R_g)),2.0)) / N;
#else
    //return 10.0 * sin(2 * M_PI * (x - L / 2.0) / L) / N;
    double halfway = delta_x * M_x / 2.0;
    const double xi = 1.0;
    if(x < halfway) {
        return 10.0 * tanh((x - M_x / 4.0 * delta_x) / xi);
    } else {
        return - 10.0 * tanh((x - 3.0 * M_x / 4.0 * delta_x) / xi);
    }
#endif
}*/

int main(int argc, char** argv) {
    //xN=13, hence flory huggins x = 13/N.
    SCFT propagator(M_x,N_s,L,N,R_g,0.5, 13.0 / N, 0.2, SimpleMixing);
    propagator.Run();
    propagator.Save("hell.csv");
    return 0;
}