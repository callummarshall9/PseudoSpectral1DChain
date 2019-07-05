#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <fftw3.h>
#define PARAMETER_SET

//Simulates single chain of homopolymer.

#ifdef PARAMETER_SET
//Fredrickson.
const double N = 10; // Chain length
const int N_s = 512; // Number of integration steps along the chain
const double delta_s = N/(N_s); // chain integration step size
const int M_x = 512; // Grid points in x direction
const double q_0 = 1.0;
const double R_g =  15 ; // chain gyration radius (unperturbed radius of gyration)
const double b = sqrt(6) * R_g / sqrt(N);//Kuhn segment length
const double L = 10.0 * R_g; // box length
const double delta_x = L / M_x;
#else
//Chen,Kim and Alexander-Kaiz
const double N = 10; // Chain length
const int N_s = 2001; // Number of integration steps along the chain
const double delta_s = N/(N_s); // chain integration step size
const int M_x = 256; // Grid points in x direction
const double q_0 = 1.0;
const double R_g =  1 ; // chain gyration radius (unperturbed radius of gyration)
const double b = sqrt(6) * R_g / sqrt(N);//Kuhn segment length
const double L = 3.2 * R_g; // box length
const double delta_x = L / M_x;
#endif
/*
 * Mathematical identity
 */

inline double sech(double x) {
    return 1.0 / cosh(x);
}

/*
 * Potential function w(x)
 */

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
}

/*
 * Designed to increment q^(n) to q^(n+1/3) for two steps in the Psuedo spectral method.
 */


void apply_q_operator(double q[M_x][N_s+1], int src_N_s, int dest_N_s) {
    for(int j = 0; j < M_x; j++) {
        double x = delta_x * (double)j;
        //Source: https://en.wikipedia.org/wiki/Pseudo-spectral_method
        q[j][dest_N_s] = q[j][src_N_s] * exp(-1.0 * w(x)   * delta_s  / 2);
        //Fig: 3.14, divided by N, as the potential function is scaled Nw(x)=function in w.
        //This leads to the unscaled potential function is therefore w(x)/N.
    }
}

/*
 * The successful method of 3 ways.
 */

void another_way(double q[M_x][N_s+1], double &Q) {
        //double x = L / (M_x+1) * j;

    //Data structures to hold the FFT transform data.

    fftw_complex* a_j = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * M_x);
    fftw_complex* h_j = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * M_x);
    fftw_complex* q_complex = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * M_x);
    fftw_plan p1 = fftw_plan_dft_1d(M_x, q_complex, a_j, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan p2 = fftw_plan_dft_1d(M_x, h_j, q_complex, FFTW_BACKWARD, FFTW_ESTIMATE);
    //Populate Q_{M_x}(0) according to Pseudo Spectral Algorithm 3.1: Continuous Gaussian Chain Step 1.
    for(int m = 0; m < M_x; m++) {
        q[m][0] = 1;
    }
    double s = delta_s;
    for(int n = 0; n < N_s; n++) {

        //Initialise M column vector.
        //Calculate next step and get q^(n+1/3) according to Pseudo Spectral Algorithm 3.1: Continuous Gaussian Chain Step 2.
        apply_q_operator(q, n,n+1);//Populate current Q from previous Q.

        for(int i = 0; i < M_x; i++) {
            q_complex[i][0] = q[i][n+1];
            q_complex[i][1] = 0;
        }
        //Populate the q complex array.
        fftw_execute(p1);
        if(n == N_s - 1) {
            Q = a_j[0][0] / M_x;
        }
        //Performed the FFT according to Pseudo Spectral Algorithm 3.1: Continuous Gaussian Chain Step 3.

        //Calculate the h_j coefficients according to Pseudo Spectral Algorithm 3.1: Continuous Gaussian Chain Step 4.
        for(int i = 0; i <= M_x/2; i++) {
            h_j[i][0] = a_j[i][0] * exp(-2.0 * pow(M_PI,2.0) * pow(b,2.0) * delta_s * pow((double)i,2.0) / (3.0 * L * L)) / M_x;
            h_j[i][1] = a_j[i][1] * exp(-2.0 * pow(M_PI,2.0) * pow(b,2.0) * delta_s * pow((double)i,2.0) / (3.0 * L * L)) / M_x;
        }
        for(int i = M_x / 2+1; i < M_x; i++) {
            double reverse_index = i - M_x;//NOT AN ARRAY INDEX.
            h_j[i][0] = a_j[i][0] * exp(-2.0 * pow(M_PI,2.0) * pow(b,2.0) * delta_s * pow(reverse_index,2.0) / (3.0 * L * L)) / M_x;
            h_j[i][1] = a_j[i][1] * exp(-2.0 * pow(M_PI,2.0) * pow(b,2.0) * delta_s * pow(reverse_index,2.0) / (3.0 * L * L)) / M_x;
        }

        //Calculate q^(n+2/3) according to Pseudo Spectral Algorithm 3.1: Continuous Gaussian Chain Step 5.
        fftw_execute(p2);

        //Rescale since inverse is scaled by M_x
        //According to FFTW documentation: http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html
        for(int i = 0; i < M_x; i++) {
            //q_complex[i][0] /= M_x;
            //q_complex[i][1] /= M_x;
            q[i][n+1] = q_complex[i][0];
        }
        //Inversed FFT to obtain q^(n+2/3)

        //Apply Q operator to get q^(n+1) according to Pseudo Spectral Algorithm 3.1: Continuous Gaussian Chain Step 6.
        apply_q_operator(q, n+1,n+1);

        //Increment s.
        s += delta_s;

    }
    //Cleaning up.
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(h_j);
    fftw_free(a_j);
    fftw_free(q_complex);
    fftw_cleanup();

}

int main(int argc, char** argv) {
    std::cout << "R_g is: " << R_g << std::endl;
    double q[M_x][N_s+1];
    //q is a 2D array to keep track of the value of q at various s points as well as to keep
    // track of q at the M_x collocation points (points of interest)
    double Q = 0.0;
    another_way(q,Q);
    std::cout << "Q: " << Q << std::endl;
    std::ofstream output("output.csv");



    output << "x,q(x),w(x),density" << '\n';
    for(int i = 0; i < M_x; i++) {
        double x = delta_x * (double)i ;
        //double x = L / (M_x+1) * i;
        double actual_value = q[i][N_s-1];


        // Simpsons 1/3 rule to integrate the density using Equation 3.243
        //Source: https://en.wikipedia.org/wiki/Simpson%27s_rule

        double sum = 0;
        double h = delta_s;
        for(int s = 0; s < N_s+1; s++) {
            if(s == 0) {
                sum = sum + h / 3 * q[i][s] * q[i][N_s - s];
            } else if(s == N_s) {
                sum = sum + h / 3 * q[i][s] * q[i][N_s - s];
            } else if(s % 2 == 1) {
                sum = sum + h / 3 * 4 * q[i][s] * q[i][N_s - s];
            } else if(s % 2 == 0) {
                sum = sum + h / 3 * 2 * q[i][s] * q[i][N_s - s];
            }
        }
        //Obtain the zero frequency component as Q(w) = a_0(N), according to FFTW
        //the zero frequency is at index 0.
        sum = sum / (Q * N);
        output << x << "," << actual_value << "," << w(x) <<  "," << sum << '\n';
    }
    output.close();
    return 0;
}