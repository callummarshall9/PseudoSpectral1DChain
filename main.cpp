#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <fftw3.h>

//Simulates single chain of homopolymer.

const double L = 10.0; // box length
const double N = 10.0; // Chain length
const int N_s = 201; // Number of integration steps along the chain
const double delta_s = N/(N_s-1.0); // chain integration step size
const double delta_x = L / (N+1);
const int M_x = 64; // Grid points in x direction
const double q_0 = 1.0;
const double R_g = 1.0; // Chain gyration radius (unperturbed radius of gyration)
const double b = 1.0;


double sech(double x) {
    return 1.0 / cosh(x);
}

double w(double x) {
    return 1.0  - 2.0 * pow(sech(3.0 * (x - L / 2.0)/ (2 * R_g)),2.0);
}

/*
 * FFTW code, though seems to have same issue as code above. Peak height is largely affected by b
 */

void apply_q_operator(double q[M_x]) {
    for(int j = 0; j < M_x; j++) {
        //double x = (double)L / (double)M_x * (double)j;
        double x = L / (M_x+1) * j;
        //Source: https://en.wikipedia.org/wiki/Pseudo-spectral_method
        q[j] = q[j] * exp(-1.0 * w(x) / N * delta_s  / 2);
        //Fig: 3.14, divided by N, as the potential function is scaled Nw(x)=function in w..
        //This leads to the unscaled potential function is therefore w(x)/N.
    }
}

void another_way(double q[M_x]) {
    for(int m = 0; m < M_x; m++) {
        q[m] = q_0;
    }

    for(int n = 0; n < N_s; n++) {
        fftw_plan p1,p2;
        fftw_complex* a_j = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * M_x);
        fftw_complex* h_j = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * M_x);
        fftw_complex* q_complex = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * M_x);
        //Initialise M column vector.
        //Calculate next step and get q^(n+1/3)
        apply_q_operator(q);
        for(int i = 0; i < M_x; i++) {
            q_complex[i][0] = q[i];
            q_complex[i][1] = 0;
        }
        p1 = fftw_plan_dft_1d(M_x, q_complex, a_j, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
        fftw_execute(p1);
        fftw_destroy_plan(p1);
        for(int i = 0; i <= M_x/2; i++) {
            h_j[i][0] = a_j[i][0] * exp(-2.0 * pow(M_PI,2.0) * pow(b,2.0) * delta_s * pow((double)i,2.0) / (3.0 * L * L));
            h_j[i][1] = a_j[i][1] * exp(-2.0 * pow(M_PI,2.0) * pow(b,2.0) * delta_s * pow((double)i,2.0) / (3.0 * L * L));
        }
        for(int i = M_x / 2 + 1; i < M_x; i++) {
            double reverse_index = M_x - i;
            h_j[i][0] = a_j[i][0] * exp(-2.0 * pow(M_PI,2.0) * pow(b,2.0) * delta_s * pow(reverse_index,2.0) / (3.0 * L * L));
            h_j[i][1] = a_j[i][1] * exp(-2.0 * pow(M_PI,2.0) * pow(b,2.0) * delta_s * pow(reverse_index,2.0) / (3.0 * L * L));
        }
        p2 = fftw_plan_dft_1d(M_x, h_j, q_complex, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
        fftw_execute(p2);
        fftw_destroy_plan(p2);
        //Rescale since inverse is scaled by M_x
        for(int i = 0; i < M_x; i++) {
            q_complex[i][0] /= M_x;
            q_complex[i][1] /= M_x;
            q[i] = q_complex[i][0];
        }
        apply_q_operator(q);
        fftw_free(h_j);
        fftw_free(a_j);
        fftw_free(q_complex);
        fftw_cleanup();
    }


}

int main(int argc, char** argv) {
    double q[M_x];
    another_way(q);
    std::ofstream output("output.csv");
    for(int i = 0; i < M_x; i++) {
        //double x = (double)L / (double)M_x * (double)i;
        double x = L / (M_x+1) * i;
        //Due to FFTReal scaling. And the fact it's been scaled N_s times as well as the fact that due to complex coefficients, the output is actually (M_x/2) for what it is divided by when
        //scaling. This results in a 1/(N_s * M_x / 2) scaling, so to fix this:
        double actual_value = q[i];

        output << x << "," << actual_value << "," << w(x) <<  '\n';
    }
    output.close();

    return 0;
}