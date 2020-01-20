//
// Created by callum on 20/01/2020.
//

#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include "Propagator.h"

Propagator::Propagator(int M_x, int N_s, double L, double N, double R_g, double f, double flory_higgs) : M_x(M_x), N_s(N_s), R_g(R_g), L(L), f(f), N(N), flory_higgs(flory_higgs) {
    //M_x collaction points.
    //N_s integration steps.
    //L Box Length (R_g).
    //N Chain Length
    delta_s = N / N_s;
    delta_x = L / M_x;
    b = sqrt(6) * R_g / sqrt(N);

    q = new double*[M_x];
    for(int i = 0; i < M_x; i++) {
        q[i] = new double[N_s];
        for(int j = 0; j < N_s; j++) {
            q[i][j] = 1.0;
        }
    }
    a_j = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * M_x);
    h_j = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * M_x);
    q_complex = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * M_x);
    p1 = fftw_plan_dft_1d(M_x, q_complex, a_j, FFTW_FORWARD, FFTW_ESTIMATE);
    p2 = fftw_plan_dft_1d(M_x, h_j, q_complex, FFTW_BACKWARD, FFTW_ESTIMATE);

}

void Propagator::Output_Parameters() {
    std::cout << "M_x: " << M_x << '\n';
    std::cout << "N_s: " << N_s << '\n';
    std::cout << "L: " << L << '\n';
    std::cout << "N: " << N << '\n';
    std::cout << "R_g: " << R_g << '\n';
    std::cout << "Delta_x: " << delta_x << '\n';
    std::cout << "Delta s: " << delta_s << '\n';
}

double Propagator::GetDeltaX() {
    return delta_x;
}

int Propagator::GetN_s() {
    return N_s;
}

double Propagator::GetF() {
    return f;
}



void Propagator::Propagate() {
    //Step along ds.
    for(int i = 0; i < N_s; i++) {
        Propagate_n(i);
    }
}

inline double sech(double x) {
    return 1.0 / cosh(x);
}


double Propagator::w(double x, int n) {
    std::cout << "Using the rubbish one." << '\n';
    return (1.0  - 2.0 * pow(sech(3.0 * (x - L / 2.0)/ (2 * R_g)),2.0)) / N;

}

void Propagator::Apply_Q_Operator(double **q, int src_N_s, int dest_N_s) {
    for(int j = 0; j < M_x; j++) {
        double x = delta_x * (double)j;
        q[j][dest_N_s] = q[j][src_N_s] * exp(-1.0 * w(x, src_N_s)   * delta_s  / 2);
        //Fig: 3.14, divided by N, as the potential function is scaled Nw(x)=function in w.
        //This leads to the unscaled potential function is therefore w(x)/N.
    }
}

double Propagator::DetermineQ() {
    double sum = 0.0;
    for(int i =0 ; i < M_x; i++) {
        sum = sum + q[i][N_s-1] * delta_x;
    }
    return 1.0 / M_x * sum;
}

double Propagator::GetN() {
    return N;
}

double Propagator::GetRg() {
    return R_g;
}

double Propagator::GetL() {
    return L;
}

void Propagator::Propagate_n(int n) {
    //Initialise M column vector.
    //Calculate next step and get q^(n+1/3) according to Pseudo Spectral Algorithm 3.1: Continuous Gaussian Chain Step 2.
    Apply_Q_Operator(q, n,n+1);//Populate current Q from previous Q.

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
    Apply_Q_Operator(q, n+1,n+1);
}

void Propagator::Save(const char *file_name) {
    remove(file_name);
    std::ofstream output;
    output.open(file_name);
    for(int i = 0; i < M_x; i++) {
        double x = delta_x * i;
        output << x << "," << q[i][N_s - 1] << "," << w(x,N_s-1) << '\n';
    }
    output.close();
}

double** Propagator::GetPropagator() {
    return q;
}

double Propagator::GetQ() {
    return Q;
    //return DetermineQ();
}