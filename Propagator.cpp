//
// Created by callum on 20/01/2020.
//

#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include "mathutils.h"
#include "Propagator.h"

Propagator::Propagator(int M_x, int N_s, double L, double N, double R_g, double f, double chiN) : M_x(M_x), N_s(N_s), R_g(R_g), L(L), f(f), N(N), chiN(chiN) {
    //M_x collaction points.
    //N_s integration steps.
    //L Box Length (R_g).
    //N Chain Length
    delta_s = N / N_s;
    delta_x = L / M_x;
    b = sqrt(6) * R_g / sqrt(N);

    q = new double*[M_x];
    for(int i = 0; i < M_x; i++) {
        q[i] = new double[N_s+1];
    }
    for(int m = 0; m < M_x; m++) {
        q[m][0] = 1;
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
    std::cout << "Box Length L: " << L << '\n';
    std::cout << "Chain Length N: " << N << '\n';
    std::cout << "Radius of gyration R_g: " << R_g << '\n';
    std::cout << "Delta_x: " << delta_x << '\n';
    std::cout << "Delta s: " << delta_s << '\n';
    std::cout << "b: " << b << '\n';
}



void Propagator::Propagate() {
    //Step along ds.
    s = 0;
    for(int i = 0; i < N_s; i++) {
        Propagate_n(i);
    }
}




double Propagator::w(double x, int n) {
    return (1.0  - 2.0 * pow(sech(3.0 * (x - L / 2.0)/ (2 * R_g)),2.0)) / N;

}

void Propagator::Apply_W_To_Q(double **q, int src_N_s, int dest_N_s) {
    for(int j = 0; j < M_x; j++) {
        double x = delta_x * (double)j;
        q[j][dest_N_s] = q[j][src_N_s] * exp(-1.0 * w(x, src_N_s)   * delta_s  / 2);
        //Fig: 3.14, divided by N, as the potential function is scaled Nw(x)=function in w.
        //This leads to the unscaled potential function is therefore w(x)/N.
    }
}

double Propagator::Find_Max() {
    double max = q[0][N_s];
    for(int i = 0; i < M_x; i++) {
        if(q[i][N_s] > max) {
            max = q[i][N_s];
        }
    }
    return max;
}

double Propagator::Find_Min() {
    double min = q[0][N_s];
    for(int i = 0; i < M_x; i++) {
        if(q[i][N_s] < min) {
            min = q[i][N_s];
        }
    }
    return min;
}


double Propagator::DetermineQ() {
    double sum = 0.0;
    for(int i =0 ; i < M_x; i++) {
        sum = sum + q[i][N_s-1] * delta_x;
    }
    return 1.0 / M_x * sum;
}


void Propagator::Propagate_n(int n) {
    //Initialise M column vector.
    //Calculate next step and get q^(n+1/3) according to Pseudo Spectral Algorithm 3.1: Continuous Gaussian Chain Step 2.
    Apply_W_To_Q(q, n,n+1);//Populate current Q from previous Q.

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
    Apply_W_To_Q(q, n+1,n+1);
    s = s + delta_s;

}

Propagator::~Propagator() {
    //Cleaning up.
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(h_j);
    fftw_free(a_j);
    fftw_free(q_complex);
    fftw_cleanup();
}

void Propagator::Save(const char *file_name) {
    remove(file_name);
    std::ofstream output;
    output.open(file_name);
    output << "x,q,w,density" << '\n';
    double phi_A[M_x];
    double phi_B[M_x];
    for(int i = 0; i < M_x; i++) {
        double x = delta_x * (double)i ;
        // Simpsons 1/3 rule to integrate the density using Equation 3.243
        //Source: https://en.wikipedia.org/wiki/Simpson%27s_rule
        double sum = 0;
        double h = delta_s;
        for(int s = 0; s < f * (N_s+1); s++) {
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
        phi_A[i] = sum;
    }
    for(int i = 0; i < M_x; i++) {
        double x = delta_x * (double)i ;
        // Simpsons 1/3 rule to integrate the density using Equation 3.243
        //Source: https://en.wikipedia.org/wiki/Simpson%27s_rule
        double sum = 0;
        double h = delta_s;
        for(int s = f*(N_s+1); s < (N_s+1); s++) {
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
        phi_B[i] = sum;
    }
    for(int i = 0; i < M_x; i++) {
        double x = delta_x * (double)i ;
        double density = phi_A[i] + phi_B[i];
        output << x << "," << q[i][N_s - 1] << "," << w(x,N_s-1) << "," << density << '\n';
    }
    output.close();
}

void Propagator::Set_Fields(double *w_A, double *w_B) {
    this->w_A = w_A;
    this->w_B = w_B;
}

double** Propagator::GetPropagator() {
    return q;
}

