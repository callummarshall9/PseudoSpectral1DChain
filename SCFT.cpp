//
// Created by callum on 20/01/2020.
//

#include <fstream>
#include "SCFT.h"

SCFT::SCFT(int M_x, int N_s, double L, double N, double R_g, double f, double flory_higgs, double mixing_parameter) : gen(rd()), dist(0.0,1.0),
    M_x(M_x), N_s(N_s), L(L), N(N), R_g(R_g), f(f), flory_higgs(flory_higgs), mixing_parameter(mixing_parameter), q_propagator(M_x, N_s, L, N, R_g, f, flory_higgs),
                                                                                                                      q_star_propagator(M_x, N_s, L, N, R_g, f, flory_higgs) {
    w_A = new double[M_x];
    w_B = new double[M_x];
    for(int i = 0; i < M_x; i++) {
        w_A[i] = dist(gen);
        w_B[i] = dist(gen);
    }
    q_propagator.Set_Fields(w_A, w_B);
    q_star_propagator.Set_Fields(w_A, w_B);
    phi_A = new double[M_x];
    phi_B = new double[M_x];
    delta_s = N / N_s;
    delta_x = L / M_x;
}

double SCFT::Determine_Error() {
    double sum = 0.0;
    for(int i = 0; i < M_x; i++) {
        double w_A_change = phi_B[i] + 0.5 * (w_A[i] + w_B[i]);
        double w_B_change = phi_A[i] + 0.5 * (w_A[i] + w_B[i]);
        sum = sum + (pow(abs(w_A_change - w_A[i]), 2.0) + pow(abs(w_B_change - w_B[i]), 2.0)) * delta_x;
    }
    sum = sqrt(1.0 / M_x * sum) * flory_higgs;
    return sum;
}

void SCFT::Update_Fields() {
    for(int i =0; i < M_x; i++) {
        double w_A_change = flory_higgs * phi_B[i] + 0.5 * (w_A[i] + w_B[i]);
        w_A[i] = (w_A[i] * (1.0 - mixing_parameter) + mixing_parameter * w_A_change);
        double w_B_change = flory_higgs * phi_A[i] + 0.5 * (w_A[i] + w_B[i]);
        w_B[i] = (w_B[i] * (1.0 - mixing_parameter) + mixing_parameter * w_B_change);
    }
}

void SCFT::Run() {
    q_propagator.Output_Parameters();
    q_propagator.Propagate();
    q_star_propagator.Propagate();
    Determine_Density_Differences();
    double error = Determine_Error();
    int index = 1.0;
    std::string file_name = "hell_" + std::to_string(index) + ".csv";
    Save(file_name);
    std::cout << "Error: " << error << std::endl;
    index++;
    while(error > pow(10.0, -6.0)) {
        std::cout << "Error: " << error << std::endl;
        q_propagator.Propagate();
        q_star_propagator.Propagate();
        Determine_Density_Differences();
        Update_Fields();
        std::string file_name = "hell_" + std::to_string(index) + ".csv";
        //Save(file_name);
        error = Determine_Error();
        index++;
    }

}

void SCFT::Save(std::string file_name) {
    std::ofstream output;
    output.open(file_name);
    output << "x,q,q*,w_A,w_B,phi_A,phi_B,density" << '\n';
    for(int i = 0; i < M_x; i++) {
        double** q = q_propagator.GetPropagator();
        double** q_star = q_propagator.GetPropagator();
        double x = delta_x * i;
        double density = phi_A[i] + phi_B[i];
        output << x << "," << q[i][N_s - 1] << "," << q_star[i][N_s-1] <<","<< w_A[i] << "," << w_B[i] << "," << phi_A[i] << "," << phi_B[i] << "," << density << '\n';
    }
    output.close();
}

void SCFT::Determine_Density_Differences() {
    double** q = q_propagator.GetPropagator();
    double** q_star = q_star_propagator.GetPropagator();
    double q_Q = q_propagator.GetQ();
    double q_star_Q = q_star_propagator.GetQ();
    if(q_Q != q_star_Q) {
        std::cout << "Warning q_Q != q_star_Q" << '\n';
        std::cout << q_Q << "," << q_star_Q << '\n';
    }
    double Q = q_Q;
    for(int i = 0; i < M_x; i++) {
        double x = delta_x * (double)i ;
        // Simpsons 1/3 rule to integrate the density using Equation 3.243
        //Source: https://en.wikipedia.org/wiki/Simpson%27s_rule
        double sum = 0;
        double h = delta_s;
        for(int s = 0; s < f * (N_s+1); s++) {
            if(s == 0) {
                sum = sum + h / 3 * q[i][s] * q_star[i][s];
            } else if(s == N_s) {
                sum = sum + h / 3 * q[i][s] * q_star[i][s];
            } else if(s % 2 == 1) {
                sum = sum + h / 3 * 4 * q[i][s] * q_star[i][s];
            } else if(s % 2 == 0) {
                sum = sum + h / 3 * 2 * q[i][s] * q_star[i][s];
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
        for(int s = f * (N_s + 1); s < N_s+1; s++) {
            if(s == 0) {
                sum = sum + h / 3 * q[i][s] * q_star[i][s];
            } else if(s == N_s) {
                sum = sum + h / 3 * q[i][s] * q_star[i][s];
            } else if(s % 2 == 1) {
                sum = sum + h / 3 * 4 * q[i][s] * q_star[i][s];
            } else if(s % 2 == 0) {
                sum = sum + h / 3 * 2 * q[i][s] * q_star[i][s];
            }
        }
        //Obtain the zero frequency component as Q(w) = a_0(N), according to FFTW
        //the zero frequency is at index 0.
        sum = sum / (Q * N);
        phi_B[i] = sum;
    }
}