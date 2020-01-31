//
// Created by callum on 20/01/2020.
//

#include <fstream>
#include "FieldMethod.h"
#include "SCFT.h"

SCFT::SCFT(int M_x, int N_s, double L, double N, double R_g, double f, double flory_huggins, double mixing_parameter, FieldMethod field_method) : gen(rd()), dist(0.0,1.0),
    M_x(M_x), N_s(N_s), L(L), N(N), R_g(R_g), f(f), flory_huggins(flory_huggins), mixing_parameter(mixing_parameter), q_propagator(M_x, N_s, L, N, R_g, f, flory_huggins, field_method),
                                                                                                                      q_star_propagator(M_x, N_s, L, N, R_g, f, flory_huggins, field_method),
                                                                                                                      gauss_noise(0, sqrt(2 * 0.00001)),
                                                                                                                      field_method(field_method) {//TODO: Add sqrt(2.0.01) as parameter.
    w_A = new double[M_x];
    w_B = new double[M_x];

    for(int i = 0; i < M_x; i++) {
        //w_A[i] = dist(gen) * 0.1;
        //w_B[i] = dist(gen) * 0.1;
        w_A[i] = 1;
        w_B[i] = 0;
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
        double compressibility_condiiton = 0.5 * (w_A[i] + w_B[i] - flory_huggins);
        double w_a_out = flory_huggins * phi_B[i] + compressibility_condiiton;
        double w_b_out = flory_huggins * phi_A[i] + compressibility_condiiton;
        sum = sum + (pow(abs(w_a_out - w_A[i]), 2.0) + pow(abs(w_b_out - w_B[i]), 2.0)) * delta_x;
    }
    sum = sqrt(1.0 / M_x * sum) * flory_huggins;
    return sum;
}

double SCFT::DetermineVariance() {
    double mean_squares = 0.0;
    for(int i = 0; i < M_x; i++) {
        mean_squares += pow(phi_B[i],2.0);
    }
    mean_squares /= M_x;
    double mean = 0.0;
    for(int i = 0; i < M_x; i++) {
        mean += phi_B[i];
    }
    mean /= M_x;
    return mean_squares - mean * mean;
}

double SCFT::DetermineVarianceTotal() {
    double mean_squares = 0.0;
    for(int i = 0; i < M_x; i++) {
        mean_squares += pow(phi_A[i] + phi_B[i],2.0);
    }
    mean_squares /= M_x;
    double mean = 0.0;
    for(int i = 0; i < M_x; i++) {
        mean += phi_A[i] + phi_B[i];
    }
    mean /= M_x;
    return mean_squares - mean * mean;
}

void SCFT::Update_Fields() {
    if(field_method == Langevindgn) {
        double delta_t = 0.1;
        double gamma = 0.01;
        double C = 20000;
        for (int i = 0; i < M_x; i++) {
            //Note: w_B = omega_- and w_A = omega_+
            w_B[i] = w_B[i] - delta_t * gamma * delta_x * (-phi_B[i] + 2.0 / flory_huggins * w_B[i]) +
                     sqrt(delta_t) * gauss_noise(norm_generator);
        }
        q_propagator.Propagate();
        q_star_propagator.Propagate();
        Determine_Density_Differences();
        while (DetermineVariance() > pow(10.0, -5)) {
            for (int i = 0; i < M_x; i++) {
                w_A[i] = w_A[i] - delta_t * gamma * C * (phi_A[i] - 1.0);
            }
            q_propagator.Propagate();
            q_star_propagator.Propagate();
            Determine_Density_Differences();
        }
    } else if(field_method == SimpleMixing) {
        for(int i =0; i < M_x; i++) {
            double compressibility_condition = 0.5 * (w_A[i] + w_B[i] - flory_huggins);
            double w_a_out = flory_huggins * phi_B[i] + compressibility_condition;
            double w_b_out = flory_huggins * phi_A[i] + compressibility_condition;
            w_A[i] = (w_A[i] * (1.0 - mixing_parameter) + mixing_parameter * w_a_out);
            w_B[i] = (w_B[i] * (1.0 - mixing_parameter) + mixing_parameter * w_b_out);
        }
    }

}

void SCFT::Run() {
    int index = 1.0;
    while(Determine_Error() > pow(10.0,-6.0) || sqrt(DetermineVarianceTotal()) > pow(10.0, -6.0)) {
        std::cout << index << " - Error: " << Determine_Error() << '\n';
        std::cout << index << " - Boolean Error: " << (Determine_Error() > pow(10.0,-6.0)) << '\n';
        std::cout << index << " - Stdev: " << sqrt(DetermineVarianceTotal()) << '\n';
        std::cout << index << " - Boolean Stdev: " << (sqrt(DetermineVarianceTotal()) > pow(10.0, -10.0)) << '\n';
        q_propagator.Propagate();
        q_star_propagator.Propagate();
        Determine_Density_Differences();
        Update_Fields();
        std::string file_name = "hell_" + std::to_string(index) + ".csv";
        Save(file_name);
        index++;
    }

    std::cout << index << " - Error: " << Determine_Error() << '\n';
    std::cout << index << " - Boolean Error: " << (Determine_Error() > pow(10.0,-6.0)) << '\n';
    std::cout << index << " - Stdev: " << sqrt(DetermineVarianceTotal()) << '\n';
    std::cout << index << " - Boolean Stdev: " << (sqrt(DetermineVarianceTotal()) > pow(10.0, -10.0)) << '\n';
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

void SCFT::Cleanup() {
    q_propagator.Cleanup();
    q_star_propagator.Cleanup();
}

void SCFT::Determine_Density_Differences() {
    double** q = q_propagator.GetPropagator();
    double** q_star = q_star_propagator.GetPropagator();
    double q_Q = q_propagator.GetQ();
    double q_star_Q = q_star_propagator.GetQ();
    if(q_Q - q_star_Q > pow(10.0, -4.0)) {
        std::cout << "Warning q_Q != q_star_Q" << '\n';
        std::cout << q_Q << "," << q_star_Q << '\n';
    }
    double Q = q_Q;

    for(int i = 0; i < M_x; i++) {
        double sum = 0.0;
        for(int s = 0; s < f * N_s; s++) {
            sum += q[i][s] * q_star[i][N_s - s] * delta_s;
        }
        phi_A[i] = sum / (Q * N);
        sum = 0.0;
        for(int s = f * N_s; s < N_s; s++) {
            sum += q[i][s] * q_star[i][N_s - s] * delta_s;
        }
        phi_B[i] = sum / (Q * N);
    }
    /*
    for(int i = 0; i < M_x; i++) {
        double x = delta_x * (double)i ;
        // Simpsons 1/3 rule to integrate the density using Equation 3.243
        //Source: https://en.wikipedia.org/wiki/Simpson%27s_rule
        double sum = 0;
        double h = delta_s;
        for(int s = 0; s < f * (N_s+1); s++) {
            if(s == 0) {
                sum = sum + h / 3 * q[i][s] * q_star[i][N_s - s];
            } else if(s == N_s) {
                sum = sum + h / 3 * q[i][s] * q_star[i][N_s - s];
            } else if(s % 2 == 1) {
                sum = sum + h / 3 * 4 * q[i][s] * q_star[i][N_s - s];
            } else if(s % 2 == 0) {
                sum = sum + h / 3 * 2 * q[i][s] * q_star[i][N_s - s];
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
        for(int s = f*(N_s+1)+1; s < (N_s+1); s++) {
            if(s == 0) {
                sum = sum + h / 3 * q[i][s] * q_star[i][N_s - s];
            } else if(s == N_s) {
                sum = sum + h / 3 * q[i][s] * q_star[i][N_s - s];
            } else if(s % 2 == 1) {
                sum = sum + h / 3 * 4 * q[i][s] * q_star[i][N_s - s];
            } else if(s % 2 == 0) {
                sum = sum + h / 3 * 2 * q[i][s] * q_star[i][N_s - s];
            }
        }
        //Obtain the zero frequency component as Q(w) = a_0(N), according to FFTW
        //the zero frequency is at index 0.
        sum = sum / (Q * N);
        phi_B[i] = sum;
    }*/
}