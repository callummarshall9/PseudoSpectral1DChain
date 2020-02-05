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
        w_A[i] = dist(gen) * 0.1;
        w_B[i] = dist(gen) * 0.1;
        //w_A[i] = 1;
        //w_B[i] = 0;
    }

    q_propagator.Set_Fields(w_A, w_B);
    q_star_propagator.Set_Fields(w_A, w_B);
    phi_A = new double[M_x];
    phi_B = new double[M_x];
    delta_s = N / N_s;
    delta_x = L / M_x;
}

double SCFT::Determine_Error(int index) {
    bool output = false;
    //std::ofstream output_fields("output_fields_" + std::to_string(index) + ".csv");

    double sum = 0.0;
    for(int i = 0; i < M_x; i++) {
        double compressibility_condiiton = 0.5 * (w_A[i] + w_B[i] - flory_huggins * N);
        double w_a_out = flory_huggins * N * phi_B[i] + compressibility_condiiton;
        double w_b_out = flory_huggins * N * phi_A[i] + compressibility_condiiton;
        //output_fields << i << "," << w_A[i] << "," <<  w_a_out << "," << w_B[i] << "," << w_b_out << '\n';
        double difference_in_w_A = fabs(pow(w_a_out - w_A[i], 2.0));
        double difference_in_w_B = fabs(pow(w_b_out - w_B[i], 2.0));
        sum = sum + (difference_in_w_A + difference_in_w_B) * delta_x;
    }
    sum = sqrt(1.0 / (M_x * delta_x) * sum) / (flory_huggins * N);
    //output_fields.close();
    return sum;
}

double SCFT::DetermineVariance() {
    double mean_squares = 0.0;
    for(int i = 0; i < M_x; i++) {
        mean_squares += pow(phi_A[i],2.0);
    }
    mean_squares /= M_x;
    double mean = 0.0;
    for(int i = 0; i < M_x; i++) {
        mean += phi_A[i];
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
        double gamma = 0.1;
        double C = 20000;
        //Perhaps source of error is by factor of N.
        for (int i = 0; i < M_x; i++) {
            //Note: w_B = omega_- and w_A = omega_+
            w_B[i] = w_B[i] - delta_t * gamma * delta_x * (-phi_B[i] + 2.0 / (flory_huggins * N) * w_B[i]) +
                     gauss_noise(norm_generator);
        }
        q_propagator.Propagate();
        q_star_propagator.Propagate();
        Determine_Density_Differences();
        while (DetermineVariance() < pow(10.0, -5)) {
            for (int i = 0; i < M_x; i++) {
                w_A[i] = w_A[i] - delta_t * gamma * C * (phi_A[i] - 1.0);
            }
            q_propagator.Propagate();
            q_star_propagator.Propagate();
            Determine_Density_Differences();
        }
    } else if(field_method == SimpleMixing) {
        for(int i = 0; i < M_x; i++) {
            double compressibility_condition = 0.5 * (w_A[i] + w_B[i] - flory_huggins * N);
            double w_a_out = flory_huggins * N *  phi_B[i] + compressibility_condition;
            double w_b_out = flory_huggins * N * phi_A[i] + compressibility_condition;
            w_A[i] = (w_A[i] * (1.0 - mixing_parameter) + mixing_parameter * w_a_out);
            w_B[i] = (w_B[i] * (1.0 - mixing_parameter) + mixing_parameter * w_b_out);
        }
    }

}

void SCFT::Run() {
    int index = 1.0;
    double field_error_threshold = pow(10.0,-3.0);
    double variance_threshold = pow(10.0,-5.0);
    double field_error = Determine_Error(index);
    double variance_error = sqrt(DetermineVarianceTotal());
    while(field_error > field_error_threshold || variance_error > variance_threshold) {
        std::cout << index << " - Error: " << Determine_Error(index) << '\n';
        std::cout << index << " - Stdev: " << sqrt(DetermineVarianceTotal()) << '\n';
        q_propagator.Propagate();
        q_star_propagator.Propagate();
        Determine_Density_Differences();
        Update_Fields();
        std::string file_name = "hell_" + std::to_string(index) + ".csv";
        //Save(file_name);
        if(index % 200 == 0) {
            Save(file_name);
        }
        index++;
        if(index > 2000) {
            break;
        }
        field_error = Determine_Error(index);
        variance_error = sqrt(DetermineVarianceTotal());
    }
    std::cout << index << " - Error: " << field_error << '\n';
    std::cout << index << " - Stdev: " << variance_error << '\n';
    std::string file_name = "hell_" + std::to_string(index) + ".csv";
    Save(file_name);
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
    double q_Q = q_propagator.Q;
    double q_star_Q = q_star_propagator.Q;
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

}