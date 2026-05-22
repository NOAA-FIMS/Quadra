#include "../core/autodiff.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>

DECLARE_ADGRAPH()

struct QuadraIFTModel
{
    template <typename T>
    T operator()(const T& theta,
                 const T& u1,
                 const T& u2) const
    {
        return
            T(0.5) * (u1 - theta) * (u1 - theta)
            + T(0.5) * (u2 + T(0.5) * theta) * (u2 + T(0.5) * theta)
            + T(0.25) * u1 * u1
            + T(0.15) * u2 * u2
            + T(0.10) * u1 * u2;
    }
};

Eigen::Vector2d solve_mode(double theta)
{
    Eigen::Matrix2d H;
    H << 1.5, 0.10,
         0.10, 1.30;

    Eigen::Vector2d b;
    b << theta,
        -0.5 * theta;

    return H.ldlt().solve(b);
}

Eigen::Vector2d finite_difference_du_dtheta(
    double theta,
    double h = 1.0e-6)
{
    return
        (solve_mode(theta + h) - solve_mode(theta - h)) /
        (2.0 * h);
}

int main()
{
    const double theta = 2.0;
    const Eigen::Vector2d uhat = solve_mode(theta);
    const Eigen::Vector2d fd = finite_difference_du_dtheta(theta);

    had::ADGraph graph;
    had::g_ADGraph = &graph;

    quadra::AD theta_ad(theta);
    quadra::AD u1_ad(uhat[0]);
    quadra::AD u2_ad(uhat[1]);

    QuadraIFTModel model;

    const quadra::AD nll =
        model(theta_ad, u1_ad, u2_ad);

    had::ZeroAdjoints(graph);
    had::SetAdjoint(nll, 1.0);
    had::PropagateAdjointDirectional();

    Eigen::Matrix2d Huu;
    Huu << had::GetAdjoint(u1_ad, u1_ad), had::GetAdjoint(u1_ad, u2_ad),
           had::GetAdjoint(u2_ad, u1_ad), had::GetAdjoint(u2_ad, u2_ad);

    Eigen::Vector2d HuTheta;
    HuTheta << had::GetAdjoint(u1_ad, theta_ad),
               had::GetAdjoint(u2_ad, theta_ad);

    const Eigen::Vector2d ift =
        -Huu.ldlt().solve(HuTheta);

    const double error =
        (ift - fd).norm();

    if (!std::isfinite(error)) {
        std::cerr << "FAIL: non-finite Quadra IFT block error\n";
        return 1;
    }

    if (error > 1.0e-8) {
        std::cerr << "FAIL: Quadra IFT block mismatch\n";
        std::cerr << "Huu:\n" << Huu << "\n";
        std::cerr << "HuTheta:\n" << HuTheta << "\n";
        std::cerr << "IFT:\n" << ift << "\n";
        std::cerr << "FD:\n" << fd << "\n";
        std::cerr << "error: " << error << "\n";
        return 1;
    }

    std::cout << "PASS: Quadra AD IFT block validation\n";
    std::cout << "Huu:\n" << Huu << "\n";
    std::cout << "HuTheta:\n" << HuTheta << "\n";
    std::cout << "IFT du/dtheta:\n" << ift << "\n";
    std::cout << "FD du/dtheta:\n" << fd << "\n";
    std::cout << "error: " << error << "\n";

    return 0;
}
