#include <Eigen/Dense>

#include <cmath>
#include <iostream>

struct VectorModel
{
    double operator()(
        double theta,
        const Eigen::Vector2d& u) const
    {
        const double u1 = u[0];
        const double u2 = u[1];

        return
            0.5 * (u1 - theta) * (u1 - theta)
            + 0.5 * (u2 + 0.5 * theta) * (u2 + 0.5 * theta)
            + 0.25 * u1 * u1
            + 0.15 * u2 * u2
            + 0.10 * u1 * u2;
    }
};

Eigen::Vector2d solve_mode(double theta)
{
    Eigen::Matrix2d H;

    H <<
        1.5, 0.10,
        0.10, 1.30;

    Eigen::Vector2d b;

    b <<
        theta,
        -0.5 * theta;

    return H.ldlt().solve(b);
}

Eigen::Vector2d finite_difference_du_dtheta(
    double theta,
    double h = 1.0e-6)
{
    const Eigen::Vector2d up =
        solve_mode(theta + h);

    const Eigen::Vector2d um =
        solve_mode(theta - h);

    return (up - um) / (2.0 * h);
}

int main()
{
    const double theta = 2.0;

    const Eigen::Vector2d uhat =
        solve_mode(theta);

    const Eigen::Vector2d fd =
        finite_difference_du_dtheta(theta);

    const double h = 1.0e-5;

    VectorModel model;

    Eigen::Matrix2d Huu;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            Eigen::Vector2d e_i =
                Eigen::Vector2d::Zero();

            Eigen::Vector2d e_j =
                Eigen::Vector2d::Zero();

            e_i[i] = h;
            e_j[j] = h;

            const double fpp =
                model(theta, uhat + e_i + e_j);

            const double fpm =
                model(theta, uhat + e_i - e_j);

            const double fmp =
                model(theta, uhat - e_i + e_j);

            const double fmm =
                model(theta, uhat - e_i - e_j);

            if (i == j)
            {
                Huu(i, j) =
                    (model(theta, uhat + e_i) -
                     2.0 * model(theta, uhat) +
                     model(theta, uhat - e_i)) /
                    (h * h);
            }
            else
            {
                Huu(i, j) =
                    (fpp - fpm - fmp + fmm) /
                    (4.0 * h * h);
            }
        }
    }

    Eigen::Vector2d HuTheta;

    for (int i = 0; i < 2; ++i)
    {
        Eigen::Vector2d e =
            Eigen::Vector2d::Zero();

        e[i] = h;

        const double fpp =
            model(theta + h, uhat + e);

        const double fpm =
            model(theta + h, uhat - e);

        const double fmp =
            model(theta - h, uhat + e);

        const double fmm =
            model(theta - h, uhat - e);

        HuTheta[i] =
            (fpp - fpm - fmp + fmm) /
            (4.0 * h * h);
    }

    const Eigen::Vector2d ift =
        -Huu.ldlt().solve(HuTheta);

    const double error =
        (ift - fd).norm();

    if (!std::isfinite(error)) {
        std::cerr << "FAIL: non-finite IFT vector error\n";
        return 1;
    }

    if (error > 1.0e-6) {
        std::cerr << "FAIL: vector IFT mismatch\n";
        std::cerr << "IFT:\n" << ift << "\n";
        std::cerr << "FD:\n" << fd << "\n";
        std::cerr << "error: " << error << "\n";
        return 1;
    }

    std::cout << "PASS: vector IFT validation\n";
    std::cout << "u_hat:\n" << uhat << "\n";
    std::cout << "IFT du/dtheta:\n" << ift << "\n";
    std::cout << "FD du/dtheta:\n" << fd << "\n";
    std::cout << "error: " << error << "\n";

    return 0;
}
