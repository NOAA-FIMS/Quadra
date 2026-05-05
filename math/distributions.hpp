#ifndef DISTRIBUTIONS_HPP
#define DISTRIBUTIONS_HPP
namespace quadra
{
    template <typename T>
    T dnorm(const T &x, const T &mean, const T &sd, bool give_log = true)
    {
        const double log_sqrt_2pi = 0.91893853320467274178;

        T z = (x - mean) / sd;
        T log_density = -0.5 * z * z - log(sd) - log_sqrt_2pi;

        return give_log ? log_density : exp(log_density);
    }
} // namespace quadra

#endif // DISTRIBUTIONS_HPP