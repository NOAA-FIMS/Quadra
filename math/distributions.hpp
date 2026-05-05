#ifndef DISTRIBUTIONS_HPP
#define DISTRIBUTIONS_HPP
namespace quadra
{
    template <typename T>
    T dnorm(T x, T mean, T sd, bool log = true)
    {
        T val = -0.5 * pow((x - mean) / sd, 2) - log(sd);
        return log ? val : exp(val);
    }
} // namespace pelagia

#endif // DISTRIBUTIONS_HPP