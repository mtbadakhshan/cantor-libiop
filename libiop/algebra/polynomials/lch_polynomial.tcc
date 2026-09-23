#include <cstddef>
#include <vector>


namespace libiop
{
    template<typename FieldT>
    FieldT lch_polynomial_evaluation_at_point(const std::vector<FieldT> &lch_coeffs,
                                            const FieldT &x)
    {
        if (lch_coeffs.empty()) {
            return FieldT::zero();
        }

        size_t n = 1;
        while (n < lch_coeffs.size()) {
            n <<= 1;
        }
        const size_t d = libff::log2(n);

        std::vector<FieldT> X(d);
        if (d > 0) {
            X[0] = x;
            for (size_t r = 1; r < d; ++r) {
                X[r] = X[r - 1].squared() + X[r - 1];
            }
        }

        std::vector<FieldT> c(n, FieldT::zero());
        std::copy(lch_coeffs.begin(), lch_coeffs.end(), c.begin());

        for (size_t r = d; r-- > 0;) {
            const size_t half = size_t(1) << r;
            for (size_t j = 0; j < half; ++j) {
                c[j] += X[r] * c[j + half];
            }
        }
        return c[0];
    }

} // namespace libiop