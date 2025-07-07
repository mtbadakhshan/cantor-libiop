<h1 align="center">cantor-libiop: A Fork of libiop with Accelerated Additive FFT</h1>

**🚀 Accelerated `libiop` by Optimizing Additive FFT Algorithms**

This repository is a **fork of [`libiop`](https://github.com/scipr-lab/libiop)** focusing on **optimizing additive Fast Fourier Transform (FFT) algorithms** used in the original repository.

This library uses the [Cantor](https://www.sciencedirect.com/science/article/pii/0097316589900204?via%3Dihub) additive FFT implemented in the [`additive-fft`](https://github.com/mtbadakhshan/additive-fft/) library and demonstrates notable performance gains in the current __Aurora__ [[BCRSVW]](https://eprint.iacr.org/2018/828) implementation over the [Gao–Mateer](https://ieeexplore.ieee.org/document/5625613) (GM) algorithm provided in [`libiop`](https://github.com/scipr-lab/libiop) across all input sizes, and over the [LCH](https://ieeexplore.ieee.org/document/6979016) algorithm re-implemented in [`additive-fft`](https://github.com/mtbadakhshan/additive-fft/) for smaller circuits, which are prevalent in many zkSNARK applications.


We achieved about **~40% performance improvement** on the  __Aurora__  while maintaining correctness and interface compatibility with upstream `libiop`.


<span style="color:red">**WARNING:**</span> This is an academic proof-of-concept prototype, and in particular has not received careful code review. <br> This implementation is NOT ready for production use.



## Installation

Please follow the [installation instructions](INSTALL.md).


## Benchmarking

| $\log⁡_2(N)$ | $\log_2(∣L∣)$ | Prover: GM | Prover: Cantor | Prover: LCH | Verifier: GM | Verifier: Cantor | Verifier: LCH |
|:-------------:|:---------------:|:--------------:|:------------------:|:-------------------:|:----------------:|:--------------------:|:---------------------:|
| 9 | 16 | 0.44 | 0.33 | 0.35 | 0.04 | 0.04 | 0.04 |
| 10 | 17 | 0.90 | 0.67 | 0.71 | 0.05 | 0.05 | 0.05 |
| 11 | 18 | 1.87 | 1.36 | 1.44 | 0.07 | 0.06 | 0.07 |
| 12 | 19 | 3.99 | 2.91 | 2.93 | 0.10 | 0.09 | 0.10 |
| 13 | 20 | 8.53 | 6.02 | 5.95 | 0.17 | 0.15 | 0.16 |
| 14 | 21 | 19.47 | 12.01 | 12.44 | 0.29 | 0.26 | 0.28 |
| 15 | 22 | 41.05 | 25.27 | 25.41 | 0.54 | 0.48 | 0.52 |
| 16 | 23 | 84.26 | 50.83 | 50.63 | 1.02 | 0.93 | 1.00 |
| 17 | 24 | 176.67 | 104.26 | 102.95 | 1.98 | 1.79 | 1.93 |
| 18 | 25 | 373.83 | 216.00 | 213.61 | 3.88 | 3.51 | 3.78 |
| 19 | 26 | 771.42 | 443.88 | 441.51 | 7.78 | 6.91 | 7.44 |

All measurements were taken with [Google Benchmark](https://github.com/google/benchmark), with a minimum **10-second warm-up period**, on:

- **CPU:** AMD Ryzen 9 9950X @ 5.7 GHz  
- **RAM:** 64 GB DDR5  
- **OS:** Debian 12 with kernel 6.12.12

The codeword length is set to $|L| = 2^{7} N$, adopting [Preon](https://csrc.nist.gov/csrc/media/Projects/pqc-dig-sig/documents/round-1/spec-files/Preon-spec-web.pdf)’s choice of $\rho = 2^{-5}$ to tighten the FRI soundness error. Two extra dimensions are added to achieve **zero-knowledge** property.

## License

This library is licensed under the [MIT License](LICENSE).

## Acknowledgements

This work was supported by [BTQ Technologies Corp.](https://www.btq.com/) and [Mitacs](https://www.mitacs.ca/).

<!-- ## 📄 How to Cite
If you use this repository or the additive FFT optimizations in your research, please cite our paper:

```
@InProceedings{badakhshan2025accelerating,
author = {Mohammadtaghi Badakhshan and Susanta Samanta and Guang Gong},
title = {Accelerating Post-quantum Secure zkSNARKs by Optimizing Additive FFT},
booktitle = {Selected Areas in Cryptography},
year = {2025}
}
``` -->
