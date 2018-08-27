#pragma once
#include <complex>

class FFT
{
	unsigned int N, which;
	unsigned int log_2_N;
	float pi2;
	unsigned int *reversed;
	std::complex<float> **W;
	std::complex<float> *c[2];

public:
	FFT(unsigned int N);
	~FFT();
	unsigned int reverse(unsigned int i) const;
	std::complex<float> w(unsigned int x, unsigned int N) const;
	void fft(std::complex<float>* input, std::complex<float>* output, int stride, int offset);
};