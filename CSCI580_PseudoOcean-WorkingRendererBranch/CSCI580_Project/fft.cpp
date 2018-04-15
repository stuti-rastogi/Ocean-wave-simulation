#include "fft.h"

#define _USE_MATH_DEFINES
#include <math.h>
#undef _USE_MATH_DEFINES

FFT::FFT(unsigned int N)
	: N(N), pi2(2 * M_PI), reversed(nullptr), W(nullptr)
{
	c[0] = c[1] = nullptr;
	
	log_2_N = log(N) / log(2);

	reversed = new unsigned int[N];
	for (int i = 0; i < N; ++i)
		reversed[i] = reverse(i);

	int pow2 = 1;
	W = new std::complex<float>*[log_2_N];
	for (int i = 0; i < log_2_N; ++i)
	{
		W[i] = new std::complex<float>[pow2];
		for (int j = 0; j < pow2; ++j)
			W[i][j] = w(j, pow2 * 2);
		pow2 *= 2;
	}

	c[0] = new std::complex<float>[N];
	c[1] = new std::complex<float>[N];
	which = 0;
}

FFT::~FFT()
{
	if (c[0]) delete[] c[0];
	if (c[1]) delete[] c[1];
	if (W)
	{
		for (int i = 0; i < log_2_N; ++i)
			if (W[i])
				delete[] W[i];
		delete[] W;
	}
	if (reversed)
		delete[] reversed;
}

unsigned int FFT::reverse(unsigned int i) const
{
	unsigned int res = 0;
	for (int j = 0; j < log_2_N; j++) {
		res = (res << 1) + (i & 1);
		i >>= 1;
	}
	return res;
}

std::complex<float> FFT::w(unsigned x, unsigned N) const
{
	return std::complex<float>(cos(pi2 * x / N), sin(pi2 * x / N));
}

void FFT::fft(std::complex<float>* input, std::complex<float>* output, int stride, int offset)
{
	for (int i = 0; i < N; i++) c[which][i] = input[reversed[i] * stride + offset];

	int loops = N >> 1;
	int size = 1 << 1;
	int size_over_2 = 1;
	int w_ = 0;
	for (int i = 1; i <= log_2_N; i++) {
		which ^= 1;
		for (int j = 0; j < loops; j++) {
			for (int k = 0; k < size_over_2; k++) {
				c[which][size * j + k] = c[which ^ 1][size * j + k] +
					c[which ^ 1][size * j + size_over_2 + k] * W[w_][k];
			}

			for (int k = size_over_2; k < size; k++) {
				c[which][size * j + k] = c[which ^ 1][size * j - size_over_2 + k] -
					c[which ^ 1][size * j + k] * W[w_][k - size_over_2];
			}
		}
		loops >>= 1;
		size <<= 1;
		size_over_2 <<= 1;
		w_++;
	}

	for (int i = 0; i < N; i++) output[i * stride + offset] = c[which][i];
}

