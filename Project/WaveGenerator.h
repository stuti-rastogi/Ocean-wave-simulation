#pragma once

// General Includes
#include <complex>
#include "fft.h"
#include <windows.h>
#include <d3d11.h>
#include <d3dx11.h>
#include <d3dcompiler.h>
#include <xnamath.h>
#include "resource.h"
#include <ctime>

struct WaveVertex
{
	XMFLOAT3 position, normal;
};

struct FFTData
{
	XMFLOAT3 htilde0, htilde0mk;
};

inline std::complex<float> complexScalarMultiply(const std::complex<float> c, const float scalar)
{
	return std::complex<float>(c.real() * scalar, c.imag() * scalar);
}

inline double dot(XMFLOAT3 A, XMFLOAT3 B)
{
	return (A.x * B.x) + (A.y * B.y) + (A.z * B.z);
}

inline double dot(XMFLOAT2 A, XMFLOAT2 B)
{
	return (A.x * B.x) + (A.y * B.y);
}

class WaveGenerator
{
	XMFLOAT3* m_pGridVertices = nullptr;
	WaveVertex* m_pWaveVertices = nullptr;
	unsigned int* m_pIndices = nullptr;
	unsigned int m_NumIndices = 0;
	clock_t timer = clock();
	double GetTime();

	int N = -1;
	int NMinus1 = -1;
	int length = -1;
	XMFLOAT2 w = XMFLOAT2(15.0f, .0f);

	// universal constants
	float g = 9.8f;

	// Wave properties
	float Q = 0.5f;
	XMFLOAT3 D = {0.9f, .0f, 0.2f};
	float L = 100.141593f;
	float A = 0.0001f;
	float W = 10.28f;
	float S = 1.f;
	float phi = 5;

	// FFT variables
	FFTData* fft_data = nullptr;
	std::complex<float> *h_tilde = nullptr,
		*h_tilde_slopex = nullptr, *h_tilde_slopez = nullptr,
		*h_tilde_dx = nullptr, *h_tilde_dz = nullptr;

	FFT* fft = nullptr;

	// FFT functions
	std::complex<float> hTilde(float t, int n_prime, int m_prime) const;
	std::complex<float> hTilde_0(int n_prime, int m_prime);
	float phillips(int n_prime, int m_prime);
	float dispersion(int j, int i) const;

public:
	WaveGenerator();
	~WaveGenerator();

	void GenerateGrid(int N, float length);

	const unsigned int* GetVertIndices() const
	{
		return m_pIndices;
	}

	const unsigned int GetNumVertices() const
	{
		return N * N;
	}

	const unsigned int GetNumIndices() const
	{
		return m_NumIndices;
	}

	const WaveVertex* GetWave();

	// New wave
	const WaveVertex* GetWaveFFT();
};
