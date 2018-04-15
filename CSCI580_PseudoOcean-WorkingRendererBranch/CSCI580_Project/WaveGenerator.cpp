#include "./WaveGenerator.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>
#include <iostream>
#undef _USE_MATH_DEFINES

#define INDEX i * N + j
#define H_INDEX i * (N - 1) + j
#define NUM_INDICES(x, y) (x - 1) * (y - 1) * 6

float uniformRandomVariable()
{
	return (float)rand() / RAND_MAX;
}

std::complex<float> gaussianRandomVariable()
{
	float x1, x2, w;
	do
	{
		x1 = 2.f * uniformRandomVariable() - 1.f;
		x2 = 2.f * uniformRandomVariable() - 1.f;
		w = x1 * x1 + x2 * x2;
	}
	while (w >= 1.f);
	w = sqrt((-2.f * log(w)) / w);
	return std::complex<float>(x1 * w, x2 * w);
}

float len(XMFLOAT2 vec)
{
	return sqrt(vec.x * vec.x + vec.y * vec.y);
}

XMFLOAT2 unit(XMFLOAT2 vec2)
{
	float l = len(vec2);
	return XMFLOAT2(vec2.x / l, vec2.y / l);
}

float len(XMFLOAT3 vec)
{
	return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

XMFLOAT3 unit(XMFLOAT3 vec3)
{
	float l = len(vec3);
	return XMFLOAT3(vec3.x / l, vec3.y / l, vec3.y / l);
}

WaveGenerator::WaveGenerator()
{

}

WaveGenerator::~WaveGenerator()
{
	if (m_pGridVertices)
		delete[] m_pGridVertices;
	if (m_pWaveVertices)
		delete[] m_pWaveVertices;

	m_pGridVertices = nullptr;
	m_pWaveVertices = nullptr;

	// Free fft related variables
	if (fft_data)
		delete[] fft_data;
	if (h_tilde)
		delete[] h_tilde;
	if (h_tilde_slopex)
		delete[] h_tilde_slopex;
	if (h_tilde_slopez)
		delete[] h_tilde_slopez;
	if (h_tilde_dx)
		delete[] h_tilde_dx;
	if (h_tilde_dz)
		delete[] h_tilde_dz;

	delete fft;

	fft_data = nullptr;
	h_tilde = nullptr;
	h_tilde_slopex = nullptr;
	h_tilde_slopez = nullptr;
	h_tilde_dx = nullptr;
	h_tilde_dz = nullptr;

	fft = nullptr;
}

void WaveGenerator::GenerateGrid(int N, float length)
{
	this->N = N;
	this->NMinus1 = N - 1;
	this->length = length;

	// FFT data initialization
	fft = new FFT(NMinus1);
	h_tilde = new std::complex<float>[NMinus1 * NMinus1];
	h_tilde_slopex = new std::complex<float>[NMinus1 * NMinus1];
	h_tilde_slopez = new std::complex<float>[NMinus1 * NMinus1];
	h_tilde_dx = new std::complex<float>[NMinus1 * NMinus1];
	h_tilde_dz = new std::complex<float>[NMinus1 * NMinus1];
	fft_data = new FFTData[N * N];

	m_pGridVertices = new XMFLOAT3[N * N];
	m_pWaveVertices = new WaveVertex[N * N];
	m_NumIndices = NUM_INDICES(N, N);
	m_pIndices = new unsigned int[m_NumIndices];

	float xOffset = length / static_cast<float>(N);
	float zOffset = length / static_cast<float>(N);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			int index = i * N + j;

			std::complex<float> htilde0 = hTilde_0(j, i);
			std::complex<float> htilde0mk_conj = conj(hTilde_0(-j, -i));

			fft_data[INDEX].htilde0.x = htilde0.real();
			fft_data[INDEX].htilde0.y = htilde0.imag();
			fft_data[INDEX].htilde0mk.x = htilde0mk_conj.real();
			fft_data[INDEX].htilde0mk.y = htilde0mk_conj.imag();


			m_pGridVertices[INDEX].x = m_pWaveVertices[INDEX].position.x = ((j - NMinus1) / 2.0f) * length / NMinus1;
			m_pGridVertices[INDEX].y = m_pWaveVertices[INDEX].position.y = 0.0f;
			m_pGridVertices[INDEX].z = m_pWaveVertices[INDEX].position.z = ((i - NMinus1) / 2.0f) * length / NMinus1;

			m_pWaveVertices[INDEX].normal.x = 0.0f;
			m_pWaveVertices[INDEX].normal.y = 1.0f;
			m_pWaveVertices[INDEX].normal.z = 0.0f;
		}
	}

	unsigned int indexCount = 0;
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = 0; j < N - 1; j++)
		{
			m_pIndices[indexCount++] = static_cast<unsigned int>(i * N + j);
			m_pIndices[indexCount++] = static_cast<unsigned int>((i + 1) * N + j);
			m_pIndices[indexCount++] = static_cast<unsigned int>(i * N + j + 1);
			m_pIndices[indexCount++] = static_cast<unsigned int>((i + 1) * N + j);
			m_pIndices[indexCount++] = static_cast<unsigned int>((i + 1) * N + 1 + j);
			m_pIndices[indexCount++] = static_cast<unsigned int>(i * N + j + 1);
		}
	}

	assert(indexCount == m_NumIndices, "The index count needs to be the same as the number of indices");
}

double WaveGenerator::GetTime()
{
	clock_t time2 = clock();
	clock_t timediff = time2 - timer;

	return ((double)timediff) / CLOCKS_PER_SEC;
}

const WaveVertex* WaveGenerator::GetWave()
{
	double t = GetTime();
	return m_pWaveVertices;

	if (!m_pGridVertices || !m_pWaveVertices)
		return nullptr;

	t = GetTime();

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			float x_input = m_pGridVertices[INDEX].x;
			float z_input = m_pGridVertices[INDEX].z;

			XMFLOAT3 x_dir = {x_input, 0, 0};
			XMFLOAT3 z_dir = {0, 0, z_input};
			double cos_sin_arg_input = (dot(D, m_pGridVertices[INDEX]) * W) + (phi * t);

			m_pWaveVertices[INDEX].position.x = x_input + (Q * A + (dot(D, x_dir) * cos(cos_sin_arg_input)));
			m_pWaveVertices[INDEX].position.y = A * sin(cos_sin_arg_input);
			m_pWaveVertices[INDEX].position.z = z_input + (Q * A + (dot(D, z_dir) * cos(cos_sin_arg_input)));

			double cos_sin_arg_normal = (dot(D, m_pWaveVertices[INDEX].position) * W) + (phi * t);
			m_pWaveVertices[INDEX].normal.x = - (dot(D, x_dir) * W * A * cos(cos_sin_arg_normal));
			m_pWaveVertices[INDEX].normal.y = 1 - (Q * W * A * sin(cos_sin_arg_normal));
			m_pWaveVertices[INDEX].normal.z = - (dot(D, z_dir) * W * A * cos(cos_sin_arg_normal));
		}
	}
	return m_pWaveVertices;
}

const WaveVertex* WaveGenerator::GetWaveFFT()
{
//	return m_pWaveVertices;

	if (!m_pGridVertices || !m_pWaveVertices)
		return nullptr;

	double t = GetTime();

	float lambda = -1.0f;

	for (int i = 0; i < NMinus1; ++i)
	{
		float kz = M_PI * (2 * i - NMinus1) / length;
		for (int j = 0; j < NMinus1; ++j)
		{
			float kx = M_PI * (2 * j - NMinus1) / length;
			float len = sqrt(kx * kx + kz * kz);

			h_tilde[H_INDEX] = hTilde(t, j, i);
			h_tilde_slopex[H_INDEX] = h_tilde[H_INDEX] * std::complex<float>(0, kx);
			h_tilde_slopez[H_INDEX] = h_tilde[H_INDEX] * std::complex<float>(0, kz);
			if (len < 0.000001f)
			{
				h_tilde_dx[H_INDEX] = std::complex<float>(0.0f, 0.0f);
				h_tilde_dz[H_INDEX] = std::complex<float>(0.0f, 0.0f);
			}
			else
			{
				h_tilde_dx[H_INDEX] = h_tilde[H_INDEX] * std::complex<float>(0, -kx / len);
				h_tilde_dz[H_INDEX] = h_tilde[H_INDEX] * std::complex<float>(0, -kz / len);
			}
		}
	}

	for (int i = 0; i < NMinus1; ++i)
	{
		fft->fft(h_tilde, h_tilde, 1, i * NMinus1);
		fft->fft(h_tilde_slopex, h_tilde_slopex, 1, i * NMinus1);
		fft->fft(h_tilde_slopez, h_tilde_slopez, 1, i * NMinus1);
		fft->fft(h_tilde_dx, h_tilde_dx, 1, i * NMinus1);
		fft->fft(h_tilde_dz, h_tilde_dz, 1, i * NMinus1);
	}

	for (int j = 0; j < NMinus1; ++j)
	{
		fft->fft(h_tilde, h_tilde, NMinus1, j);
		fft->fft(h_tilde_slopex, h_tilde_slopex, NMinus1, j);
		fft->fft(h_tilde_slopez, h_tilde_slopez, NMinus1, j);
		fft->fft(h_tilde_dx, h_tilde_dx, NMinus1, j);
		fft->fft(h_tilde_dz, h_tilde_dz, NMinus1, j);
	}

	float signs[] = {1.0f, -1.0f};
	for (int i = 0; i < NMinus1; i++)
	{
		for (int j = 0; j < NMinus1; j++)
		{
			int sign = signs[(j + i) & 1];

			h_tilde[H_INDEX] = complexScalarMultiply(h_tilde[H_INDEX], sign);

			// height
			m_pWaveVertices[INDEX].position.y = h_tilde[H_INDEX].real();

			// displacement
			h_tilde_dx[H_INDEX] = complexScalarMultiply(h_tilde_dx[H_INDEX], sign);
			h_tilde_dz[H_INDEX] = complexScalarMultiply(h_tilde_dz[H_INDEX], sign);
			m_pWaveVertices[INDEX].position.x = m_pGridVertices[INDEX].x + h_tilde_dx[H_INDEX].real() * lambda;
			m_pWaveVertices[INDEX].position.z = m_pGridVertices[INDEX].z + h_tilde_dz[H_INDEX].real() * lambda;

			// normal
			h_tilde_slopex[H_INDEX] = complexScalarMultiply(h_tilde_slopex[H_INDEX], sign);
			h_tilde_slopez[H_INDEX] = complexScalarMultiply(h_tilde_slopez[H_INDEX], sign);
			XMFLOAT3 n = {0.0f - h_tilde_slopex[H_INDEX].real(), 1.0f, 0.0f - h_tilde_slopez[H_INDEX].real()};
			//.unit();
			m_pWaveVertices[INDEX].normal.x = n.x;
			m_pWaveVertices[INDEX].normal.y = n.y;
			m_pWaveVertices[INDEX].normal.z = n.z;

			// for tiling
			if (j == 0 && i == 0)
			{
				// N = Nplus1
				// N
				int tilingIndex = INDEX + NMinus1 + N * NMinus1;
				m_pWaveVertices[tilingIndex].position.y = h_tilde[H_INDEX].real();
				m_pWaveVertices[tilingIndex].position.x = m_pGridVertices[tilingIndex].x + h_tilde_dx[H_INDEX].real() * lambda;
				m_pWaveVertices[tilingIndex].position.z = m_pGridVertices[tilingIndex].z + h_tilde_dz[H_INDEX].real() * lambda;

				m_pWaveVertices[tilingIndex].normal.x = n.x;
				m_pWaveVertices[tilingIndex].normal.y = n.y;
				m_pWaveVertices[tilingIndex].normal.z = n.z;
			}
			if (j == 0)
			{
				int tilingIndex = INDEX + NMinus1;
				m_pWaveVertices[tilingIndex].position.y = h_tilde[H_INDEX].real();
				m_pWaveVertices[tilingIndex].position.x = m_pGridVertices[tilingIndex].x + h_tilde_dx[H_INDEX].real() * lambda;
				m_pWaveVertices[tilingIndex].position.z = m_pGridVertices[tilingIndex].z + h_tilde_dz[H_INDEX].real() * lambda;

				m_pWaveVertices[tilingIndex].normal.x = n.x;
				m_pWaveVertices[tilingIndex].normal.y = n.y;
				m_pWaveVertices[tilingIndex].normal.z = n.z;
			}
			if (i == 0)
			{
				int tilingIndex = INDEX + N * NMinus1;
				m_pWaveVertices[tilingIndex].position.y = h_tilde[H_INDEX].real();
				m_pWaveVertices[tilingIndex].position.x = m_pGridVertices[tilingIndex].x + h_tilde_dx[H_INDEX].real() * lambda;
				m_pWaveVertices[tilingIndex].position.z = m_pGridVertices[tilingIndex].z + h_tilde_dz[H_INDEX].real() * lambda;

				m_pWaveVertices[tilingIndex].normal.x = n.x;
				m_pWaveVertices[tilingIndex].normal.y = n.y;
				m_pWaveVertices[tilingIndex].normal.z = n.z;
			}
		}
	}

	return m_pWaveVertices;
}

std::complex<float> WaveGenerator::hTilde(float t, int j, int i) const
{
	std::complex<float> htilde0(fft_data[INDEX].htilde0.x, fft_data[INDEX].htilde0mk.y);
	std::complex<float> htilde0mkconj(fft_data[INDEX].htilde0mk.x, fft_data[INDEX].htilde0mk.y);

	float omegat = dispersion(j, i) * t;

	float cos_ = cos(omegat);
	float sin_ = sin(omegat);

	std::complex<float> c0(cos_, sin_);
	std::complex<float> c1(cos_, -sin_);

	return htilde0 * c0 + htilde0mkconj * c1;
}

float WaveGenerator::phillips(int n_prime, int m_prime)
{
	XMFLOAT2 k(M_PI * (2 * n_prime - N) / length, M_PI * (2 * m_prime - N) / length);
	float k_length = len(k);
	if (k_length < 0.000001) return 0.0;

	float k_length2 = k_length * k_length;
	float k_length4 = k_length2 * k_length2;

	float k_dot_w = dot(unit(k), unit(w));
	float k_dot_w2 = k_dot_w * k_dot_w * k_dot_w * k_dot_w * k_dot_w * k_dot_w;

	float w_length = len(w);
	float L = w_length * w_length / g;
	float L2 = L * L;

	float damping = 0.001;
	float l2 = L2 * damping * damping;

	return A * exp(-1.0f / (k_length2 * L2)) / k_length4 * k_dot_w2 * exp(-k_length2 * l2);
}

std::complex<float> WaveGenerator::hTilde_0(int n_prime, int m_prime)
{
	std::complex<float> r = gaussianRandomVariable();
	return r * sqrt(phillips(n_prime, m_prime) / 2.0f);
}

float WaveGenerator::dispersion(int n_prime, int m_prime) const
{
	float w_0 = 2.0f * M_PI / 200.0f;
	float kx = M_PI * (2 * n_prime - N + 1) / length;
	float kz = M_PI * (2 * m_prime - N + 1) / length;
	return floor(sqrt(g * sqrt(kx * kx + kz * kz)) / w_0) * w_0;
}
