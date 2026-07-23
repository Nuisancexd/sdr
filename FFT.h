#ifndef _FFT_H_
#define _FFT_H_

#include "types.h"
#include <fftw3.h>

typedef struct fft_t
{
    fftwf_complex* in;
    fftwf_complex* out;
    fftwf_plan plan;
} FFT_t, *PFFT;

typedef struct phase_bin_t
{
    float re;
    float im;

    float phase;
    float amplitude;

} PHASE_BIN, *PPHASE_BIN;

typedef struct phase_t
{
    PHASE_BIN bin[1024];

} PHASE, *PPHASE;

namespace FFT
{
    void phase_diff(PHASE *phase, PFFT rx1, PFFT rx2, int FFT_size);
    void fft_exec(PFFT fft, PCOMPLEX in, int FFT_size);
    bool fft_init(PFFT fft, int FFT_size);
    float fft_amplitude(PFFT fft, int m);
    float fft_phase(PFFT fft, int m);
    float fft_re(PFFT fft, int m);
    float fft_im(PFFT fft, int m);
    double* window_Hamming_init(int N);
}

#endif