#include "FFT.h"

#include "math.h"
#include "types.h"
#include <cstdio>
#include <cstring>
#include <fftw3.h>


bool FFT::fft_init(PFFT fft, int FFT_size)
{
    fft->in = (fftwf_complex*)fftwf_malloc(FFT_size * sizeof(fftwf_complex));
    fft->out = (fftwf_complex*)fftwf_malloc(FFT_size * sizeof(fftwf_complex));
    if(!fft->in || !fft->out)
        return false;

    fft->plan = fftwf_plan_dft_1d(FFT_size, fft->in, fft->out, FFTW_FORWARD, FFTW_ESTIMATE);
    if(!fft->plan)
        return false;

    return true;
}

void FFT::fft_exec(PFFT fft, PCOMPLEX in, int FFT_size)
{
    for(int i = 0; i < FFT_size; ++i)
    {
        fft->in[i][0] = in[i].i;
        fft->in[i][1] = in[i].q;
    }
    fftwf_execute(fft->plan);
}

float FFT::fft_amplitude(PFFT fft, int m)
{
    float re = fft->out[m][0];
    float im = fft->out[m][1];

    return sqrtf(re*re + im*im);
}

float FFT::fft_phase(PFFT fft, int m)
{
    return atan2f( fft->out[m][1], fft->out[m][0]);
}

float FFT::fft_re(PFFT fft, int m)
{
    return fft->out[m][0];
}

float FFT::fft_im(PFFT fft, int m)
{
    return fft->out[m][1];
}

void FFT::phase_diff(PHASE *phase, PFFT rx1, PFFT rx2, int FFT_size)
{
    float rx1_re;
    float rx1_im;
    float rx2_re;
    float rx2_im;

    for(int m = 0; m < FFT_size; ++m)
    {
        rx1_re = rx1->out[m][0];
        rx1_im = -rx1->out[m][1];
        rx2_re = rx2->out[m][0];
        rx2_im = rx2->out[m][1];

        /* 
        rx1_conj * rx2        
        (rx1_re + rx1_im) * (rx2_re * rx2_im)
        rx1_re * rx2_re + 
        rx1_re * rx2_im +
        rx1_im * rx2_re +
        rx1_im * rx2_im
        */
        phase->bin[m].re = rx1_re * rx2_re + rx1_im * rx2_im;
        phase->bin[m].im = rx1_re * rx2_im - rx1_im * rx2_re; 

        phase->bin[m].phase = atan2f(phase->bin[m].im, phase->bin[m].re);
        phase->bin[m].amplitude = sqrtf(phase->bin[m].re * phase->bin[m].re + phase->bin[m].im * phase->bin[m].im);
    }
}


double* FFT::window_Hamming_init(int FFT_size)
{
    double* window = (double*)malloc(FFT_size * sizeof(double));
    memset(window, 0x00, FFT_size * sizeof(double));
    double tPI = 2 * M_PI;
    double N = FFT_size - 1;
    for(int i = 0; i < FFT_size; ++i)
    {
        window[i] = 0.54 - 0.46*cosf((tPI * i)/N);
    }
    return window;
}