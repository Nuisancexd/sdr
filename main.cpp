#include <fftw3.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <csignal>
#include "types.h"
#include "sdr.h"
#include "math.h"
#include "FFT.h"
#include <cstring>

#define DEBUG
//#define SUM_SAMPL_SPEC

#define M 1000000
static ssize_t FREQ = 2437000000;
//static ssize_t FREQ = 2412000000;
#define SAMPLE_RATE 20000000
//#define SAMPLE_RATE 61440000
#define FFT_size 1024

sig_atomic_t doneman = 1;
void signal_handler(int)
{
    doneman = 0;
}

int main()
{
    CONFIG sdr = {};
#ifndef DEBUG
    if(sdr::init_sdr(&sdr, "ip:192.168.2.1", FREQ, SAMPLE_RATE) && sdr::free_config(&sdr))
        return EXIT_FAILURE;
#endif
    
    PCOMPLEX rx1 = (PCOMPLEX)malloc(FFT_size * sizeof(COMPLEX));
    PCOMPLEX rx2 = (PCOMPLEX)malloc(FFT_size * sizeof(COMPLEX));
    FFT_t fft1;
    FFT_t fft2;
    if(!FFT::fft_init(&fft1, FFT_size) || !FFT::fft_init(&fft2, FFT_size))
    { printf("failed fft init\n"); return EXIT_FAILURE; }

    PCOMPLEX fft_sum_1 = (PCOMPLEX)malloc(FFT_size * sizeof(COMPLEX));
    PCOMPLEX fft_sum_2 = (PCOMPLEX)malloc(FFT_size * sizeof(COMPLEX));
    fftwf_complex* fft = (fftwf_complex*)fftwf_malloc(FFT_size * sizeof(fftwf_complex));
    int padding = 2048;
    int size = FFT_size * padding;
    fftwf_complex* padded = (fftwf_complex*)fftwf_malloc(size * sizeof(fftwf_complex));
    fftwf_complex* OBPF = (fftwf_complex*)fftwf_malloc(size * sizeof(fftwf_complex));
    fftwf_plan OBPF_plan = fftwf_plan_dft_1d(size, padded, OBPF, FFTW_BACKWARD, FFTW_ESTIMATE);

    unsigned int size_padd;
    long double maxPower;
    int peak;
    long double p;
    double phase_rad;
    double s;
    double angle;
    double threshold;

    float lambda = 299792458.0f / FREQ;
    float d = lambda / 2.0f;
    float tau = (float)d / 299792458.0f;
    
    int blocks = 100;

    std::signal(SIGINT, signal_handler);
    while(doneman)
    {
        memset(fft, 0x00, FFT_size * sizeof(COMPLEX));
        memset(fft_sum_1, 0x00, FFT_size * sizeof(COMPLEX));
        memset(fft_sum_2, 0x00, FFT_size * sizeof(COMPLEX));

#ifdef SUM_SAMPL
        for(int i = 0; i < blocks; ++i)
        {
            sdr::sdr_receive(&sdr, rx1, rx2, FFT_size);
            for(int j = 0; j < FFT_size; ++j)
            {
                fft_sum_1[j].i += rx1[j].i;
                fft_sum_1[j].q += rx1[j].q;
                fft_sum_2[j].i += rx2[j].i;
                fft_sum_2[j].q += rx2[j].q;
            }
        }
        FFT::fft_exec(&fft1, fft_sum_1, FFT_size);
        FFT::fft_exec(&fft2, fft_sum_2, FFT_size);
        for(int i = 0; i < FFT_size; ++i)
        {
            fft[i][0] = fft1.out[i][0] * fft2.out[i][0] + fft1.out[i][1] * fft2.out[i][1];
            fft[i][1] = fft1.out[i][1] * fft2.out[i][0] - fft1.out[i][0] * fft2.out[i][1];
        }
#elif defined SUM_SAMPL_SPEC
        for (int i = 0; i < blocks; ++i)
        {
            if (!sdr::sdr_receive(&sdr, rx2, rx1, FFT_size))
                break;
        
            FFT::fft_exec(&fft1, rx1, FFT_size);
            FFT::fft_exec(&fft2, rx2, FFT_size);
        
            for (int j = 0; j < FFT_size; ++j)
            {
                fft[j][0] += fft2.out[j][0] * fft1.out[j][0] + fft2.out[j][1] * fft1.out[j][1];
                fft[j][1] += fft2.out[j][1] * fft1.out[j][0] - fft2.out[j][0] * fft1.out[j][1];
            }
        }
#else

        for(double i = -0.205; i <= 0.215; i += 0.01)
        {

        memset(fft, 0x00, FFT_size * sizeof(COMPLEX));
        float delay_tau = i * 1e-9;
        float phase;
        float phase_delay;

        /* f_bin = m * fs/N */
        float f0 = (30.0f * (float)SAMPLE_RATE/FFT_size);

        /*
        exp(j*2M_PI * f0 * n/fs)
        exp(j*2M_PI * f0 * (n/fs - delay))
        */            
        for (int i = 0; i < blocks; ++i)
        {
            memset(rx1, 0x00, FFT_size * sizeof(COMPLEX));
            memset(rx2, 0x00, FFT_size * sizeof(COMPLEX));
            for(int j = 0; j < FFT_size; ++j)
            {
                for(int k = 1; k <= 5; ++k)
                {
                    phase = 2.0f * M_PI * (f0 + (float)k * 2.0f) * (float)j / SAMPLE_RATE;
                    rx1[j].i += cosf(phase);
                    rx1[j].q += sinf(phase);

                    //phase_delay = 2.0f * M_PI * f0 * ((float)n / SAMPLE_RATE - delay_samples);
                    phase_delay = 2.0f * M_PI * (f0 + (float)k * 2.0f) * ((float)j / SAMPLE_RATE) + (-2.0f * M_PI * FREQ * delay_tau);
                    rx2[j].i += cosf(phase_delay);
                    rx2[j].q += sinf(phase_delay);
                }
            }
            FFT::fft_exec(&fft1, rx1, FFT_size);
            FFT::fft_exec(&fft2, rx2, FFT_size);

            for (int j = 0; j < FFT_size; ++j)
            {
                fft[j][0] += fft2.out[j][0] * fft1.out[j][0] + fft2.out[j][1] * fft1.out[j][1];
                fft[j][1] += fft2.out[j][1] * fft1.out[j][0] - fft2.out[j][0] * fft1.out[j][1];
            }
        }
#endif

        memset(padded, 0x00, size * sizeof(fftwf_complex));
        memset(OBPF, 0x00, size * sizeof(fftwf_complex));
        
        for(int i = 0; i <= FFT_size / 2; ++i)
        {
            padded[i][0] = fft[i][0];
            padded[i][1] = fft[i][1];
        }
        for(int i = FFT_size/2 + 1; i < FFT_size; ++i)
        {
            padded[size - (FFT_size - i)][0] = fft[i][0];
            padded[size - (FFT_size - i)][1] = fft[i][1];
        }

        fftwf_execute(OBPF_plan);

        size_padd = (unsigned int)ceilf((double)tau * (double)SAMPLE_RATE * padding) + 2;

        maxPower = 0;
        peak = 0;
        for(int i = 0; i <= size_padd; ++i)
        {
            p = OBPF[i][0]*OBPF[i][0] + OBPF[i][1]*OBPF[i][1];
            if(p > maxPower) { maxPower = p; peak = i; }
        }
        for(int i = size - size_padd; i < size; ++i)
        {
            p = OBPF[i][0]*OBPF[i][0] + OBPF[i][1]*OBPF[i][1];
            if(p > maxPower) { maxPower = p; peak = i; }
        }
        
        phase_rad = atan2f(OBPF[peak][1],OBPF[peak][0]);
#ifdef DEBUG
        double time_delay = phase_rad / (2.0 * M_PI * FREQ);
#else
        phase_rad -= 3.62f;
#endif
        if (phase_rad > M_PI)  phase_rad -= 2.0f * M_PI;
        if (phase_rad < -M_PI) phase_rad += 2.0f * M_PI;

        s = phase_rad * lambda/ (2.0 * M_PI * d);
        if(s > 1.0) s = 1.0;
        if(s < -1.0) s = - 1.0;
        angle = asinf(s) * 180.0 / M_PI;
#ifndef DEBUG
        printf("%-14s  phase %+7.2f rad   Angle %+8.2f\n", " ", phase_rad, angle);
#else
        
        printf("peak %8d  phase %+7.4f rad  angle %+7.2f delay_peak_phase %+7.2f ns  delay_time%+7.2f ns  delay_samp %.3f\n", 
            peak, phase_rad, angle, time_delay * 1e9, delay_tau * 1e9f, delay_tau * SAMPLE_RATE);
        } 
        break;
#endif
    }

    if(fft_sum_1)
        free(fft_sum_1);
    if(fft_sum_2)
        free(fft_sum_2);
    if(rx1)
        free(rx1);
    if(rx2)
        free(rx2);
    sdr::free_config(&sdr);
    fftwf_free(padded);
    fftwf_free(OBPF);
    fftwf_destroy_plan(OBPF_plan);
    return EXIT_SUCCESS;
}