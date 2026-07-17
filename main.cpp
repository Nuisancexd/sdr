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

#define FREQ 2412000000
#define SAMPLE_RATE 20000000
#define FFT_size 1024

sig_atomic_t doneman = 1;
void signal_handler(int)
{
    doneman = 0;
}

int main()
{
    PCOMPLEX rx1 = (PCOMPLEX)malloc(FFT_size * sizeof(COMPLEX));
    PCOMPLEX rx2 = (PCOMPLEX)malloc(FFT_size * sizeof(COMPLEX));
    CONFIG sdr = {};
    if(sdr::init_sdr(&sdr, "ip:192.168.2.1", FREQ, SAMPLE_RATE) && sdr::free_config(&sdr))
        return EXIT_FAILURE;
    
    PCOMPLEX polyph_1 = (PCOMPLEX)malloc(FFT_size * sizeof(COMPLEX));
    PCOMPLEX polyph_2 = (PCOMPLEX)malloc(FFT_size * sizeof(COMPLEX));
    PHASE phase;
    FFT_t fft1;
    FFT_t fft2;

    if(!FFT::fft_init(&fft1, FFT_size) || !FFT::fft_init(&fft2, FFT_size))
    { printf("failed fft init\n"); return EXIT_FAILURE; }

    float lambda = 299792458.0f / FREQ;
    float d = lambda / 2.0f;
    int max_bin;
    float max_ampl, delta_phi, sin_theta, theta_rad, theta_deg;

    //int blocks = SAMPLE_RATE / FFT_size;
    double wk;
    int blocks = 500;
    int win_len = blocks * FFT_size;
    double* w = FFT::window_Hamming_init(win_len);
    std::signal(SIGINT, signal_handler);
    while(doneman)
    {
        memset(polyph_1, 0x00, FFT_size * sizeof(COMPLEX));
        memset(polyph_2, 0x00, FFT_size * sizeof(COMPLEX));
        for(int i = 0; i < blocks; ++i)
        {
            if(!sdr::sdr_receive(&sdr, rx1, rx2, FFT_size))
                break;

            for(int k = 0; k < FFT_size; ++k)
            {
                wk = w[k + i * FFT_size]; 
                polyph_1[k].i += rx1[k].i * wk;
                polyph_1[k].q += rx1[k].q * wk;
                polyph_2[k].i += rx2[k].i * wk;
                polyph_2[k].q += rx2[k].q * wk;
            }
        }

        FFT::fft_exec(&fft1, polyph_1, FFT_size);
        FFT::fft_exec(&fft2, polyph_2, FFT_size);

        FFT::phase_diff(&phase, &fft1, &fft2, FFT_size);

        max_ampl = 0.0f;
        max_bin = 0;
        for (int m = 1; m < FFT_size; ++m)
        {
            float amplitude = FFT::fft_amplitude(&fft1, m) + FFT::fft_amplitude(&fft2, m);
            if (amplitude > max_ampl)
            {
                max_ampl = amplitude;
                max_bin = m;
            }
        }

        static float PHASE_OFFSET = -2.7f;
        delta_phi = phase.bin[max_bin].phase;
        delta_phi -= PHASE_OFFSET;

        delta_phi = fmodf(delta_phi + M_PI, 2.0 * M_PI);
        if (delta_phi < 0)
            delta_phi += 2.0 * M_PI;
        delta_phi -= M_PI;

        sin_theta = (delta_phi * lambda) / (2.0f * M_PI * d);
        if (sin_theta >  1.0f) sin_theta =  1.0f;
        if (sin_theta < -1.0f) sin_theta = -1.0f;

        theta_rad = asinf(sin_theta);
        theta_deg = theta_rad * 180.0f / M_PI;

        printf("Signal at bin %d, amplitude = %.1f\n", max_bin, max_ampl / FFT_size);
        printf("Phase difference = %.3f rad (%.1f deg)\n", delta_phi, delta_phi * 180.0f / M_PI);
        printf("angle = %.1f degrees\n", theta_deg);
    }

    free(polyph_1);
    free(polyph_2);
    free(rx1);
    free(rx2);
    free(w);    
    sdr::free_config(&sdr);
    return EXIT_SUCCESS;
}