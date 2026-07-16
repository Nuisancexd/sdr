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
#define SAMPLE_RATE 10000000
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

    //int blocks = SAMPLE_RATE / FFT_size;
    double wk;
    int blocks = 1000;
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
    }

    free(polyph_1);
    free(polyph_2);
    free(rx1);
    free(rx2);
    free(w);    
    sdr::free_config(&sdr);
    return EXIT_SUCCESS;
}