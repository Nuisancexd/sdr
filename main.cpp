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

#define M 1000000

//#define FREQ 2412000000
static ssize_t FREQ = 2412000000;
#define BASEBAND 20 * M
#define SAMPLE_RATE 20000000
#define FFT_size 1024

constexpr float TIME_SIGNAL = (float)FFT_size/(float)SAMPLE_RATE;
float TIME_SIGNAL_WIFI = 0.02 / TIME_SIGNAL;
int WIFI_INTERVAL = (int)TIME_SIGNAL_WIFI;

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
    
    PCOMPLEX fft_sum_1 = (PCOMPLEX)malloc(FFT_size * sizeof(COMPLEX));
    PCOMPLEX fft_sum_2 = (PCOMPLEX)malloc(FFT_size * sizeof(COMPLEX));
    PHASE phase;
    FFT_t fft1;
    FFT_t fft2;

    // iio_buffer* rx_scan_buf = NULL;
    // if(!sdr::create_buffer(&sdr, &rx_scan_buf, 1 << 16))
    //     { printf("failed create buffer scan\n"); return EXIT_FAILURE; }
    PCOMPLEX rx_scan = (PCOMPLEX)malloc(FFT_size * WIFI_INTERVAL * sizeof(COMPLEX));
    int samples_scan = FFT_size * WIFI_INTERVAL;
    float power_scan = 0;
    float power_scan_db;
    float porog_H1_db = 28.0;
    int ch = 1;
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
        for(; doneman; power_scan = 0)
        {
            sdr::change_channel_freq(&sdr, FREQ + ((ch - 1) * 5 * M));
            usleep(10000);
            sdr.buf_pos = sdr.buf_end;
            for(int i = 0; i < 2; ++i)
                sdr::sdr_receive_one_channel(&sdr, sdr.rxbuf, rx_scan, FFT_size);    
            
            sdr::sdr_receive_one_channel(&sdr, sdr.rxbuf, rx_scan, samples_scan);
                
            for(int j = 0; j < samples_scan; ++j)
                power_scan += (rx_scan[j].i * rx_scan[j].i) + (rx_scan[j].q * rx_scan[j].q);

            power_scan /= samples_scan;
            power_scan_db = 10.0 * log10(power_scan + 1e-10);
            printf("Channel %d\t avg_sum_iq %f\t power_db %f\n", ch, power_scan, power_scan_db);
            if(power_scan_db >= porog_H1_db)
            {
                printf("Peleng for channel %d power_signal >= porog_H1\n", ch);
                ++ch;
                power_scan = 0;
                break;
            }
            if(++ch > 13)
                ch = 1;
        }

        memset(fft_sum_1, 0x00, FFT_size * sizeof(COMPLEX));
        memset(fft_sum_2, 0x00, FFT_size * sizeof(COMPLEX));
        for(int i = 0; i < blocks; ++i)
        {
            if(!sdr::sdr_receive(&sdr, rx1, rx2, FFT_size))
                break;

            for(int k = 0; k < FFT_size; ++k)
            {
                wk = w[k + i * FFT_size]; 
                fft_sum_1[k].i += rx1[k].i * wk;
                fft_sum_1[k].q += rx1[k].q * wk;
                fft_sum_2[k].i += rx2[k].i * wk;
                fft_sum_2[k].q += rx2[k].q * wk;
            }
        }

        FFT::fft_exec(&fft1, fft_sum_1, FFT_size);
        FFT::fft_exec(&fft2, fft_sum_2, FFT_size);

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
        printf("angle = %.1f degrees\n\n", theta_deg);
    }

    if(fft_sum_1)
        free(fft_sum_1);
    if(fft_sum_2)
        free(fft_sum_2);
    if(rx1)
        free(rx1);
    if(rx2)
        free(rx2);
    if(w)
        free(w);    
    if(rx_scan)
        free(rx_scan);
    //sdr::free_buf_rx(&rx_scan_buf);
    sdr::free_config(&sdr);
    return EXIT_SUCCESS;
}