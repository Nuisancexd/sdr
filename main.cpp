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
#include <vector>
#include <algorithm>

#include <matplot/matplot.h>
using namespace matplot;

//#define DEBUG

#define M 1000000

static ssize_t FREQ = 2440000000;
//static ssize_t FREQ = 2412000000;
//static ssize_t FREQ = 2462000000;
#define SAMPLE_RATE 20000000
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
    
    fftwf_complex* fft = (fftwf_complex*)fftwf_malloc(FFT_size * sizeof(fftwf_complex));
    float* phase_fft = (float*)malloc(FFT_size * sizeof(float));
    float* amplitude = (float*)malloc(FFT_size * sizeof(float));
    int size_padd;
    long double maxPower;
    double phase_rad;
    double s;
    double angle;
    double threshold;

    float speed_light = 299792458.0f;
    float lambda =  speed_light / FREQ;
    float d = lambda / 2.0f;
    float tau = d / speed_light;
    int blocks = 0;
    float scale = (blocks + 1); //* FFT_size * FFT_size;
    //float dt_s = 1.0f/(SAMPLE_RATE * padding);

    float max_correlation = -1.0f, correlation_angle;
    float bin = SAMPLE_RATE / TO_FLOAT(FFT_size);
    float freq_signal = 0.0f;

    std::vector<double> angle_vec(601);
    std::vector<double> angle_th(601);
    std::vector<std::vector<double>> phase_model_vec(601, std::vector<double>(FFT_size, 0.0f));
    int indx = 0;
    for (float ang = -90.0f; ang <= 90.0f; ang += 0.3f) 
    {
        angle_vec[indx] = ang;
        angle_th[indx++] = d * sinf(ang * M_PI/180.0f);
    }

    std::vector<double> corr_vec = {0.0};

    std::vector<double> time_axis = {-1.0, 1.0};
    std::vector<double> time_corr = {0.0, 1.0};

    auto fig_corr_time = matplot::figure(true);
    fig_corr_time->x_position(2000);
    fig_corr_time->y_position(500);
    fig_corr_time->size(512, 768);
    auto fig_time_ax = fig_corr_time->current_axes();
    fig_time_ax->clear();
    fig_time_ax->hold(matplot::on);
    fig_time_ax->xlabel("Задержка, нс");
    fig_time_ax->ylabel("Мощность");
    fig_time_ax->title("Корреляция во временной");
    fig_time_ax->ylim(matplot::automatic);
    auto plt_time = fig_time_ax->plot(time_axis, time_corr);
    plt_time->line_width(2);

    auto fig_corr = matplot::figure(true);
    fig_corr->x_position(1500);
    fig_corr->y_position(500);
    fig_corr->size(512, 768);
    auto fig_ax = fig_corr->current_axes();
    fig_ax->clear();
    fig_ax->hold(matplot::on);
    fig_ax->xlabel("Угол");
    fig_ax->ylabel("Корреляция");
    fig_ax->title("Корреляция");
    fig_ax->xlim({-90.0, 90.0});
    fig_ax->ylim({0.0, 1.1}); 
    auto plt = fig_ax->plot(angle_vec, corr_vec);
    plt->line_width(2);

    auto fig_polar = matplot::figure(true);
    fig_polar->x_position(1000);
    fig_polar->y_position(500);
    fig_polar->size(512, 768);
    auto fig_plr_ax = fig_polar->current_axes();
    fig_plr_ax->hold(matplot::on);
    fig_plr_ax->axis(false);
    fig_plr_ax->xlim({-1.2, 1.2});
    fig_plr_ax->ylim({-0.2, 1.2});
    auto plt_plr = fig_plr_ax->plot(angle_vec, corr_vec);
    plt_plr->line_width(2);
    std::vector<double> arc_x, arc_y;
    for (double a = -90.0; a <= 90.0; a += 2.0) 
    {
        arc_x.push_back(sinf(a * M_PI / 180.0));
        arc_y.push_back(cosf(a * M_PI / 180.0));
    }
    fig_plr_ax->plot(arc_x, arc_y)->color("lightgray");

    fig_plr_ax->plot(std::vector<double>{-1.0, 1.0}, std::vector<double>{0.0, 0.0})->color("black");
    fig_plr_ax->text(-1.15, 0.0, "-90");
    fig_plr_ax->text(1.05,  0.0, "90");
    fig_plr_ax->text(-0.88, 0.6, "-45");
    fig_plr_ax->text( 0.88, 0.6, "45");
    fig_plr_ax->text(-0.05, 1.05, "0");
    
    auto peleng_line_1 = fig_plr_ax->plot(std::vector<double>{0.0, 0.0}, std::vector<double>{0.0, 0.0});
    peleng_line_1->color("blue");
    peleng_line_1->line_width(3);
    auto peleng_line_2 = fig_plr_ax->plot(std::vector<double>{0.0, 0.0}, std::vector<double>{0.0, 0.0});
    peleng_line_2->color("blue");
    peleng_line_2->line_width(3);

    float* spectrum_acc = (float*)calloc(FFT_size, sizeof(float));
    std::vector<double> freq_axis(FFT_size);
    for (int i = 0; i < FFT_size; ++i)
    {
        int k = i - FFT_size / 2;
        freq_axis[i] = (FREQ + k * bin) / 1e6;
    }
    std::vector<double> spectrum_db(FFT_size, -120.0);
    auto fig_spectrum = matplot::figure(true);
    fig_spectrum->x_position(500);
    fig_spectrum->y_position(500);
    fig_spectrum->size(550*2, 768);
    auto fig_spec_ax = fig_spectrum->current_axes();
    fig_spec_ax->clear();
    fig_spec_ax->hold(matplot::on);
    fig_spec_ax->xlabel("Частота, МГц");
    fig_spec_ax->ylabel("Амплитуда, дБ");
    fig_spec_ax->title("Спектр");
    fig_spec_ax->xlim({(FREQ - SAMPLE_RATE/2.0)/1e6, (FREQ + SAMPLE_RATE/2.0)/1e6});
    fig_spec_ax->ylim(matplot::automatic);
    auto plt_spectrum = fig_spec_ax->plot(freq_axis, spectrum_db);
    plt_spectrum->line_width(2);

    std::vector<double> phase_deg(FFT_size, 0.0);
    auto fig_phase = matplot::figure(true);
    fig_phase->x_position(0);
    fig_phase->y_position(500);
    fig_phase->size(550*2, 768);
    auto fig_phase_ax = fig_phase->current_axes();
    fig_phase_ax->clear();
    fig_phase_ax->hold(matplot::on);
    fig_phase_ax->xlabel("Частота, МГц");
    fig_phase_ax->ylabel("Фаза, град");
    fig_phase_ax->title("ФЧХ");
    fig_phase_ax->xlim({(FREQ - SAMPLE_RATE/2.0)/1e6, (FREQ + SAMPLE_RATE/2.0)/1e6});
    fig_phase_ax->ylim({-200.0, 200});
    auto plt_phase = fig_phase_ax->plot(freq_axis, phase_deg);
    plt_phase->line_width(2);

    std::signal(SIGINT, signal_handler);
    while(doneman)
    {
#ifdef DEBUG
        for(double k = -0.2048; k <= 0.15 && doneman; k += 0.01)
        {
            //sleep(1);
            memset(fft, 0x00, FFT_size * sizeof(COMPLEX));
            memset(spectrum_acc, 0, FFT_size * sizeof(float));
            float delay_tau = k * 1e-9;
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
                    float time = (float)j / SAMPLE_RATE;
                    for(int t = -50; t <= 50; ++t)
                    {   
                        float subf = f0 + (float)t * bin;
                        phase = 2.0f * M_PI * subf * time;
                        rx1[j].i += cosf(phase);
                        rx1[j].q += sinf(phase);

                        phase_delay = 2.0f * M_PI * subf * time + (-2.0f * M_PI * FREQ * delay_tau);
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
                    spectrum_acc[j] += fft1.out[j][0] * fft1.out[j][0] + fft1.out[j][1] * fft1.out[j][1];   
                }
            }
#else
        sdr::sdr_receive(&sdr, rx1, rx2, FFT_size);
        float energy = 0.0f;
        for(int i = 0; i < FFT_size; ++i)
            energy += rx1[i].i * rx1[i].i + rx1[i].q * rx1[i].q;
        // printf("%f\n", energy);
        // continue;
        if((energy) < 500000.0f)
           continue;

        memset(fft, 0x00, FFT_size * sizeof(COMPLEX));
        memset(spectrum_acc, 0, FFT_size * sizeof(float));
        FFT::fft_exec(&fft1, rx1, FFT_size);
        FFT::fft_exec(&fft2, rx2, FFT_size);

        for (int j = 0; j < FFT_size; ++j)
        {
            fft[j][0] += fft2.out[j][0] * fft1.out[j][0] + fft2.out[j][1] * fft1.out[j][1];
            fft[j][1] += fft2.out[j][1] * fft1.out[j][0] - fft2.out[j][0] * fft1.out[j][1];
            spectrum_acc[j] += fft1.out[j][0] * fft1.out[j][0] + fft1.out[j][1] * fft1.out[j][1];
        }

        for(int i = 0; i < blocks; ++i)
        {
            sdr::sdr_receive(&sdr, rx1, rx2, FFT_size);
            FFT::fft_exec(&fft1, rx1, FFT_size);
            FFT::fft_exec(&fft2, rx2, FFT_size);
            for (int j = 0; j < FFT_size; ++j)
            {
                fft[j][0] += fft2.out[j][0] * fft1.out[j][0] + fft2.out[j][1] * fft1.out[j][1];
                fft[j][1] += fft2.out[j][1] * fft1.out[j][0] - fft2.out[j][0] * fft1.out[j][1];
                spectrum_acc[j] += fft1.out[j][0] * fft1.out[j][0] + fft1.out[j][1] * fft1.out[j][1];
            }
        }
#endif
        float max_amp = 0.0f;
        for(int i = 0; i < FFT_size; ++i)
        {
            fft[i][0] /= scale;
            fft[i][1] /= scale;
            phase_fft[i] = atan2f(fft[i][1], fft[i][0]); 
            amplitude[i] = sqrt(fft[i][0] * fft[i][0] + fft[i][1] * fft[i][1]);
            if(amplitude[i] > max_amp)
                max_amp = amplitude[i];

            spectrum_acc[i] /= (float)(blocks + 1);
            int shifted = (i + FFT_size / 2) % FFT_size;
            spectrum_db[shifted] = 10.0 * log10((double)spectrum_acc[i] + 1e-20);
            phase_deg[shifted] = phase_fft[i] * 180.0 / M_PI;
        }    

        float amp_threshold = max_amp * .3f;
        corr_vec.clear();
        float max_corr = -2.0f, corr_angle = 0.0f;
        int j = 0;
        for(float ang = -90.0f; ang <= 90.0f; ang += 0.3f)
        {
            float correlation = 0.0f;
            float weight_sum = 0.0f;
            int bins_count = 0;
            for(int i = 1; i < FFT_size; ++i)
            {
                float freq_bin = TO_FLOAT(i) * bin;
                if(i > FFT_size / 2)
                    freq_bin -= SAMPLE_RATE;
                if(amplitude[i]< amp_threshold)
                    continue;
                float phase_model = ((2.0f * M_PI * (FREQ + freq_bin) * angle_th[j]) / speed_light);
                correlation += std::cos(phase_fft[i] - (phase_model - 2.443f)) * amplitude[i];
                weight_sum += amplitude[i];
            }
            if(weight_sum > 1e-12f)
               correlation /= weight_sum;
            else correlation = 0;
            corr_vec.push_back(correlation);
            if(max_corr < correlation)
            {
                max_corr = correlation;
                corr_angle = ang;
            }
            ++j;
        }

        if (max_corr < 0.65f)
        {
            peleng_line_1->x_data(std::vector<double>{0.0, 0.0});
            peleng_line_1->y_data(std::vector<double>{0.0, 0.0});
        }
        else
        {
            float angle_rad = corr_angle * M_PI / 180.0f;
            peleng_line_1->x_data(std::vector<double>{0.0, sinf(angle_rad)});
            peleng_line_1->y_data(std::vector<double>{0.0, cosf(angle_rad)});
        }
        fig_plr_ax->draw();

        plt->x_data(angle_vec);
        plt->y_data(corr_vec);
        plt->touch();
        fig_corr->draw();

        plt_spectrum->y_data(spectrum_db);
        plt_spectrum->touch();
        fig_spectrum->draw();

        // plt_phase->y_data(phase_deg);
        // plt_phase->touch();
        // fig_phase->draw();

        printf("angle %5.1f corr %5.2f\n", corr_angle, max_corr);
        usleep(200000);
#ifdef DEBUG
    }
#endif
    }

    if(rx1)
        free(rx1);
    if(rx2)
        free(rx2);
    sdr::free_config(&sdr);
    return EXIT_SUCCESS;
}