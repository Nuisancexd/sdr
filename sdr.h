#ifndef _SDR_H_
#define _SDR_H_

#include "types.h"
#include <cstdint> 
#include <cstdlib>

typedef struct iio_config
    {
        struct iio_context* ctx;
        struct iio_channel* rx1_i;
        struct iio_channel* rx1_q;
        struct iio_channel* rx2_i;
        struct iio_channel* rx2_q;
        struct iio_buffer* rxbuf;
        struct iio_device* dev;
        struct iio_channel* chn;
        struct iio_channel* rx1_cfg;
        struct iio_channel* rx2_cfg;
        struct iio_device* rx;
    }CONFIG, *PCONFIG;

namespace sdr
{
    bool init_sdr(PCONFIG sdr, const char* uri, size_t freq, size_t sample_rate);
    bool free_config(PCONFIG sdr);
    bool sdr_receive(PCONFIG sdr, PCOMPLEX rx1, PCOMPLEX rx2, size_t samples);
    bool receive_block(PCONFIG sdr, PCOMPLEX rx1, PCOMPLEX rx2, size_t count);
    void change_channel_freq(PCONFIG sdr, size_t freq);
    bool sdr_receive_one_channel(PCONFIG sdr, PCOMPLEX rx, size_t samples);
}

#endif