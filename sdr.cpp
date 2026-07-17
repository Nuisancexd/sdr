#include "sdr.h"

#include <iio.h>
#include <cstdio>
#include <unistd.h>

#define EXECUTE_OR_GOTO(label, msg, ...) \
    do                                   \
    {                                    \
        if (!(__VA_ARGS__))              \
        {                                \
            printf(msg "\n");            \
            return_code = EXIT_FAILURE;  \
            goto label;                  \
        }                                \
    } while(0)


bool sdr::init_sdr(PCONFIG sdr, const char* uri, size_t freq, size_t sample_rate)
{
    int return_code = EXIT_SUCCESS;
    EXECUTE_OR_GOTO(end, "[context] Failed connect to SDR", sdr->ctx = iio_create_context_from_uri(uri));
    EXECUTE_OR_GOTO(end, "[device] Failed connect to SDR", sdr->dev = iio_context_find_device(sdr->ctx, "ad9361-phy"));
    EXECUTE_OR_GOTO(end, "[device_find_ch] Failed find channel", sdr->chn = iio_device_find_channel(sdr->dev, "altvoltage0", true));
    iio_channel_attr_write_longlong(sdr->chn, "frequency", freq);

    EXECUTE_OR_GOTO(end, "[device_find_channel] failed", sdr->rx1_cfg = iio_device_find_channel(sdr->dev, "voltage0", false));
    iio_channel_attr_write(sdr->rx1_cfg, "gain_control_mode", "manual");
    iio_channel_attr_write_longlong(sdr->rx1_cfg, "hardwaregain", 10);
    iio_channel_attr_write(sdr->rx1_cfg, "rf_port_select", "A_BALANCED");
    iio_channel_attr_write_longlong(sdr->rx1_cfg, "sampling_frequency", sample_rate);

    EXECUTE_OR_GOTO(end, "[device_find_channel] failed", sdr->rx2_cfg = iio_device_find_channel(sdr->dev, "voltage1", false));
    iio_channel_attr_write(sdr->rx2_cfg, "gain_control_mode", "manual");
    iio_channel_attr_write_longlong(sdr->rx2_cfg, "hardwaregain", 10);
    iio_channel_attr_write(sdr->rx2_cfg, "rf_port_select", "A_BALANCED");
    iio_channel_attr_write_longlong(sdr->rx2_cfg, "sampling_frequency", sample_rate);

    EXECUTE_OR_GOTO(end, "[ctx_find_deviec] failed", sdr->rx = iio_context_find_device(sdr->ctx, "cf-ad9361-lpc"));
    
    sdr->rx1_i = iio_device_find_channel(sdr->rx, "voltage0", false);
    sdr->rx1_q = iio_device_find_channel(sdr->rx, "voltage1", false);
    sdr->rx2_i = iio_device_find_channel(sdr->rx, "voltage2", false);
    sdr->rx2_q = iio_device_find_channel(sdr->rx, "voltage3", false);

    iio_channel_enable(sdr->rx1_i);
    iio_channel_enable(sdr->rx1_q);
    iio_channel_enable(sdr->rx2_i);
    iio_channel_enable(sdr->rx2_q);

    EXECUTE_OR_GOTO(end, "[create_buff] failed", sdr->rxbuf = iio_device_create_buffer(sdr->rx, 4096, false));

end:
    return return_code;
}

bool sdr::free_config(PCONFIG sdr)
{
    if(sdr->rxbuf)
        iio_buffer_destroy(sdr->rxbuf);
    if(sdr->ctx)
        iio_context_destroy(sdr->ctx);
    return true;
}

bool sdr::sdr_receive(PCONFIG sdr, PCOMPLEX rx1, PCOMPLEX rx2, size_t samples)
{
    size_t received = 0;
    while(received < samples)
    {
        ssize_t bytes = iio_buffer_refill(sdr->rxbuf);
        if(bytes < 0) return false;

        ptrdiff_t step = iio_buffer_step(sdr->rxbuf);
        char* first = (char*)iio_buffer_first(sdr->rxbuf, sdr->rx1_i);
        char* end = (char*)iio_buffer_end(sdr->rxbuf);

        while (first < end && received < samples)
        {
            int16_t* iq = (int16_t*)first;
            rx1[received].i = (float)iq[0];
            rx1[received].q = (float)iq[1];
            rx2[received].i = (float)iq[2];
            rx2[received].q = (float)iq[3];
            first += step;
            ++received;
        }
    }

    return received != 0 ? true : false;
}
