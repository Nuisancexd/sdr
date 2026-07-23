#include "sdr.h"
#include "types.h"

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
    iio_channel_attr_write_longlong(sdr->rx1_cfg, "hardwaregain", 50);
    iio_channel_attr_write(sdr->rx1_cfg, "rf_port_select", "A_BALANCED");
    iio_channel_attr_write_longlong(sdr->rx1_cfg, "sampling_frequency", sample_rate);

    EXECUTE_OR_GOTO(end, "[device_find_channel] failed", sdr->rx2_cfg = iio_device_find_channel(sdr->dev, "voltage1", false));
    iio_channel_attr_write(sdr->rx2_cfg, "gain_control_mode", "manual");
    iio_channel_attr_write_longlong(sdr->rx2_cfg, "hardwaregain", 50);
    iio_channel_attr_write(sdr->rx2_cfg, "rf_port_select", "A_BALANCED");
    iio_channel_attr_write_longlong(sdr->rx2_cfg, "sampling_frequency", sample_rate);
    
    //iio_channel_attr_write_longlong(sdr->rx1_cfg, "rf_bandwidth", 2000000);
    //iio_channel_attr_write_longlong(sdr->rx2_cfg, "rf_bandwidth", 2000000);

    EXECUTE_OR_GOTO(end, "[ctx_find_deviec] failed", sdr->rx = iio_context_find_device(sdr->ctx, "cf-ad9361-lpc"));
    
    sdr->rx1_i = iio_device_find_channel(sdr->rx, "voltage0", false);
    sdr->rx1_q = iio_device_find_channel(sdr->rx, "voltage1", false);
    sdr->rx2_i = iio_device_find_channel(sdr->rx, "voltage2", false);
    sdr->rx2_q = iio_device_find_channel(sdr->rx, "voltage3", false);

    iio_channel_enable(sdr->rx1_i);
    iio_channel_enable(sdr->rx1_q);
    iio_channel_enable(sdr->rx2_i);
    iio_channel_enable(sdr->rx2_q);

    EXECUTE_OR_GOTO(end, "[create_buff] failed", sdr->rxbuf = iio_device_create_buffer(sdr->rx, 1 << 16, false));

    iio_device_attr_write_bool(sdr->dev, "quadrature_tracking_en", true);
    iio_device_attr_write_bool(sdr->dev, "rf_dc_offset_tracking_en", true);
    iio_device_attr_write_bool(sdr->dev, "bb_dc_offset_tracking_en", true);
    iio_device_attr_write(sdr->dev, "calib_mode", "auto");
    usleep(200000);

    sdr->buf_pos = NULL;
    sdr->buf_end = NULL;

end:
    return return_code;
}

void sdr::change_channel_freq(PCONFIG sdr, size_t freq)
{
    iio_channel_attr_write_longlong(sdr->chn, "frequency", freq);
}

bool sdr::create_buffer(PCONFIG sdr, struct iio_buffer** buff_rx_out, size_t samples_count)
{
    if(!(*buff_rx_out = iio_device_create_buffer(sdr->rx, samples_count, false)))
    {
        printf("[create_buff] failed\n");
        return false;
    }
    return true;
}

void sdr::free_buf_rx(struct iio_buffer** buff_rx_out)
{
    if(*buff_rx_out)
        iio_buffer_destroy(*buff_rx_out);
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
    ptrdiff_t step = iio_buffer_step(sdr->rxbuf);
    while(received < samples)
    {
        if(sdr->buf_pos >= sdr->buf_end)
        {
            ssize_t bytes = iio_buffer_refill(sdr->rxbuf);
            if(bytes < 0) return false;

            sdr->buf_pos = (char*)iio_buffer_first(sdr->rxbuf, sdr->rx1_i);
            sdr->buf_end = (char*)iio_buffer_end(sdr->rxbuf);
        }

        while (sdr->buf_pos < sdr->buf_end && received < samples)
        {
            int16_t* iq = (int16_t*)sdr->buf_pos;
            rx1[received].i = (float)iq[0];
            rx1[received].q = (float)iq[1];
            rx2[received].i = (float)iq[2];
            rx2[received].q = (float)iq[3];
            sdr->buf_pos += step;
            ++received;
        }
    }

    return received != 0 ? true : false;
}

bool sdr::sdr_receive_one_channel(PCONFIG sdr, struct iio_buffer* buff_rx, PCOMPLEX rx, size_t samples)
{
    size_t received = 0;
    ptrdiff_t step = iio_buffer_step(buff_rx);

    while(received < samples)
    {
        if(sdr->buf_pos >= sdr->buf_end)
        {
            ssize_t bytes = iio_buffer_refill(buff_rx);
            if(bytes < 0) return false;

            sdr->buf_pos = (char*)iio_buffer_first(buff_rx, sdr->rx1_i);
            sdr->buf_end = (char*)iio_buffer_end(buff_rx);
        }

        while (sdr->buf_pos < sdr->buf_end && received < samples)
        {
            int16_t* iq = (int16_t*)sdr->buf_pos;
            rx[received].i = (float)iq[0];
            rx[received].q = (float)iq[1];
            sdr->buf_pos += step;
            ++received;
        }
    }

    return received != 0 ? true : false;
}