#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <csignal>
#include "types.h"
#include "sdr.h"


#define SAMPLE_RATE 4096

COMPLEX rx1[SAMPLE_RATE];
COMPLEX rx2[SAMPLE_RATE];

sig_atomic_t doneman = 1;
void signal_handler(int)
{
    doneman = 0;
}

int main()
{
    CONFIG sdr;
    if(sdr::init_sdr(&sdr, "ip:192.168.2.1", 2412000000, SAMPLE_RATE) && sdr::free_config(&sdr))
        return EXIT_FAILURE;
    
    std::signal(SIGINT, signal_handler);
    while(doneman)
    {
        if(!sdr::sdr_receive(&sdr, rx1, rx2, SAMPLE_RATE))
            break;
        printf("%f %f | %f %f\n",
               rx1[0].i,
               rx1[0].q,
               rx2[0].i,
               rx2[0].q);
        
        usleep(100000);
    }

    sdr::free_config(&sdr);
    return EXIT_SUCCESS;
}