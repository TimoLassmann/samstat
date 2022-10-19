#include "tld.h"

#define PLOT_INTERVAL_IMPORT
#include "plot_interval.h"

int get_plot_interval(int len,int bin_start,int target_n_bin, int *interval)
{
        int base[3] = {2,5,10};
        int ival = 0;
        int igroup = 0;
        int m = 1;
        int run = 1;
        fprintf(stderr,"LEN: %d %d\n", len,bin_start);
        while(run){
                for(int i = 0; i < 3;i++){
                        ival = base[i] * m;
                        igroup = bin_start + ((len-bin_start) / ival);

                        if((len-bin_start) % ival != 0){
                                igroup += 1;
                        }
                        fprintf(stderr,"group count :  %d (int: %d) all : %d targer n_greoup: %d\n", igroup,ival, bin_start + (igroup - bin_start)*ival ,target_n_bin );
                        if(igroup <  target_n_bin){
                                run = 0;
                                break;
                        }
                }
                m *= 10;
                if (m == 10000000) {
                        ERROR_MSG("No interval found");
                }
        }
        fprintf(stderr,"val: %d\n", ival);
        *interval = ival;

        exit(0);
        return OK;
ERROR:
        return FAIL;
}

