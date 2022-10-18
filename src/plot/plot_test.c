#include "tld.h"

#include "plot.h"

int random_plot_data(struct plot_data**p, struct rng_state* rng);

int main(void)
{
        struct plot_data* pd = NULL;
        struct rng_state* rng = NULL;
        tld_strbuf* out = NULL;

        tld_strbuf_alloc(&out, 1024);
        RUN(init_rng(&rng,0));
        RUN(random_plot_data(&pd,rng));

        html_header(out,"TEST");
        plot_add(out,pd);
        html_end(out);
        fprintf(stdout,"%s",TLD_STR(out));
        plot_data_free(pd);


        tld_strbuf_free(out);
        free_rng(rng);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


int random_plot_data(struct plot_data **p, struct rng_state *rng)
{
        struct plot_data *pd = NULL;
        char buf[17];
        int x = 10;
        int y = 4 + tl_random_int(rng, 12);
        int viz = 2 + tl_random_int(rng,2);

        plot_data_alloc(&pd, x,  y);

        for(int i = 0; i < 16;i++){
                buf[i] = tl_random_int(rng, 26) + 65;
        }
        buf[16] = 0;

        RUN(tld_append(pd->id, buf));
        RUN(tld_append(pd->title, "RandomPlot"));
        RUN(tld_append(pd->xlabel, "RandomXlabel"));
        RUN(tld_append(pd->ylabel, "RandomYlabel"));
        RUN(tld_append(pd->save_file_name, "Myplottosave"));
        int group =0;
        for(int i = 0; i < pd->L;i++){
                if(i % viz == 0){
                        group++;
                }
                snprintf(pd->series_label[i], 256, "%d_%c",i, tl_random_int(rng, 26) + 65);
                snprintf(pd->group_label[i], 256, "%d", group);
                fprintf(stderr,"Group label %d -> %d (viz = %d)\n", i, group,viz);
                for(int j = 0; j < x;j++){
                        pd->data[i][j] =  tl_random_int(rng, 100);
                }
        }
        pd->type = (viz << PLOT_TYPE_SHIFT) | tl_random_int(rng,3);



        *p = pd;

        return OK;
ERROR:
        return FAIL;
}
