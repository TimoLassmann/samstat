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
        pd->id->str[0] = 'A';
        plot_add(out,pd);
        pd->id->str[0] = 'B';
        pd->mod =  PLOT_MOD_NORMAL;
        pd->viz = PLOT_VIZ_FIRSTGROUP;

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

        int group_size = 2 + tl_random_int(rng,6);
        int y = group_size * (tl_random_int(rng, 3)+1);


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
                if(i % group_size == 0){
                        group++;
                }
                snprintf(pd->series_label[i], 256, "%d_%c",i, tl_random_int(rng, 26) + 65);
                snprintf(pd->group_label[i], 256, "%d", group);
                fprintf(stderr,"Group label %d -> %d (groupsize = %d)\n", i, group,group_size);
                for(int j = 0; j < x;j++){
                        pd->data[i][j] =  tl_random_int(rng, 100);
                }
        }
        /* pd->type = (viz << PLOT_TYPE_SHIFT) | tl_random_int(rng,3); */

        pd->mod = PLOT_MOD_ERROR_BAR;
        pd->type = tl_random_int(rng,3);
        fprintf(stderr,"TPYE: %d",pd->type);

        double t = tl_random_double(rng);
        if(t < 0.5){
                fprintf(stderr,"Viz group!\n");
                pd->viz = PLOT_VIZ_FIRSTGROUP;
        }else{
                fprintf(stderr,"Viz all!\n");
                pd->viz = PLOT_VIZ_ALL;
        }
        pd->viz = PLOT_VIZ_ALL;
        pd->group_size =group_size;

        *p = pd;

        return OK;
ERROR:
        return FAIL;
}
