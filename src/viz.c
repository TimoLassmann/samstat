
#include "viz.h"


#ifdef ITEST
int main (int argc,char * argv[])
{
	fprintf(stderr,"Running hmm sanity tests\n");
	int i;
	
	struct plot_data* pd = 0;
	pd = malloc_plot_data(4, 32);
	
	for(i = 0; i < 32;i++){
		pd->labels[i][0] = (char) (55+i);
		pd->data[0][i] = (float)(i +2) / 10.0;
		pd->data[1][i] = (float)(10 - i) / 10.0-5;
		
		pd->data[2][i] = (float)i*(float)i/ 10;
		pd->data[3][i] = (float)(i - 1) + 2;
	}
	pd->series_labels[0][0] = 'A';
	pd->series_labels[1][0] = 'C';
	pd->series_labels[2][0] = 'G';
	pd->series_labels[3][0] = 'T';

	
	
	print_html5_header(stdout,"TEST");
	pd->plot_type = LINE_PLOT;
	
	
	print_html5_chart(stdout, pd);
	print_html_table(stdout, pd);
	pd->plot_type = BAR_PLOT;
	print_html5_chart(stdout, pd);
	pd->plot_type = PIE_PLOT;
	print_html5_chart(stdout, pd);
	pd->plot_type = RADAR_PLOT;
	print_html5_chart(stdout, pd);
	print_html5_footer(stdout);
	
	free_plot_data(pd);
	
	exit(EXIT_SUCCESS);
}

#endif


struct plot_data* malloc_plot_data(int num_series,int num_points)
{
	struct plot_data* pd = NULL;
	int i,j;
	MMALLOC(pd, sizeof(struct plot_data));
	pd->num_series = num_series;
	pd->num_points = num_points;
	pd->labels = NULL;
	pd->series_labels = NULL;
	pd->data = NULL;
	pd->description = NULL;
	pd->plot_title = NULL;
	
	pd->plot_type = 0;
	pd->width = 700;
	pd->height = 300;
	
	MMALLOC(pd->description, sizeof(char) * 100000);
	MMALLOC(pd->plot_title, sizeof(char) * 100000);
	
	MMALLOC(pd->labels , sizeof(char* ) * num_points);
	
	MMALLOC(pd->series_labels , sizeof(char* ) * num_series);
	
	MMALLOC(pd->data , sizeof(float* ) * num_series);
	for(i = 0; i < num_points;i++){
		pd->labels[i] = 0;
		MMALLOC(pd->labels[i], sizeof(char)* MAXLABEL_LEN);
		for(j = 0; j < MAXLABEL_LEN;j++){
			pd->labels[i][j] = 0;
		}
		
	}
	for(i = 0; i < num_series;i++){
		pd->data[i] = 0;
		pd->series_labels[i] = 0;
		MMALLOC(pd->series_labels[i], sizeof(char)* MAXLABEL_LEN);
		for(j = 0; j < MAXLABEL_LEN;j++){
			pd->series_labels[i][j] = 0;
		}
		MMALLOC(pd->data[i], sizeof(float)*  num_points);
		for(j = 0; j < num_points;j++){
			pd->data[i][j] = 0.0f;

		}
	}
	return pd;
}

void free_plot_data(struct plot_data* pd)
{
	int i;
	for(i = 0; i < pd->num_series;i++){
		MFREE(pd->data[i]);
	}
	for(i = 0; i < pd->num_points;i++){
		MFREE(pd->labels[i]);// = 0;
	}
	MFREE(pd->data);
	MFREE(pd->labels);
	MFREE(pd->description);
	MFREE(pd->plot_title);
	MFREE(pd);
	
}





void print_html5_header(FILE* out,struct plot_data* pd)
{
	int i;
	
	fprintf(out,"<!doctype html>\n");
	fprintf(out,"<html>\n");
	fprintf(out,"<head>\n");
	fprintf(out,"<title>%s</title>\n", pd->plot_title);
	
	fprintf(out,"<script>");
	for(i = 0; i < Chartjs_len;i++ ){
		fprintf(out,"%c",Chartjs[i]);
	}
	fprintf(out,"</script>");
	fprintf(out,"<meta name = \"viewport\" content = \"initial-scale = 1, user-scalable = no\">\n");
	fprintf(out,"<style>\n");
	fprintf(out,"canvas{\n");
	fprintf(out,"}\n");
	
	fprintf(out,"body {\n");
	fprintf(out,"margin: 0 auto;\n");
	fprintf(out,"padding: 22px 0;\n");
	fprintf(out,"width: 940px;\n");
	//fprintf(out,"font-family: \"Open Sans\";\n");
	fprintf(out,"font: 14px Arial,Helvetica, sans-serif;\n");
	fprintf(out,"background: #F0F0F0;\n");
	fprintf(out,"}\n");
	
	fprintf(out,"#intro {\n");
	fprintf(out,"position: relative;\n");
	fprintf(out,"	margin-top: 66px;\n");
	fprintf(out,"padding: 44px;\n");
	fprintf(out,"background: #3366AA;\n");
	
	
	
	
	/* Border-radius not implemented yet */
	fprintf(out,"	-moz-border-radius: 22px;\n");
	fprintf(out,"	-webkit-border-radius: 22px;\n");
	fprintf(out,"}\n");
	
	fprintf(out,"#intro h1, #intro p {\n");
	fprintf(out,"position: relative;\n");
	fprintf(out,"z-index: 9999;\n");
	fprintf(out,"width: 620px;\n");
	fprintf(out,"}\n");
	
	fprintf(out,"#intro h1 {\n");
	fprintf(out,"padding: 0 0 0 0;\n");
	fprintf(out,"font-weight: normal;\n");
	fprintf(out,"font-size: 300%%;\n");
	fprintf(out,"color: #fff;\n");
	fprintf(out,"}\n");
	
	fprintf(out,"#intro p {\n");
	fprintf(out,"padding: 0;\n");
	fprintf(out,"color:#CCC;\n");
	fprintf(out,"}\n");
	


	fprintf(out,"table.simple {background-color: #FFFFFF;border-collapse:collapse;text-align :right; border: 2px solid #000000; float:none;}\n");
	fprintf(out,"table.simple td  {border: 1px solid #000000;width:auto;padding:5px;font-family:\"Courier New\" Courier 	monospace;font-size: 12pt;}\n");
	fprintf(out,"table.simple th  {border: 2px solid #000000;width:auto;padding:5px;font-family:\"Courier New\" Courier 	monospace;font-size: 12pt;}\n");
	        
	fprintf(out,"footer {\n");
	fprintf(out,"position: absolute;\n");
	fprintf(out,"left: 0;\n");
	fprintf(out,"width: 100%%;\n");
	fprintf(out,"background: #222;\n");
	fprintf(out,"}\n");
	
	fprintf(out,"footer div {\n");
	fprintf(out,"display: table;\n");
	fprintf(out,"margin: 0 auto;\n");
	fprintf(out,"padding: 44px 0;\n");
	fprintf(out,"width: 940px;\n");
	fprintf(out,"color: #777;\n");
	fprintf(out,"}\n");
	
	fprintf(out,"footer div section {\n");
	fprintf(out,"display: table-cell;\n");
	fprintf(out,"width: 300px;\n");
	fprintf(out,"}\n");
	
	fprintf(out,"footer div #about, footer div #blogroll {\n");
	fprintf(out,"	padding-right: 20px;\n");
	fprintf(out,"}\n");
	
	fprintf(out,"footer h3 {\n");
	fprintf(out,"color: #FFF;\n");
	fprintf(out,"}\n");
	
	fprintf(out,"footer a {\n");
	fprintf(out,"color: #999;\n");
	fprintf(out,"}\n");
	
	fprintf(out,"footer a:hover {\n");
	fprintf(out,"color: #FFF;\n");
	fprintf(out,"	text-decoration: none;\n");
	fprintf(out,"}\n");
	
	fprintf(out,"footer ul {\n");
	fprintf(out,"margin: 0 0 0 40px;\n");
	fprintf(out,"	list-style: square;\n");
	fprintf(out,"color: #565656;\n");
	fprintf(out,"}\n");
	
	fprintf(out,"footer ul li a {\n");
	fprintf(out,"display: block;\n");
	fprintf(out,"}\n");
	
	
	
	fprintf(out,"nav {\n");
	//fprintf(out,"position: absolute;\n");
	fprintf(out,"left: 0;\n");
	fprintf(out,"width: 100%%;\n");
	fprintf(out,"background: #000000;\n");
	fprintf(out,"color: #FFF;\n");
	fprintf(out,"}\n");
	
	
	
	
	fprintf(out,"</style>\n");
	fprintf(out,"</head>\n");
	fprintf(out,"<body>\n");
	
	
	//fprintf(out,"<ul>\n");
	//fprintf(out,"<li  style=\"font-size: 300%%;\">%s</li>\n",pd->plot_title);
	//fprintf(out,"<li style=\"color: #999\" >%s</li>\n",pd->description);
	//fprintf(out,"<li><a href=\"#contact\">Contact</a></li>\n");
	//fprintf(out,"</ul>\n");
	
	
	fprintf(out,"<nav>\n");
	fprintf(out,"<table style=\"margin: 0 auto;width: 940px;\">\n");
	fprintf(out,"<tr>\n");
	fprintf(out,"<td style=\"font-size: 200%%;\">%s</td>\n",pd->plot_title);
	fprintf(out,"<td style=\"vertical-align: bottom; \">%s</td>\n",pd->description);
	fprintf(out,"</tr>\n");
	        fprintf(out,"</table>\n");
	fprintf(out,"</nav>\n");
	
	
	//fprintf(out,"<section id=\"intro\">\n");
	//fprintf(out,"<header>\n");
	//fprintf(out,"<h1>%s</h1>\n",pd->plot_title);
	//fprintf(out,"</header>\n");
	//fprintf(out,"<p>%s</p>\n",pd->description);
	
	//fprintf(out,"</section>\n");
	
}

void print_html5_footer(FILE* out)
{
	//fprintf(out,"http://www.chartjs.org");
	fprintf(out,"<footer>\n");
	fprintf(out,"<div>\n");
	fprintf(out,"<section id=\"about\">\n");
	fprintf(out,"<h3>About</h3>\n");
	fprintf(out,"<p>The vast amount of data produced by next-generation sequencing machines necesitates the development of efficient visualization tools. SAMStat addresses the basic need to display information about the obtained reads both quickly and in a concise form.</p> <p>The plots on this page were draw using the excellent <a href=\"http://www.chartjs.org\">Chartjs library</a>.</p>\n");
	fprintf(out,"</section>\n");
	fprintf(out,"<section id=\"contact\">\n");
	
	fprintf(out,"<h3>Contact</h3>\n");
	
	fprintf(out,"<p>SAMStat was developed by Timo Lassmann (timolassmann at gmail dot com).</p>\n");
	fprintf(out,"<h3>Please cite:</h3>\n");
	
	fprintf(out,"<p>Lassmann et al. (2010) \"SAMStat: monitoring biases in next generation sequencing data.\" Bioinformatics doi:10.1093/bioinformatics/btq614 [<a href =\"http://www.ncbi.nlm.nih.gov/pubmed/21088025/\">PMID: 21088025</a>] </p>\n");
	fprintf(out,"</section>\n");
	fprintf(out,"<section id=\"blogroll\">\n");
	fprintf(out,"<h3>Links</h3>\n");
	fprintf(out,"<ul>\n");
	fprintf(out,"<li><a href=\"http://www.yokohama.riken.jp/english/\">RIKEN</a></li>\n");
	fprintf(out,"<li><a href=\"http://www.osc.riken.jp/english/\">OMICS Science Center</a></li>\n");
	fprintf(out,"<li><a href=\"http://samtools.sourceforge.net/\">samtools</a></li>\n");
	fprintf(out,"<li><a href=\"http://code.google.com/p/bedtools/\">BEDtools</a></li>\n");
	fprintf(out,"<li><a href=\"http://msa.sbc.su.se/\">Kalign</a></li>\n");
	fprintf(out,"<li><a href=\"http://ngsview.sourceforge.net/\">NGSview</a></li>\n");
	
	
	fprintf(out,"</ul>\n");
	fprintf(out,"</section>\n");
	
	fprintf(out,"</div>\n");
	fprintf(out,"</footer>\n");
	
	
	
	fprintf(out,"</body>\n");
	fprintf(out,"</html>\n");
}

void print_html_table(FILE* out,struct plot_data* pd)
{
	
	int i,j;
	int start,stop;
	
	
	start = 0;
	
	stop = 15;
	
	if(pd->num_points < stop){
		stop = pd->num_points;
	}
	while(1){
		
		fprintf(out,"<table  class=\"simple\" >\n");
		fprintf(out,"<tr>\n");
		
		//head row..
		
		if(pd->series_labels){
			if(pd->num_points < 5){
				fprintf(out,"<td style=\"width:100px \"></td>\n");
			}else{
				fprintf(out,"<td style=\"width:50px \"></td>\n");
			}
		}
		for(i = start; i < stop;i++){
			if(pd->num_points < 5){
				fprintf(out,"<td style=\"width:90px\">%s</td>\n",pd->labels[i]);
			}else{
				fprintf(out,"<td style=\"width:45px\">%s</td>\n",pd->labels[i]);
			}
			
		}
		
		fprintf(out,"</tr>\n");
		for(j = 0;j < pd->num_series;j++){
			fprintf(out,"<tr>\n");
			if(pd->series_labels){
				fprintf(out,"<td style=\"background-color: %s;\" >%s</td>\n",colors[j & 3], pd->series_labels[j]);
			}
			for(i =start ; i < stop ;i++){
				fprintf(out,"<td>%0.1f</td>\n",pd->data[j][i]);
				
			}
			fprintf(out,"</tr>\n");
			
		}
		fprintf(out,"</table>\n");
		
		if(stop == pd->num_points){
			break;
		}
		fprintf(out,"<br>\n");
		start += 15;
		stop += 15;
		if(pd->num_points < stop){
			stop = pd->num_points;
		}
		
	}
}


void print_html5_chart(FILE* out,struct plot_data* pd)
{
	static int id = 1;
	int i,j;
	fprintf(out,"<section>\n");
	fprintf(out,"<h2>%s</h2>\n", pd->plot_title);
	int start,stop;
	int points_shown =15;
	
	start = 0;
	
	switch (pd->plot_type) {
		case RADAR_PLOT:
		case PIE_PLOT:
			stop =pd->num_points;
			points_shown =pd->num_points;
			break;
		default:
			stop = points_shown;
			
			if(pd->num_points < stop){
				stop = pd->num_points;
			}
			break;
	}
	
	while(1){
		fprintf(out,"<canvas id=\"canvas%d\" height=\"300\" width=\"%f\"></canvas>\n",id , 900.0 * (float)(stop-start )/ points_shown);
		fprintf(out,"<script>\n");
		
		fprintf(out,"var ChartData%d = {\n",id);
		fprintf(out,"labels : [");
		for(i =start ; i < stop;i++){
			if(i != start){
				fprintf(out,",");
			}
			fprintf(out,"\"%s\"",pd->labels[i]);
		}
		fprintf(out,"],\n");
		
		//labels : ["January","February","March","April","May","June","July"],
		fprintf(out,"datasets : [\n");
		
		for(j = 0;j < pd->num_series;j++){
			if(j){
				fprintf(out,",\n");
			}
			fprintf(out,"{\n");
			fprintf(out,"fillColor : \"%s\",\n",colors[j]);
			fprintf(out,"strokeColor : \"%s\",\n",colors[j]);
			fprintf(out,"pointColor : \"%s\",\n",colors[j]);
			fprintf(out,"data : [");
			for(i =start; i < stop;i++){
				if(i!= start){
					fprintf(out,",");
				}
				fprintf(out,"%f",pd->data[j][i]);
			}
			fprintf(out,"]\n");
			fprintf(out,"}");
			
		}
		fprintf(out,"\n]\n}\n");
		
		
		fprintf(out,"var pieData%d = [\n",id);
		for(i = 0; i < pd->num_series; i++){//start ; i < stop;i++){
			if(i!=0){
				fprintf(out,",");
			}
			fprintf(out,"{ value: %f, color:\"%s\"}",pd->data[i][0],colors[i & 3] );
		}
		fprintf(out,"];\n");
		
		
		
		fprintf(out,"var Chartopt = {datasetFill : false}\n");
		
		
		switch (pd->plot_type) {
			case LINE_PLOT:
				fprintf(out,"var myLine = new Chart(document.getElementById(\"canvas%d\").getContext(\"2d\")).Line(ChartData%d,Chartopt);\n",id,id);
				break;
			case BAR_PLOT:
				fprintf(out,"var myLine = new Chart(document.getElementById(\"canvas%d\").getContext(\"2d\")).Bar(ChartData%d);\n",id,id);
				break;
			case RADAR_PLOT:
				fprintf(out,"var myLine = new Chart(document.getElementById(\"canvas%d\").getContext(\"2d\")).Radar(ChartData%d);\n",id,id);
				break;
			case PIE_PLOT:
				fprintf(out,"var myLine = new Chart(document.getElementById(\"canvas%d\").getContext(\"2d\")).Pie(pieData%d);\n",id,id);
				break;
			default:
				fprintf(out,"var myLine = new Chart(document.getElementById(\"canvas%d\").getContext(\"2d\")).Bar(ChartData%d);\n",id,id);
				break;
		}
		
		
		
		
		fprintf(out,"</script>\n");
		if(stop == pd->num_points){
			break;
		}
		fprintf(out,"<br>\n");
		start += 15;
		stop += 15;
		if(pd->num_points < stop){
			stop = pd->num_points;
		}
		id++;
	}
	fprintf(out,"<br>");
	fprintf(out,"%s\n",  pd->description);
	id++;
}


