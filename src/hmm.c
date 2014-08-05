

#include "hmm.h"
#include <pthread.h>

static float logsum_lookup[LOGSUM_SIZE];
//static float logsub_lookup[LOGSUM_SIZE];

struct hmm* run_EM_iterations(struct hmm* hmm,struct hmm_data* data)
{
	int i;
	data->run_mode = MODE_BAUM_WELCH;
	
	for(i = 0; i < data->iterations;i++ ){
		fprintf(stderr,"ITER:%d\n", i);
		hmm = run_pHMM(hmm, data);
		hmm= reestimate_hmm_parameters(hmm);
	}
	return hmm;
}

struct hmm* run_pHMM(struct hmm* hmm,struct hmm_data* data)
{
	struct thread_data* thread_data = 0;
	int num_threads = data->num_threads;
	int numseq = data->num_seq;
	int mode = data->run_mode;
	
	MMALLOC(thread_data,sizeof(struct thread_data)* num_threads);
	pthread_t threads[num_threads];
	pthread_attr_t attr;
	int t;
	
	int interval = 0;
	int rc;
	
	/* to be safe */
	init_logsum();
	
	interval =  (int)((double)numseq /(double)num_threads);
	
	for(t = 0;t < num_threads ;t++) {
		thread_data[t].data = data;
		thread_data[t].hmm = copy_hmm(hmm);
		thread_data[t].start = t*interval;
		thread_data[t].end = t*interval + interval;
		
	}
	thread_data[num_threads-1].end = numseq;
	//unsigned int seed = (unsigned int) (time(NULL) * ( 42));
	
	rc = pthread_attr_init(&attr);
	if(rc){
		//sprintf(param->buffer,"ERROR; return code from pthread_attr_init() is %d\n", rc);
		//param->messages = append_message(param->messages, param->buffer);
		
		//free_param(param);
		exit(EXIT_FAILURE);
	}
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	for(t = 0;t < num_threads;t++) {
		switch (mode) {
			case MODE_BAUM_WELCH:
				rc = pthread_create(&threads[t], &attr, do_baum_welch, (void *) &thread_data[t]);
				break;
			case MODE_FORWARD:
				rc = pthread_create(&threads[t], &attr, do_forward, (void *) &thread_data[t]);
				break;
		}
		
		if (rc) {
			//sprintf(param->buffer,"ERROR; return code from pthread_create() is %d\n", rc);
			//param->messages = append_message(param->messages, param->buffer);
			//free_param(param);
			exit(EXIT_FAILURE );
		}
	}
	
	pthread_attr_destroy(&attr);
	
	for (t = 0;t < num_threads;t++){
		rc = pthread_join(threads[t], NULL);
		if (rc){
			//sprintf(param->buffer,"ERROR; return code from pthread_join()is %d\n", rc);
			//param->messages = append_message(param->messages, param->buffer);
			//free_param(param);
			exit(EXIT_FAILURE );
		}
	}
	
	for (t = 0;t < num_threads;t++){
		hmm =  copy_estimated_parameters(hmm, thread_data[t].hmm);
	}
	
	
	
	for(t = 0;t < num_threads;t++) {
		free_hmm(thread_data[t].hmm);
	}
	
	MFREE(thread_data);
	return hmm;
}

void* do_baum_welch(void *threadarg)
{
	struct thread_data *data;
	data = (struct thread_data *) threadarg;
	
	struct hmm_data* hmm_data =  data->data;
	struct hmm* hmm = data->hmm;
	
	int start = data->start;
	int end = data->end;
	
	int i;
	
	for(i = start; i < end;i++){
		hmm = forward(hmm,hmm_data->string[i], hmm_data->length[i]);
		hmm = backward(hmm,hmm_data->string[i], hmm_data->length[i]);
		hmm = collect_estimated(hmm,hmm_data->string[i], hmm_data->weight[i], hmm_data->length[i]);
	}
	pthread_exit((void *) 0);
}

void* do_forward(void *threadarg)
{
	struct thread_data *data;
	data = (struct thread_data *) threadarg;
	
	struct hmm_data* hmm_data =  data->data;
	struct hmm* hmm = data->hmm;
	
	int start = data->start;
	int end = data->end;
	
	int i;
	
	for(i = start; i < end;i++){
		hmm = forward(hmm,hmm_data->string[i], hmm_data->length[i]);
		hmm_data->score[i] = hmm->f_score;
		
		
	}
	pthread_exit((void *) 0);
	
}

struct hmm* malloc_hmm(int num_states, int alphabet_len, int max_seq_len)
{
	struct hmm* hmm = NULL;
	int i = 0;
	MMALLOC(hmm, sizeof(struct hmm));
	hmm->alphabet_len =  alphabet_len;
	hmm->num_states = num_states;
	hmm->max_seq_len = max_seq_len;
	hmm->emissions = NULL;
	hmm->transitions = NULL;
	
	hmm->tindex = NULL;
	
	hmm->emissions_e = NULL;
	hmm->transitions_e = NULL;
	
	hmm->F = NULL;
	hmm->B = NULL;
	
	hmm->F_memory = NULL;
	hmm->B_memory = NULL;
	
	MMALLOC(hmm->emissions,sizeof(float*) * num_states);
	MMALLOC(hmm->transitions,sizeof(float*) * num_states);
	MMALLOC(hmm->tindex,sizeof(float*) * num_states);
	MMALLOC(hmm->emissions_e,sizeof(float*) * num_states);
	MMALLOC(hmm->transitions_e,sizeof(float*) * num_states);
	
	
	
	for(i = 0; i < num_states;i++){
		hmm->emissions[i] = NULL;
		hmm->transitions[i] = NULL;
		hmm->emissions_e[i] = NULL;
		hmm->transitions_e[i] = NULL;
		hmm->tindex[i] = NULL;
		//hmm->F[i] = NULL;
		//hmm->B[i] = NULL;
		
		MMALLOC(hmm->transitions[i],sizeof(float) * num_states);
		MMALLOC(hmm->transitions_e[i],sizeof(float) * num_states);
		MMALLOC(hmm->tindex[i],sizeof(float*) * (num_states+1));
		
		//MMALLOC(hmm->F[i],sizeof(float) * max_seq_len);
		//MMALLOC(hmm->B[i],sizeof(float) * max_seq_len);
		
		
	}
	MMALLOC(hmm->F_memory ,sizeof(float*) * num_states *max_seq_len );
	MMALLOC(hmm->B_memory, sizeof(float*) * num_states * max_seq_len);
	
	MMALLOC(hmm->F ,sizeof(float*) * max_seq_len);
	MMALLOC(hmm->B,sizeof(float*) * max_seq_len);
	
	for(i = 0; i < max_seq_len;i++){
		hmm->F[i] = (float*) (hmm->F_memory + i *num_states);
		hmm->B[i] = (float*) (hmm->B_memory + i *num_states);
	}
	//MMALLOC(hmm->F ,sizeof(float*) * num_states);
	//MMALLOC(hmm->B,sizeof(float*) * num_states);
	
	
	return hmm;
}

void free_hmm(struct hmm* hmm)
{
	int i;
	for(i = 0; i < hmm->num_states;i++){
		if(hmm->emissions[i]){
			MFREE(hmm->emissions[i]);
			MFREE(hmm->emissions_e[i]);
		}
		MFREE(hmm->transitions[i]);
		MFREE(hmm->transitions_e[i]);
		MFREE(hmm->tindex[i]);
		
		//MFREE(hmm->F[i]);
		//MFREE(hmm->B[i]);
		
	}
	
	MFREE(hmm->F);
	MFREE(hmm->B);
	MFREE(hmm->F_memory);
	MFREE(hmm->B_memory);
	
	MFREE(hmm->transitions);
	MFREE(hmm->tindex);
	MFREE(hmm->transitions_e);
	MFREE(hmm->emissions);
	MFREE(hmm->emissions_e);
	MFREE(hmm);
}


struct hmm* init_hmm_simple_ID(struct hmm* hmm)
{
	int i = 0;
	int j = 0;
	//transitions
	for(i = 0; i < hmm->num_states;i++){
		for(j = 0; j < hmm->num_states;j++){
			hmm->transitions[i][j] = -INFINITY;// prob2scaledprob(0.0);
			hmm->transitions_e[i][j] = prob2scaledprob(0.5f);
			
		}
	}
	// first state and seconf... (doubles as start and end
	for(j = 0; j < hmm->num_states;j++){
		if(j!= 0){
			hmm->transitions[STARTSTATE][j] =prob2scaledprob (1.0f / (float) (hmm->num_states-1));
			
		}
		if(j != 1){
			hmm->transitions[2][j] =prob2scaledprob (1.0f / (float) (hmm->num_states-1));
		}
		hmm->transitions[ENDSTATE][j] =prob2scaledprob (0.0f);
		
	}
	// remaining states..
	
	for(i = 3;i < hmm->num_states;i++){
		if(i == hmm->num_states-1){
			hmm->transitions[i][ENDSTATE] = prob2scaledprob(0.5f);
			hmm->transitions[i][2] = prob2scaledprob(0.5f);
		}else{
			hmm->transitions[i][ENDSTATE] = prob2scaledprob(0.05f);
			hmm->transitions[i][2] = prob2scaledprob(0.05f);
			hmm->transitions[i][i+1] = prob2scaledprob(0.9f);
		}
	}
	
	
	//emission (remember to malloc!
	// state 0 is delete state - all others emit!
	for(i =1 ; i < hmm->num_states;i++){
		MMALLOC(hmm->emissions[i], sizeof(float) * hmm->alphabet_len);
		MMALLOC(hmm->emissions_e[i], sizeof(float) * hmm->alphabet_len);
		for(j = 0;j < hmm->alphabet_len;j++){
			hmm->emissions[i][j] = prob2scaledprob(1.0f / (float) hmm->alphabet_len);
			hmm->emissions_e[i][j] = prob2scaledprob(0.5f);
		}
	}
	return hmm;
}

void print_hmm_parameters(struct hmm* hmm)
{
	int i,j;
	fprintf(stderr,"\nEmission:\n");
	
	
	for(i = 0; i < hmm->num_states ;i++){
		fprintf(stderr,"%d:",i);
		if(hmm->emissions[i]){
			for(j = 0; j < hmm->alphabet_len;j++){
				fprintf(stderr," %f", scaledprob2prob(hmm->emissions[i][j]));
			}
			fprintf(stderr," \n");
		}else{
			fprintf(stderr," NULL\n");
		}
	}
	fprintf(stderr,"\nTransition:\n");
	
	for(i = 0; i < hmm->num_states ;i++){
		//fprintf(stdout,"%d:",i);
		for(j = 1; j < hmm->tindex[i][0];j++){
			fprintf(stderr,"from %d to %d:%0.3f\n",i,hmm->tindex[i][j],scaledprob2prob(hmm->transitions[i][hmm->tindex[i][j]] ));
			
		}
	}
	
}

struct hmm* reestimate_hmm_parameters(struct hmm* hmm)
{
	float sum = 0.0f;
	int i,j;
	
	for(i = 2; i < hmm->num_states ;i++){
		//if(hmm->emissions[i]){
		sum = -INFINITY;// prob2scaledprob(0.0f);
		for(j = 0; j < hmm->alphabet_len;j++){
			sum = logsum(sum, hmm->emissions_e[i][j] );
		}
		for(j = 0; j < hmm->alphabet_len;j++){
			hmm->emissions[i][j] = hmm->emissions_e[i][j] -sum;
		}
		//}
	}
	
	for(i = 0; i < hmm->num_states ;i++){
		sum = -INFINITY;// prob2scaledprob(0.0f);
		for(j = 0; j < hmm->num_states;j++){
			if(hmm->transitions_e[i][j] != -INFINITY){
				sum = logsum(sum, hmm->transitions_e[i][j]);
			}
 		}
		
		for(j = 0; j < hmm->num_states;j++){
			if(hmm->transitions_e[i][j] != -INFINITY){
				hmm->transitions[i][j] = hmm->transitions_e[i][j] - sum;
			}else{
				hmm->transitions[i][j] = -INFINITY;// prob2scaledprob(0.0);
			}
 		}
	}
	
	for(i = 0; i < hmm->num_states;i++){
		for(j = 0; j < hmm->num_states;j++){
			//hmm->transitions[i][j] = prob2scaledprob(0.0);
			if(hmm->transitions[i][j] != -INFINITY){
				hmm->transitions_e[i][j] = prob2scaledprob(0.5f);
			}else{
				hmm->transitions_e[i][j] = -INFINITY;// prob2scaledprob(0.0);
			}
			
		}
	}
	
	
	//emission (remember to malloc!
	// state 0 is delete state - all others emit!
	for(i =2 ; i < hmm->num_states;i++){
		for(j = 0;j < hmm->alphabet_len;j++){
			//hmm->emissions[i][j] = prob2scaledprob(1.0f / (float) hmm->alphabet_len);
			hmm->emissions_e[i][j] = prob2scaledprob(0.5f);
		}
	}
	
	
	
	return hmm;
}


void print_hmm_estimated_parameters(struct hmm* hmm)
{
	int i,j;
	fprintf(stderr,"\nEmission:\n");
	
	
	for(i = 0; i < hmm->num_states ;i++){
		fprintf(stderr,"%d:",i);
		if(hmm->emissions[i]){
			for(j = 0; j < hmm->alphabet_len;j++){
				fprintf(stderr," %0.3f", scaledprob2prob(hmm->emissions_e[i][j]));
			}
			fprintf(stderr," \n");
		}else{
			fprintf(stderr," NULL\n");
		}
	}
	fprintf(stderr,"\nTransition:\n");
	for(i = 0; i < hmm->num_states ;i++){
		//fprintf(stdout,"%d:",i);
		for(j = 1; j < hmm->tindex[i][0];j++){
			fprintf(stderr,"from %d to %d:%0.3f\n",i,hmm->tindex[i][j],scaledprob2prob(hmm->transitions_e[i][hmm->tindex[i][j]] ));
			
		}
		//f = hmm->tindex[j][c];
		
		//for(j = 0; j < hmm->num_states;j++){
		//	fprintf(stdout," %0.3f",scaledprob2prob (hmm->transitions_e[i][j]));
		//}
		//fprintf(stdout," \n");
	}
}


void print_dyn_matrix(struct hmm* hmm, int seq_len)
{
	int i,j;
	seq_len = seq_len+2;
	
	if(seq_len < 10){
		fprintf(stdout,"FORWARD:\n");
		
		for(j = 0; j < seq_len ;j++){
			fprintf(stdout,"%d ",j);
			fprintf(stdout,"M:");
			for(i = 0; i < hmm->num_states ;i++){
				fprintf(stdout," %0.3f", scaledprob2prob(hmm->F[j][i]));
			}
			
			fprintf(stdout,"\n");
		}
		
		
	}
	fprintf(stdout,"BACKWARD:\n");
	//for(i = 0; i < hmm->num_states ;i++){
	for(j = 0; j < seq_len ;j++){
		fprintf(stdout,"%d ",j);
		fprintf(stdout,"M:");
		for(i = 0; i < hmm->num_states ;i++){
			fprintf(stdout," %0.3f", scaledprob2prob(hmm->B[j][i]));
		}
		
		fprintf(stdout,"\n");
	}
	//}
	
	/*float sum = 0;
	 fprintf(stdout,"SANITY:\n");
	 for(j = 0; j < seq_len ;j++){
	 sum = prob2scaledprob(0.0);
	 for(i = 0; i < hmm->num_states ;i++){
	 sum =logsum(sum,hmm->F[i][j] + hmm->B[i][j]  );
	 }
	 fprintf(stdout," %0.3f", (sum));
	 
	 }
	 fprintf(stdout,"\n");
	 }*/
	fprintf(stdout,"F:%f\n", (hmm->f_score));
	fprintf(stdout,"B:%f\n", (hmm->b_score));
	
	
	
}

struct hmm* copy_estimated_parameters(struct hmm* target,struct hmm* source )
{
	int i,j;
	
	
	for(i = 2; i < target->num_states ;i++){
		for(j = 0; j < target->alphabet_len;j++){
			target->emissions_e[i][j] = logsum(target->emissions_e[i][j], source->emissions_e[i][j]);
		}
	}
	
	for(i = 0; i < target->num_states ;i++){
		for(j = 0; j < target->num_states;j++){
			target->transitions_e[i][j] = logsum(target->transitions_e[i][j] , source->transitions_e[i][j] );
		}
	}
	return target;
}


struct hmm* copy_hmm(struct hmm* org )
{
	struct hmm* new = 0;
	int i,j;
	MMALLOC(new, sizeof(struct hmm));
	new->alphabet_len =  org->alphabet_len;
	new->num_states = org->num_states;
	new->max_seq_len = org->max_seq_len;
	new->emissions = NULL;
	new->transitions = NULL;
	
	new->emissions_e = NULL;
	new->transitions_e = NULL;
	
	new->F = NULL;
	new->B = NULL;
	new->F_memory = NULL;
	new->B_memory = NULL;
	
	new->tindex = NULL;
	
	MMALLOC(new->emissions,sizeof(float*) * org->num_states);
	MMALLOC(new->transitions,sizeof(float*) * org->num_states);
	MMALLOC(new->emissions_e,sizeof(float*) * org->num_states);
	MMALLOC(new->transitions_e,sizeof(float*) * org->num_states);
	
	MMALLOC(new->F_memory ,sizeof(float*) * org->num_states *org->max_seq_len );
	MMALLOC(new->B_memory, sizeof(float*) * org->num_states * org->max_seq_len);
	
	MMALLOC(new->F ,sizeof(float*) * org->max_seq_len);
	MMALLOC(new->B,sizeof(float*) * org->max_seq_len);
	
	for(i = 0; i < org->max_seq_len;i++){
		new->F[i] = (float*) (new->F_memory + i *org->num_states);
		new->B[i] = (float*) (new->B_memory + i *org->num_states);
	}
	
	
	
	MMALLOC(new->tindex,sizeof(float*) * org->num_states);
	
	
	for(i = 0; i < org->num_states;i++){
		new->emissions[i] = NULL;
		new->transitions[i] = NULL;
		new->emissions_e[i] = NULL;
		new->transitions_e[i] = NULL;
		new->tindex[i] = NULL;
		MMALLOC(new->tindex[i],sizeof(float*) * (org->num_states+1));
		MMALLOC(new->transitions[i],sizeof(float) * org->num_states);
		MMALLOC(new->transitions_e[i],sizeof(float) * org->num_states);
		
		for(j = 0; j < org->tindex[i][0];j++){
			new->tindex[i][j]= org->tindex[i][j];
		}
		
		
	}
	for(i =2 ; i < org->num_states;i++){
		MMALLOC(new->emissions[i], sizeof(float) * org->alphabet_len);
		MMALLOC(new->emissions_e[i], sizeof(float) * org->alphabet_len);
		for(j = 0; j < org->alphabet_len;j++){
			new->emissions[i][j] = org->emissions[i][j];
			new->emissions_e[i][j] = org->emissions_e[i][j];
		}
	}
	
	for(i = 0; i < org->num_states ;i++){
		for(j = 0; j < org->num_states;j++){
			new->transitions[i][j] = org->transitions[i][j];
			new->transitions_e[i][j] = org->transitions_e[i][j] ;
		}
	}
	
	
	
	return new;
}





struct hmm* forward(struct hmm* hmm, char* a, int len)
{
	int i,j,c,f;
	
	float** matrix = hmm->F;
	float* last= 0;
	float* cur = 0;
	const float* trans = 0;
	
	float tmp = 0;
	
	cur = matrix[0];
	
	//for(i = 0; i < hmm->num_states;i++){
	for(j = 0; j < hmm->num_states;j++){
		cur[j]  = -INFINITY;
	}
	cur[STARTSTATE] = 0.0f;
	
	for(i = 1; i < len+1;i++){
		last = cur;
		cur = matrix[i];
		for(j = 0; j < hmm->num_states;j++){
			cur[j] = -INFINITY;
		}
		
		for(j = 0; j < hmm->num_states;j++){
			tmp = last[j];
			trans = hmm->transitions[j];
			for(c = 1; c < hmm->tindex[j][0];c++){
				f = hmm->tindex[j][c];
				cur[f] = logsum(cur[f], tmp + trans[f] );//+ hmm->emissions[c][(int)a[i-1]]);
			}
			
		}
		for(c = 2;c < hmm->num_states;c++){
			cur[c] += hmm->emissions[c][(int)a[i-1]];
		}
	}
	
	//All goes to 1.
	last = cur;//matrix[len];
	cur = matrix[len+1];
	
	
	for(j = 0; j < hmm->num_states;j++){
		cur[j] = -INFINITY;// prob2scaledprob(0.0);
	}
	
	
	for(j = 2; j < hmm->num_states;j++){
		cur[ENDSTATE] = logsum(cur[ENDSTATE],last[j] + hmm->transitions[j][ENDSTATE]);
	}
	hmm->f_score = cur[ENDSTATE];// matrix[ENDSTATE][i];
	return hmm;
}


struct hmm* backward(struct hmm* hmm, char* a, int len)
{
	int i,j,c,f;
	
	float** matrix = hmm->B;
	
	float* next= 0;
	float* cur = 0;
	const float* trans = 0;
	
	cur = matrix[len+1];
	
	for(j = 0; j < hmm->num_states;j++){
		cur[j] = -INFINITY;
	}
	
	cur[ENDSTATE] = 0.0f;
	
	next = cur;
	
	cur = matrix[len];
	for(j = 0; j < hmm->num_states;j++){
		cur[j] =  hmm->transitions[j][ENDSTATE] + next[ENDSTATE];
	}
	for(c = 2;c < hmm->num_states;c++){
		cur[c] += hmm->emissions[c][(int)a[len-1]];
	}
	
	// backward recursion...
	for(i = len-1; i > 0; i -- ){
		next = cur;
		cur = matrix[i];
		for(j = 0; j < hmm->num_states;j++){
			trans = hmm->transitions[j];
			cur[j] = -INFINITY;
			for(c = 1; c < hmm->tindex[j][0];c++){
				f = hmm->tindex[j][c];
				cur[j] = logsum(cur[j],trans[f] + next[f]);// hmm->emissions[c][(int)a[i]]);// + next[c]);
			}
		}
		for(j = 2; j < hmm->num_states;j++){
			cur[j] += hmm->emissions[j][(int)a[i-1]];
		}
	}
	
	cur = matrix[0];
	next = matrix[1];
	
	for(j = 0; j < hmm->num_states;j++){
		cur[j] = -INFINITY;// prob2scaledprob(0.0f);
	}
	for(i = 0; i < hmm->num_states;i++){
		cur[0] = logsum(cur[0], hmm->transitions[0][i] + next[i]);//  + hmm->emissions[i][(int)a[0]]  );
	}
	hmm->b_score = cur[0];
	return hmm;
}

struct hmm* collect_estimated(struct hmm* hmm, char* a,float weight, int len)
{
	int i,j,c,f;
	float** Fmatrix = hmm->F;
	float** Bmatrix = hmm->B;
	
	float* last_F = 0;
	float* this_F = 0;
	float* this_B = 0;
	
	const float* trans = 0;
	
	float total = hmm->f_score + weight;

	for(i = 1; i < len+1;i++){
		last_F = Fmatrix[i-1];
		this_F = Fmatrix[i];
		this_B = Bmatrix[i];
		for(j = 0; j < hmm->num_states;j++){
			trans =hmm->transitions[j];
			for(c = 1; c < hmm->tindex[j][0];c++){
				f = hmm->tindex[j][c];
				hmm->transitions_e[j][f] = logsum(hmm->transitions_e[j][f], last_F[j] +trans[f] + this_B[f]- total);// +hmm->emissions[c][(int)a[i-1]]   -total);
			}
			if(j > 1){
				hmm->emissions_e[j][(int)a[i-1]] = logsum(hmm->emissions_e[j][(int)a[i-1]] ,(this_F[j]  +( this_B[j] -hmm->emissions[j][(int)a[i-1]] ))   -total  );
			}
		}
		
	}
	
	//All goes to 1.
	i = len+1;
	last_F = Fmatrix[len];
	//this_F = Fmatrix[i];
	this_B = Bmatrix[len+1];
	
	for(j = 0; j < hmm->num_states;j++){
		for(c = 0; c < hmm->num_states;c++){
			hmm->transitions_e[j][c] = logsum(hmm->transitions_e[j][c], last_F[j] +hmm->transitions[j][c] + this_B[c] -total);
			//hmm->transitions_e[c][j] = logsum(hmm->transitions_e[c][j], Fmatrix[c][i-1] + hmm->transitions[c][j] + Bmatrix[j][i]  - total);
		}
	}
	return hmm;
}

void init_logsum()
{
	static int called = 0;
	int i;
	if(!called){
		called = 1;
		for(i = 0; i < LOGSUM_SIZE;i++){
			logsum_lookup[i] = log(1.0 +exp((double) -i / SCALE));
			//		logsub_lookup[i]  = log(1.0 - exp((double) -i / SCALE));
		}
	}
}

float logsum(const float a,const float b)
{
	register const float max = MAX(a, b);
	register const float min = MIN(a, b);
	
	if(min == -INFINITY){
		return max;
	}
	if( (max-min) >= 15.7){
		return max;
	}
	return  max+ logsum_lookup[(int)((max-min)*SCALE)];
}

/*
 float logsub(const float a,const float b)
 {
 if(b == -INFINITY ){
 return a;
 }
 if( (a-b) >= 15.7){
 return a;
 }
 return  logsub_lookup[(int)((a-b )*SCALE)] + a;
 }
 
 float logsub_error(float a, float b)
 {
 float approx = logsub(a,b);
 float exact  = log(exp(a) - exp(b));
 return (exp(approx) - exp(exact));
 }
 
 float logsum_error(float a, float b)
 {
 float approx = logsum(a,b);
 float exact  = log(exp(a) + exp(b));
 return (exp(approx) - exp(exact));
 }*/

float prob2scaledprob(float p)
{
	if(p == 0.0){
		return -INFINITY;
	}else{
		return  log(p);
	}
}


float scaledprob2prob(float p)
{
	if(p == -INFINITY){
		return 0.0;
	}else{
		return exp(p);
	}
}

#ifdef ITEST
int main (int argc,char * argv[])
{
	fprintf(stderr,"Running hmm sanity tests\n");
	int len = 10;
	clock_t start, end;
	double cpu_time_used;
	
	
	//struct hmm* hmm = malloc_hmm(len,4);
	
	int i;
	
	init_logsum();
	
	struct hmm* hmm = complicated_hmm();
	//print_hmm_parameters(hmm);
	
	
	char* test_seq = 0;
	
	MMALLOC(test_seq, sizeof(char)*(len+500) );
	
	len = 500;
	for(i =0; i < 250;i++){
		test_seq[i] =  0;
	}
	
	for(i =250; i < 500;i++){
		test_seq[i] =  3;
	}
	test_seq[500] = 0;
	start = clock();
	for(i = 0; i < 10000;i++){
		
		hmm = backward(hmm,test_seq,len);
		hmm = forward(hmm,test_seq,len);
		
		
		
		
		//print_dyn_matrix(hmm,len);
		
		//fprintf(stderr,"Forward:\t%f\nBackward:\t%f\n", hmm->f_score, hmm->b_score);
		hmm = collect_estimated(hmm,test_seq,1.0, len);
	}
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	fprintf(stderr,"Forward:\t%f\nBackward:\t%f\n", hmm->f_score, hmm->b_score);
	
	fprintf(stderr,"%f time\n", cpu_time_used);
	
	print_hmm_estimated_parameters(hmm);
	exit(EXIT_SUCCESS);
	
	//
	//print_max_posterior(hmm,test_seq,len);
	/*float a = prob2scaledprob(-0.4);
	 float b =  prob2scaledprob (0.1);
	 
	 
	 double average_sumerror = 0.0;
	 double average_suberror= 0.0;
	 
	 
	 for(i = 0; i < 1000;i++){
	 
	 a=(float)rand()/((float)RAND_MAX / 1.0);
	 b=(float)rand()/((float)RAND_MAX/ (a)) ;
	 a = prob2scaledprob(a);
	 b = prob2scaledprob(b);
	 fprintf(stdout,"%f	%f	%f	%f	%g	%g = %g\n", scaledprob2prob( a),scaledprob2prob( b),exp(a) + exp(b), exp(logsub(a,b)), logsub_error(a,b), logsub_error(a,b),logsum_error(a,b)+logsub_error(a,b));
	 average_sumerror +=logsum_error(a,b);
	 average_suberror +=logsub_error(a,b);
	 
	 }
	 
	 fprintf(stderr," %g	%g\n",average_sumerror / (double)i, average_suberror / (double)i);
	 */
	//fprintf (stderr," %f	%g\n",logsum_error(-0.4,0.5),logsum_error(-0.4,0.5));
	
	exit(EXIT_SUCCESS);
	hmm = reestimate_hmm_parameters(hmm);
	
	print_hmm_parameters(hmm);
	
	
	print_dyn_matrix(hmm,len);
	//exit(0);
	len = 2;
	for(i =0; i < len;i++){
		test_seq[i] =  i & 3;
	}
	test_seq[len] = 0;
	
	hmm = forward(hmm,test_seq,len);
	
	hmm = backward(hmm,test_seq,len);
	print_dyn_matrix(hmm,len);
	//fprintf(stderr,"Forward:\t%f\nBackward:\t%f\n", hmm->f_score, hmm->b_score);
	
	free_hmm(hmm);
	
	MFREE(test_seq);
}



struct hmm* complicated_hmm(void)
{
	struct hmm* hmm = malloc_hmm(14, 4,1000);
	int i,j,c;
	for(i = 0; i < hmm->num_states;i++){
		for(j = 0; j < hmm->num_states;j++){
			hmm->transitions[i][j] = -INFINITY;//prob2scaledprob(0.0);
			hmm->transitions_e[i][j] = prob2scaledprob(0.5f);
			
		}
	}
	
	hmm->transitions[0][2] = prob2scaledprob(0.5f);
	hmm->transitions[0][8] = prob2scaledprob(0.5f);
	
	hmm->transitions[2][3] = prob2scaledprob(1.0f / 3.0f);
	hmm->transitions[2][4] = prob2scaledprob(1.0f / 3.0f);
	hmm->transitions[2][6] = prob2scaledprob(1.0f / 3.0f);
	
	hmm->transitions[3][3] = prob2scaledprob(0.5f);
	hmm->transitions[3][4] = prob2scaledprob(0.5f);
	
	hmm->transitions[4][5] = prob2scaledprob(0.5f);
	hmm->transitions[4][6] = prob2scaledprob(0.5f);
	
	hmm->transitions[5][5] = prob2scaledprob(0.5f);
	hmm->transitions[5][6] = prob2scaledprob(0.5f);
	
	hmm->transitions[6][7] = prob2scaledprob(0.5f);
	hmm->transitions[6][8] = prob2scaledprob(0.5f);
	
	hmm->transitions[7][7] = prob2scaledprob(0.5f);
	hmm->transitions[7][8] = prob2scaledprob(0.5f);
	
	hmm->transitions[8][9] = prob2scaledprob(1.0f / 3.0f);
	hmm->transitions[8][10] = prob2scaledprob(1.0f / 3.0f);
	hmm->transitions[8][12] = prob2scaledprob(1.0f / 3.0f);
	
	
	
	hmm->transitions[9][9] = prob2scaledprob(0.5f);
	hmm->transitions[9][10] = prob2scaledprob(0.5f);
	
	
	
	hmm->transitions[10][11] = prob2scaledprob(0.5f);
	hmm->transitions[10][12] = prob2scaledprob(0.5f);
	
	
	hmm->transitions[11][11] = prob2scaledprob(0.5f);
	hmm->transitions[11][12] = prob2scaledprob(0.5f);
	
	hmm->transitions[12][13] = prob2scaledprob(0.5f);
	hmm->transitions[12][ENDSTATE] = prob2scaledprob(0.5f);
	
	hmm->transitions[13][13] = prob2scaledprob(0.5f);
	hmm->transitions[13][ENDSTATE] = prob2scaledprob(0.5f);
	
	
	
	
	for(i =1 ; i < hmm->num_states;i++){
		MMALLOC(hmm->emissions[i], sizeof(float) * hmm->alphabet_len);
		MMALLOC(hmm->emissions_e[i], sizeof(float) * hmm->alphabet_len);
		for(j = 0;j < hmm->alphabet_len;j++){
			hmm->emissions[i][j] = prob2scaledprob(1.0f / (float) hmm->alphabet_len);
			hmm->emissions_e[i][j] = prob2scaledprob(0.5f);
		}
	}
	
	
	i = 2;
	hmm->emissions[i][0] = prob2scaledprob(0.7f);
	hmm->emissions[i][1] = prob2scaledprob(0.1f);
	hmm->emissions[i][2] = prob2scaledprob(0.1f);
	hmm->emissions[i][3] = prob2scaledprob(0.1f);
	
	i = 4;
	hmm->emissions[i][0] = prob2scaledprob(0.7f);
	hmm->emissions[i][1] = prob2scaledprob(0.1f);
	hmm->emissions[i][2] = prob2scaledprob(0.1f);
	hmm->emissions[i][3] = prob2scaledprob(0.1f);
	
	
	i = 6;
	hmm->emissions[i][0] = prob2scaledprob(0.7f);
	hmm->emissions[i][1] = prob2scaledprob(0.1f);
	hmm->emissions[i][2] = prob2scaledprob(0.1f);
	hmm->emissions[i][3] = prob2scaledprob(0.1f);
	
	i = 8;
	hmm->emissions[i][0] = prob2scaledprob(0.1f);
	hmm->emissions[i][1] = prob2scaledprob(0.1f);
	hmm->emissions[i][2] = prob2scaledprob(0.1f);
	hmm->emissions[i][3] = prob2scaledprob(0.7f);
	
	i = 10;
	hmm->emissions[i][0] = prob2scaledprob(0.1f);
	hmm->emissions[i][1] = prob2scaledprob(0.1f);
	hmm->emissions[i][2] = prob2scaledprob(0.1f);
	hmm->emissions[i][3] = prob2scaledprob(0.7f);
	
	
	i = 12;
	hmm->emissions[i][0] = prob2scaledprob(0.1f);
	hmm->emissions[i][1] = prob2scaledprob(0.1f);
	hmm->emissions[i][2] = prob2scaledprob(0.1f);
	hmm->emissions[i][3] = prob2scaledprob(0.7f);
	
	for(i = 0; i < hmm->num_states;i++){
		c = 0;
		for(j = 0; j < hmm->num_states;j++){
			if(hmm->transitions[i][j]  != -INFINITY){
				hmm->tindex[i][c+1] = j;
				c++;
			}
		}
		hmm->tindex[i][0] = c+1;
	}
	
	return hmm;
}


void print_max_posterior(struct hmm* hmm, char* a, int len)
{
	int i,j,c;
	float** Fmatrix = hmm->F;
	float** Bmatrix = hmm->B;
	
	float total = hmm->f_score;
	float max = 0.0f;
	c = -1;
	for(i = 1; i < len+1;i++){
		max = -INFINITY;// prob2scaledprob(0.0f);
		c = -1;
		for(j = 0; j < hmm->num_states;j++){
			//add emission
			if(Fmatrix[i][j]  + Bmatrix[i][j]-total  > max){
				max =Fmatrix[i][j]  + Bmatrix[i][j]-total  ;
				c = j;
			}
		}
		fprintf(stdout,"seq:%d State:%d	%f\n",a[i-1] , c, scaledprob2prob(max));
	}
}

#endif

