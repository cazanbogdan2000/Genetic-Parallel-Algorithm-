#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "genetic_algorithm.h"

// this function merges two sorted vectors of individuals
individual* merge(individual* source, int start, int mid, int end, individual* destination) {
	int iA = start;
	int iB = mid;
	int i;

	for (i = start; i < end; i++) {
		if (end == iB || (iA < mid && cmpfunc(&source[iA], &source[iB]) <= 0)) {
			destination[i] = source[iA];
			iA++;
		} else {
			destination[i] = source[iB];
			iB++;
		}
	}
	return destination;
}


int read_input(sack_object **objects, int *object_count, int *sack_capacity, int *generations_count, int argc, char *argv[])
{
	FILE *fp;

	if (argc < 3) {
		fprintf(stderr, "Usage:\n\t./tema1 in_file generations_count\n");
		return 0;
	}

	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		return 0;
	}

	if (fscanf(fp, "%d %d", object_count, sack_capacity) < 2) {
		fclose(fp);
		return 0;
	}

	if (*object_count % 10) {
		fclose(fp);
		return 0;
	}

	sack_object *tmp_objects = (sack_object *) calloc(*object_count, sizeof(sack_object));

	for (int i = 0; i < *object_count; ++i) {
		if (fscanf(fp, "%d %d", &tmp_objects[i].profit, &tmp_objects[i].weight) < 2) {
			free(objects);
			fclose(fp);
			return 0;
		}
	}

	fclose(fp);

	*generations_count = (int) strtol(argv[2], NULL, 10);
	
	if (*generations_count == 0) {
		free(tmp_objects);

		return 0;
	}

	*objects = tmp_objects;

	return 1;
}

void print_objects(const sack_object *objects, int object_count)
{
	for (int i = 0; i < object_count; ++i) {
		printf("%d %d\n", objects[i].weight, objects[i].profit);
	}
}

void print_generation(const individual *generation, int limit)
{
	for (int i = 0; i < limit; ++i) {
		for (int j = 0; j < generation[i].chromosome_length; ++j) {
			printf("%d ", generation[i].chromosomes[j]);
		}

		printf("\n%d - %d\n", i, generation[i].fitness);
	}
}

void print_best_fitness(const individual *generation)
{
	printf("%d\n", generation[0].fitness);
}

void compute_fitness_function(const sack_object *objects, individual *generation, int object_count, int sack_capacity, int start)
{
	int weight;
	int profit;

	for (int i = start; i < object_count; ++i) {
		weight = 0;
		profit = 0;

		for (int j = 0; j < generation[i].chromosome_length; ++j) {
			if (generation[i].chromosomes[j]) {
				weight += objects[j].weight;
				profit += objects[j].profit;
			}
		}
		generation[i].fitness = (weight <= sack_capacity) ? profit : 0;
	}
}

int cmpfunc(const void *a, const void *b)
{
	int i;
	individual *first = (individual *) a;
	individual *second = (individual *) b;

	int res = second->fitness - first->fitness; // decreasing by fitness
	if (res == 0) {
		int first_count = 0, second_count = 0;

		for (i = 0; i < first->chromosome_length && i < second->chromosome_length; ++i) {
			first_count += first->chromosomes[i];
			second_count += second->chromosomes[i];
		}

		res = first_count - second_count; // increasing by number of objects in the sack
		if (res == 0) {
			return second->index - first->index;
		}
	}

	return res;
}

void mutate_bit_string_1(const individual *ind, int generation_index)
{
	int i, mutation_size;
	int step = 1 + generation_index % (ind->chromosome_length - 2);

	if (ind->index % 2 == 0) {
		// for even-indexed individuals, mutate the first 40% chromosomes by a given step
		mutation_size = ind->chromosome_length * 4 / 10;
		for (i = 0; i < mutation_size; i += step) {
			ind->chromosomes[i] = 1 - ind->chromosomes[i];
		}
	} else {
		// for even-indexed individuals, mutate the last 80% chromosomes by a given step
		mutation_size = ind->chromosome_length * 8 / 10;
		for (i = ind->chromosome_length - mutation_size; i < ind->chromosome_length; i += step) {
			ind->chromosomes[i] = 1 - ind->chromosomes[i];
		}
	}
}

void mutate_bit_string_2(const individual *ind, int generation_index)
{
	int step = 1 + generation_index % (ind->chromosome_length - 2);

	// mutate all chromosomes by a given step
	for (int i = 0; i < ind->chromosome_length; i += step) {
		ind->chromosomes[i] = 1 - ind->chromosomes[i];
	}
}

void crossover(individual *parent1, individual *child1, int generation_index)
{
	individual *parent2 = parent1 + 1;
	individual *child2 = child1 + 1;
	int count = 1 + generation_index % parent1->chromosome_length;

	memcpy(child1->chromosomes, parent1->chromosomes, count * sizeof(int));
	memcpy(child1->chromosomes + count, parent2->chromosomes + count, (parent1->chromosome_length - count) * sizeof(int));

	memcpy(child2->chromosomes, parent2->chromosomes, count * sizeof(int));
	memcpy(child2->chromosomes + count, parent1->chromosomes + count, (parent1->chromosome_length - count) * sizeof(int));
}

void copy_individual(const individual *from, const individual *to)
{
	memcpy(to->chromosomes, from->chromosomes, from->chromosome_length * sizeof(int));
}

void free_generation(individual *generation)
{
	int i;

	for (i = 0; i < generation->chromosome_length; ++i) {
		free(generation[i].chromosomes);
		generation[i].chromosomes = NULL;
		generation[i].fitness = 0;
	}
}

int min(int a, int b) {
	return a < b ? a : b;
}

void *thread_function(void *arg)
{
	Thread_arguments* t_args = (Thread_arguments*) arg;

		// set initial generation (composed of object_count individuals with a single item in the sack)
		// calculate the start and end points, in order to parallelize the input  
		int start = t_args->id * (double)t_args->object_count / t_args->P;
		int end = min((t_args->id + 1) * (double)t_args->object_count / t_args->P, t_args->object_count);
		for (int i = start; i < end; ++i) {
			t_args->current_generation[i].fitness = 0;
			t_args->current_generation[i].chromosomes = (int*) calloc(t_args->object_count, sizeof(int));
			t_args->current_generation[i].chromosomes[i] = 1;
			t_args->current_generation[i].index = i;
			t_args->current_generation[i].chromosome_length = t_args->object_count;

			t_args->next_generation[i].fitness = 0;
			t_args->next_generation[i].chromosomes = (int*) calloc(t_args->object_count, sizeof(int));
			t_args->next_generation[i].index = i;
			t_args->next_generation[i].chromosome_length = t_args->object_count;
		}
		// all threads wait to finish each other
		pthread_barrier_wait(t_args->barrier);

		// iterate for each generation
		for (int k = 0; k < t_args->generations_count; ++k) {
			t_args->cursor = 0;

			// compute fitness and sort by it
			// calculate the start and end points, in order to parallelize the
			// compute fitness function (which is also modified, in order to be
			// able to parallelize it) and also the qsort
			start = t_args->id * (double)t_args->object_count / t_args->P;
			end = min((t_args->id + 1) * (double)t_args->object_count / t_args->P, t_args->object_count);
			compute_fitness_function(t_args->objects, t_args->current_generation, end, t_args->sack_capacity, start);
			// keep first 30% children (elite children selection)
			// parralelized qsort
			qsort(t_args->current_generation + start, end - start, sizeof(individual), cmpfunc);
			// all threads wait to finish each other
			pthread_barrier_wait(t_args->barrier);

			// the following lines are used to merge all the sorted subvectors
			// we made only one merge (we can't paralellize merge)
			if (t_args->id == 0) {
				individual* aux = (individual*)calloc(t_args->object_count, sizeof(individual));
				for (int i = 1; i < t_args->P; i++) {
					int new_start = i * (double)t_args->object_count / t_args->P;
					int new_end = min((i + 1) * (double)t_args->object_count / t_args->P, t_args->object_count);
					aux = merge(t_args->current_generation, 0, new_start, new_end, aux);
					memcpy(t_args->current_generation, aux, new_end * sizeof(individual));
					memset(aux, 0, t_args->object_count * sizeof(individual));
				}				
			}
			// all threads wait to finish each other
			pthread_barrier_wait(t_args->barrier);

			// keep first 30% (elite selection) + parallel
			t_args->count = t_args->object_count * 3 / 10;
			start = t_args->id * (double)t_args->count / t_args->P;
			end = min((t_args->id + 1) * (double)t_args->count / t_args->P, t_args->count);
			for (int i = start; i < end; ++i) {
				copy_individual(t_args->current_generation + i, t_args->next_generation + i);
			}
			t_args->cursor = t_args->count;
			// all threads wait to finish each other
			pthread_barrier_wait(t_args->barrier);

			// mutate first 20% children with the first version of bit string
			// mutation + parallel
			t_args->count = t_args->object_count * 2 / 10;
			start = t_args->id * (double)t_args->count / t_args->P;
			end = min((t_args->id + 1) * (double)t_args->count / t_args->P, t_args->count);
			for (int i = start; i < end; ++i) {
				copy_individual(t_args->current_generation + i, t_args->next_generation + t_args->cursor + i);
				mutate_bit_string_1(t_args->next_generation + t_args->cursor + i, k);
			}
			t_args->cursor += t_args->count;
			// all threads wait to finish each other
			pthread_barrier_wait(t_args->barrier);
			

			// mutate next 20% children with the second version of bit string
			// mutation + parallel
			t_args->count = t_args->object_count * 2 / 10;
			start = t_args->id * (double)t_args->count / t_args->P;
			end = min((t_args->id + 1) * (double)t_args->count / t_args->P, t_args->count);
			for (int i = start; i < end; ++i) {
				copy_individual(t_args->current_generation + i + t_args->count, t_args->next_generation + t_args->cursor + i);
				mutate_bit_string_2(t_args->next_generation + t_args->cursor + i, k);
			}
			t_args->cursor += t_args->count;
			pthread_barrier_wait(t_args->barrier);

			if(t_args->id == 0) {
				// crossover first 30% parents with one-point crossover
				// (if there is an odd number of parents, the last one is kept as such)
				t_args->count = t_args->object_count * 3 / 10;

				if (t_args->count % 2 == 1) {
					copy_individual(t_args->current_generation + t_args->object_count - 1, t_args->next_generation + t_args->cursor + t_args->count - 1);
					t_args->count--;
				}

				for (int i = 0; i < t_args->count; i += 2) {
					crossover(t_args->current_generation + i, t_args->next_generation + t_args->cursor + i, k);
				}

				// switch to new generation
				t_args->tmp = t_args->current_generation;
				t_args->current_generation = t_args->next_generation;
				t_args->next_generation = t_args->tmp;
			}
			
			start = t_args->id * (double)t_args->object_count / t_args->P;
			end = min((t_args->id + 1) * (double)t_args->object_count / t_args->P, t_args->object_count);
			for (int i = start; i < end; ++i) {
				t_args->current_generation[i].index = i;
			}
				
			if(t_args->id == 0) {
				//print_generation(t_args->current_generation, t_args->generations_count);
				if (k % 5 == 0) {
					print_best_fitness(t_args->current_generation);
				}
			}
			// all threads wait to finish each other
			pthread_barrier_wait(t_args->barrier);
		}
		// all threads wait to finish each other
		pthread_barrier_wait(t_args->barrier);

		// We did another parallelized compute fitness function and sorting
		// using qsort
		start = t_args->id * (double)t_args->object_count / t_args->P;
		end = min((t_args->id + 1) * (double)t_args->object_count / t_args->P, t_args->object_count);
		compute_fitness_function(t_args->objects, t_args->current_generation, end, t_args->sack_capacity, start);

		// keep first 30% children (elite children selection)
		qsort(t_args->current_generation + start, end - start, sizeof(individual), cmpfunc);
		pthread_barrier_wait(t_args->barrier);

		// merge the resulted subvectors, just like we did before
		if (t_args->id == 0) {
			individual* aux = (individual*)calloc(t_args->object_count, sizeof(individual));
			for (int i = 1; i < t_args->P; i++) {
				int new_start = i * (double)t_args->object_count / t_args->P;
				int new_end = min((i + 1) * (double)t_args->object_count / t_args->P, t_args->object_count);
				aux = merge(t_args->current_generation, 0, new_start, new_end, aux);
				memcpy(t_args->current_generation, aux, new_end * sizeof(individual));
				memset(aux, 0, t_args->object_count * sizeof(individual));
			}				
		}

		// all threads wait to finish each other
		pthread_barrier_wait(t_args->barrier);
		if(t_args->id == 0) {
			print_best_fitness(t_args->current_generation);
		}
	pthread_exit(NULL);
}


void run_genetic_algorithm(const sack_object *objects, int object_count, int generations_count, int sack_capacity, int P)
{
	int count = 0, cursor = 0;
	individual *current_generation = (individual*) calloc(object_count, sizeof(individual));
	individual *next_generation = (individual*) calloc(object_count, sizeof(individual));
	individual *tmp = NULL;

	// create the threads and their arguments 
	Thread_arguments* thread_args = (Thread_arguments *) calloc(P, sizeof(Thread_arguments));
	pthread_t* threads = (pthread_t*) calloc(P, sizeof(pthread_t));

	// create barrier
	pthread_barrier_t barrier;
	pthread_barrier_init(&barrier, NULL, P);

	int r;

	// initialize the threads and their corresponding arguments
	// we also save important variables, such as count, cursor, object_count, etc
	for (int i = 0; i < P; i++) {
		thread_args[i].id = i;
		thread_args[i].count = count;
		thread_args[i].cursor = cursor;
		thread_args[i].object_count = object_count;
		thread_args[i].generations_count = generations_count;
		thread_args[i].sack_capacity = sack_capacity;
		thread_args[i].objects = objects;
		thread_args[i].current_generation = current_generation;
		thread_args[i].next_generation = next_generation;
		thread_args[i].tmp = tmp;
		thread_args[i].barrier = &barrier;
		thread_args[i].P = P;

		// create the threads
		r = pthread_create(&threads[i], NULL, thread_function, &thread_args[i]);
	}

	// join threads
	for (int i = 0; i < P; i++) {
		r = pthread_join(threads[i], NULL);

		if (r) {
			printf("Eroare la asteptarea thread-ului %d\n", i);
			exit(-1);
		}
	}
	// destroy the barrier
	pthread_barrier_destroy(thread_args->barrier);

	// free resources for old generation
	free_generation(current_generation);
	free_generation(next_generation);

	// free resources
	free(current_generation);
	free(next_generation);
}