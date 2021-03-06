/*
MRPrimes - a program to generate large prime numbers using the Miller-Rabin probabalistic
primality test implemented with the GNU Multiple Precision math library and POSIX threads.
Copyright (C) 2012, 2013 Evan Brown

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include "gmp.h"
#include "timer.h"

/* Create boolean type. */
enum boolean {FALSE, TRUE};
/* Set constants. */
enum constants {BASE = 10, THREAD_STACKSIZE = 1024};
static uint32_t *offset_primes = NULL; // to be initialized in main before threads are created

#define VERSION_NUMBER_STRING "1.0.7"

/* The following structure contains the necessary arguments to allow the threads to perform their function. */
struct thread_data_t
{
	const long num_digits;
	const long precision;
	const long num_offsets;
	long current_num_primes;
	char *out_file_name_pointer;
	/* Two randstates are necessary when multithreading to ensure the same primes are found when the same seed is used. */
	gmp_randstate_t random1; // for generating starting points
	gmp_randstate_t random2; // for generating random numbers for the test itself
	pthread_mutex_t current_num_primes_mutex;
	pthread_mutex_t out_file_mutex;
	pthread_mutex_t rand_state_mutex1;
	pthread_mutex_t rand_state_mutex2;
};


/* This function generates low odd primes for offsets at start of program run. */
static void
init_offsets (const long num_offsets)
{
	uint32_t np = 0, n = 3; // number of primes, first odd prime
	while (np < num_offsets)
	{
		enum boolean prime = TRUE;
		for (int i = 0; i < np; ++i)
		{
			if (n % offset_primes[i] == 0)
			{
				prime = FALSE;
				break;
			}
		}
		if (prime)
			offset_primes[np++] = n;
		n += 2; // next odd integer
	}
}

/* This function takes the remainder of a product and assigns the result to the first parameter (mpz_t acts as a reference). */
static void
mul_mod (mpz_t result, const mpz_t factor0, const mpz_t factor1, const mpz_t mod)
{
	mpz_mul (result, factor0, factor1);
	mpz_mod (result, result, mod);
}

/* This function uses the Miller-Rabin method to test the primality of an odd integer n with precision variable k. */
static enum boolean
miller_rabin (const mpz_t n, int k, gmp_randstate_t random, pthread_mutex_t *rand_state_mutex)
{
	/* Write n - 1 as 2^s*d with d odd by factoring powers of 2 from n - 1. */
	uint64_t s = 0;
	mpz_t d, a, x, tmp;
	mpz_inits (d, a, x, tmp, NULL);
	mpz_sub_ui (d, n, 1); // d = n - 1
	
	while (mpz_even_p (d))
	{
		mpz_tdiv_q_2exp (d, d, 1); // d /= 2
		++s;
	}
	
	for (int i = 0; i < k; ++i)
	{
		uint64_t r;
		
		/* Generate random number a in range [2, n - 2]. */
		mpz_set (tmp, n); // tmp = n
		mpz_sub_ui (tmp, tmp, 4); // tmp -= 4
		pthread_mutex_lock (rand_state_mutex);
		mpz_urandomm (a, random, tmp); // a in [0, n - 4]
		pthread_mutex_unlock (rand_state_mutex);
		mpz_add_ui (a, a, 2); // a in [2, n - 2]
		
		/* x = a^d % n */
		mpz_powm (x, a, d, n);
		
		mpz_add_ui (tmp, tmp, 3); // tmp = n - 1
		
		if (!mpz_cmp_ui (x, 1) || !mpz_cmp (x, tmp)) // if (x == 1 || x == n - 1)
			continue;
		for (r = 1; r < s; ++r)
		{
			mul_mod (x, x, x, n); // x = x * x % n
			if (!mpz_cmp_ui (x, 1)) // if (x == 1)
			{
				mpz_clears (d, a, x, tmp, NULL);
				return FALSE;
			}
			if (!mpz_cmp (x, tmp)) // if (x == n - 1)
				break;
		}
		if (r == s) // if the for loop completed
		{
			mpz_clears (d, a, x, tmp, NULL);
			return FALSE;
		}
	}
	mpz_clears (d, a, x, tmp, NULL);
	return TRUE;
}

/* This method generates a random integer with a certain number of digits in the parameter n to be used as a starting point. */
static void
gen_start (mpz_t n, int num_digits, gmp_randstate_t random, pthread_mutex_t *rand_state_mutex)
{
	char s[num_digits + 1];
	mpz_t tmp;
	mpz_init (tmp);
	
	/* Set tmp = 45 followed by (num_digits - 2) zeroes. */
	s[0] = '4';
	s[1] = '5';
	int i;
	for (i = 2; i < num_digits; ++i)
		s[i] = '0';
	s[num_digits] = '\0';
	mpz_set_str (tmp, (const char *)&s, BASE);
	
	pthread_mutex_lock (rand_state_mutex);
	mpz_urandomm (n, random, tmp); // n in [0, 44999...]
	pthread_mutex_unlock (rand_state_mutex);
	mpz_mul_2exp (n, n, 1); // n *= 2 -> n in [0, 8999... - 1]
	
	/* Set tmp = 1 followed by (num_digits - 2) zeroes. */
	s[0] = '1';
	s[1] = '0';
	s[i] = '\0'; // shorten string by one character
	mpz_set_str (tmp, (const char *)&s, BASE);
	
	mpz_add (n, n, tmp); // n in [1000... , 999... - 1] (set of all even integers with specified number of digits)
	mpz_add_ui (n, n, 1); // ++n
	mpz_clear (tmp);
}

/* This function initializes the offsets from the starting point (a random odd integer with the specified number of digits). */
static void
offset_init (const mpz_t start_point, const long num_offsets, int *offsets)
{
	for (int i = 0; i < num_offsets; ++i)
	{
		/* First take the starting point mod each low prime, then use this value to find the offset. See readme for explanation. */
		offsets[i] = (int)mpz_tdiv_ui (start_point, offset_primes[i]);
		offsets[i] = (offsets[i] + (offsets[i] % 2) * offset_primes[i]) / 2;
	}
}

/* This function updates the offsets after the test value has been incremented. */
static void
update_offsets (const long num_offsets, int *offsets)
{
	/* If the incremented offset is equal to the corresponding low prime, then it must be reset to 0. */
	for (int i = 0; i < num_offsets; ++i)
		if (++offsets[i] == offset_primes[i])
			offsets[i] = 0;
}

/* This function, based on strlen, tests whether any of the offsets is equal to 0. */
static enum boolean
any_offset_equals_zero (const long num_offsets, const int *start_offset)
{
	int *current_offset = (int *)start_offset;
	while (*current_offset) // loop until the value of current_offset is zero
		++current_offset;
	/* (current_offset - start_offset) will equal num_offsets if and only if none of the offsets equals 0. */
	return ((current_offset - start_offset) != num_offsets);
}

/* This function finds the next odd number which should be tested. */
static void
next_test (mpz_t test_value, const long num_offsets, int *offsets)
{
	while (any_offset_equals_zero (num_offsets, offsets))
	{
		mpz_add_ui (test_value, test_value, 2); // next odd integer
		update_offsets (num_offsets, offsets);
	}
}

/* This method defines the behavior of each thread: find one prime and print it to the output file. */
static void *
find_prime (void *thread_args)
{
	struct thread_data_t *data = (struct thread_data_t *)thread_args;
	enum boolean probably_prime;
	mpz_t test_value;
	mpz_init (test_value);
	FILE *out_file;
	
	/* Creating the offsets as an int array and storing a 0 at the end allows
	for the use of the function to test if any offsets are equal to 0. */
	int offsets[data->num_offsets + 1];
	offsets[data->num_offsets] = 0;
	
	/* Generate random starting position for search from the set of odd integers with the specified number of digits. */
	gen_start (test_value, data->num_digits, data->random1, &data->rand_state_mutex1);
	
	/* Keeping track of the offsets from odd integers divisible by low prime numbers allows for skipping the
	testing of odd numbers divisible by these low primes.  See readme for explanation of this principle. */
	offset_init (test_value, data->num_offsets, offsets);
	next_test (test_value, data->num_offsets, offsets);
	probably_prime = miller_rabin (test_value, data->precision, data->random2, &data->rand_state_mutex2);
	
	while (!probably_prime)
	{
		mpz_add_ui (test_value, test_value, 2); // next odd integer
		update_offsets (data->num_offsets, offsets);
		next_test (test_value, data->num_offsets, offsets);
		/* Moving the test to the end of the loop (which involves running the test once before the loop begins)
		ensures that no part of the loop executes unnecessarily once a result of TRUE has already been returned. */
		probably_prime = miller_rabin (test_value, data->precision, data->random2, &data->rand_state_mutex2);
	}
	
	/* Increment and print current number of primes found. */
	pthread_mutex_lock (&data->current_num_primes_mutex);
	printf ("Prime #%ld found\n", ++data->current_num_primes);
	pthread_mutex_unlock (&data->current_num_primes_mutex);
	
	/* Open output file, append prime, and close output file.  Opening and closing the output file each
	time allows the program to be aborted without losing the primes which have already been found. */
	pthread_mutex_lock (&data->out_file_mutex);
	out_file = fopen (data->out_file_name_pointer, "a");
	if (!out_file)
	{
		fprintf (stderr, "Error: failure to open output file.\n");
		exit (EXIT_FAILURE);
	}
	gmp_fprintf (out_file, "%Zd\n", test_value);
	fclose (out_file);
	pthread_mutex_unlock (&data->out_file_mutex);
	
	mpz_clear (test_value);
	pthread_exit (EXIT_SUCCESS);
}

/* Thus method prints the version number and a copyright message. */
static void
print_version ()
{
	printf ("\tMRPrimes %s\n", VERSION_NUMBER_STRING);
	printf ("\tCopyright (C) 2012, 2013 Evan Brown\n");
	printf ("\tLicense GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n");
	printf ("\tThis is free software: you are free to change and redistribute it.\n");
	printf ("\tThere is NO WARRANTY, to the extent permitted by law.\n");
}

/* This method prints a usage message. */
static void
print_help ()
{
	printf ("usage:\n");
	printf ("\t-o set output file\n");
	printf ("\t-n set number of primes to generate\n");
	printf ("\t-d set number of digits of primes to generate\n");
	printf ("\t-p set number of rounds of Miller-Rabin test to perform\n");
	printf ("\t-O set number of offset primes to generate\n");
	printf ("\t-s set random seed\n");
	printf ("\t-a set whether to append output to an existing file\n");
	printf ("\t-h print this help information\n");
	printf ("\t-v print program version information\n");
}

/* The main method generates prime numbers based on the command line arguments. */
int
main (int argc, char *argv[])
{
	/* Set default argument values. */
	char out_file_name[11] = "primes.txt";
	char *out_file_name_pointer = out_file_name; // output file name (-o)
	long num_digits = 300; // number of digits of primes to generate (-d)
	long num_primes = 10; // number of primes to generate (-n)
	long num_offsets = 10000; // number of offset primes (-O)
	int precision = 8; // rounds of Miller-Rabin test to perform (-p)
	uint64_t seed = (uint64_t)time (NULL); // random seed (-s)
	enum boolean append = FALSE; // whether the program should try to append output to an existing file (-a)
	
	/* Set argument values and check for validity. */
	if (argc > 1)
	{
		char **invalid_int = NULL; // for checking whether integer arguments are valid
		
		for (int i = 1; i < argc; ++i)
		{
			if (strcmp (argv[i], "-o") == 0 || strcmp (argv[i], "--output") == 0)
			{
				++i;
				if (i < argc)
				{
					out_file_name_pointer = argv[i];
				}
				else
				{
					fprintf (stderr, "Error: %s takes an argument. See readme for usage.\n", argv[i - 1]);
					return EXIT_FAILURE;
				}
			}
			else if (strcmp (argv[i], "-O") == 0 || strcmp (argv[i], "--numoffsets") == 0)
			{
				++i;
				if (i < argc)
				{
					num_offsets = strtol (argv[i], invalid_int, BASE);
					if (num_offsets <= 0 || invalid_int)
					{
						fprintf (stderr, "Error: number of offset primes must be a valid integer greater than 0.\n");
						return EXIT_FAILURE;
					}
				}
				else
				{
					fprintf (stderr, "Error: %s takes an argument. See readme for usage.\n", argv[i - 1]);
					return EXIT_FAILURE;
				}
			}
			else if (strcmp (argv[i], "-n") == 0 || strcmp (argv[i], "--numprimes") == 0)
			{
				++i;
				if (i < argc)
				{
					num_primes = strtol (argv[i], invalid_int, BASE);
					if (num_primes <= 0 || invalid_int)
					{
						fprintf (stderr, "Error: number of primes must be a valid integer greater than 0.\n");
						return EXIT_FAILURE;
					}
				}
				else
				{
					fprintf (stderr, "Error: %s takes an argument. See readme for usage.\n", argv[i - 1]);
					return EXIT_FAILURE;
				}
			}
			else if (strcmp (argv[i], "-d") == 0 || strcmp (argv[i], "--numdigits") == 0)
			{
				++i;
				if (i < argc)
				{
					num_digits = strtol (argv[i], invalid_int, BASE);
					if (num_digits < 10 || invalid_int)
					{
						fprintf (stderr, "Error: number of digits must be a valid integer greater than or equal to 10.\n");
						return EXIT_FAILURE;
					}
				}
				else
				{
					fprintf (stderr, "Error: %s takes an argument. See readme for usage.\n", argv[i - 1]);
					return EXIT_FAILURE;
				}
			}
			else if (strcmp (argv[i], "-s") == 0 || strcmp (argv[i], "--seed") == 0)
			{
				++i;
				if (i < argc)
				{
					/* Create a temporary long int to test for a negative value. */
					long seed_temp = strtol (argv[i], invalid_int, BASE);
					if (seed_temp < 0 || invalid_int)
					{
						fprintf (stderr, "Error: seed value must be a valid long integer greater than or equal to 0.\n");
						return EXIT_FAILURE;
					}
					seed = seed_temp;
				}
				else
				{
					fprintf (stderr, "Error: %s takes an argument. See readme for usage.\n", argv[i - 1]);
					return EXIT_FAILURE;
				}
			}
			else if (strcmp (argv[i], "-p") == 0 || strcmp (argv[i], "--precision") == 0)
			{
				++i;
				if (i < argc)
				{
					precision = strtol (argv[i], invalid_int, BASE);
					if (precision <= 0 || precision >= 200 || invalid_int)
					{
						fprintf (stderr, "Error: Miller Rabin test precision must be a valid integer greater than 0 and less than 200.\n");
						return EXIT_FAILURE;
					}
				}
				else
				{
					fprintf (stderr, "Error: %s takes an argument. See readme for usage.\n", argv[i - 1]);
					return EXIT_FAILURE;
				}
			}
			else if (strcmp (argv[i], "-a") == 0 || strcmp (argv[i], "--append") == 0)
			{
				append = TRUE;
			}
			else if (strcmp (argv[i], "-v") == 0 || strcmp (argv[i], "--version") == 0)
			{
				print_version ();
				return EXIT_SUCCESS;
			}
			else if (strcmp (argv[i], "-h") == 0 || strcmp (argv[i], "--help") == 0)
			{
				print_help ();
				return EXIT_SUCCESS;
			}
			else
			{
				fprintf (stderr, "Error: invalid arguments. See readme for usage.\n");
				return EXIT_FAILURE;
			}
		}
	}
	
	/* Get start time. */
	timer ();
	
	/* Initialize offset primes. */
	offset_primes = malloc (num_offsets * sizeof (uint32_t));
	init_offsets (num_offsets);
	
	/* Delete contents of output file unless user specified otherwise. */
	if (!append)
	{
		FILE *out_file = fopen (out_file_name_pointer, "w");
		fclose (out_file);
	}
	
	/* Initialize and set thread detached and stacksize attributes. */
	pthread_attr_t attr;
	pthread_attr_init (&attr);
	pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
	pthread_attr_setstacksize (&attr, THREAD_STACKSIZE);
	
	/* Initialize thread arguments. */
	struct thread_data_t thread_args = {num_digits, precision, num_offsets}; // num_digits, precision, num_offsets must be initialized immediately because they are const
	thread_args.out_file_name_pointer = out_file_name_pointer;
	gmp_randinit_mt (thread_args.random1);
	gmp_randseed_ui (thread_args.random1, seed);
	gmp_randinit_mt (thread_args.random2);
	gmp_randseed_ui (thread_args.random2, seed);
	thread_args.current_num_primes = 0;
	pthread_mutex_init (&thread_args.out_file_mutex, NULL);
	pthread_mutex_init (&thread_args.current_num_primes_mutex, NULL);
	pthread_mutex_init (&thread_args.rand_state_mutex1, NULL);
	pthread_mutex_init (&thread_args.rand_state_mutex2, NULL);
	
	/* Print initialization time. */
	if (CLOCK_PRECISION == 9)
		printf ("Initialization time: %.9lf seconds.\n", timer ());
	else
		printf ("Initialization time: %.6lf seconds.\n", timer ());
	
	/* Create one thread to find each prime. */
	pthread_t threads[num_primes];
	int return_code;
	for (int i = 0; i < num_primes; ++i)
	{
		return_code = pthread_create (&threads[i], &attr, find_prime, (void *)&thread_args);
		if (return_code)
		{
			fprintf (stderr, "Error: return code from pthread_create is %d\n", return_code);
			return EXIT_FAILURE;
		}
	}
	
	/* Free attribute and wait for the threads to finish. */
	pthread_attr_destroy (&attr);
	for (int i = 0; i < num_primes; ++i)
	{
		return_code = pthread_join (threads[i], NULL);
		if (return_code)
		{
			fprintf (stderr, "Error: return code from pthread_join is %d\n", return_code);
			return EXIT_FAILURE;
		}
	}
	
	/* Get end time and print time taken. */
	if (CLOCK_PRECISION == 9)
		printf ("Execution time: %.9lf seconds.\n", timer ());
	else
		printf ("Execution time: %.6lf seconds.\n", timer ());
	
	/* Cleanup and exit. */
	free (offset_primes);
	gmp_randclear (thread_args.random1);
	gmp_randclear (thread_args.random2);
	pthread_mutex_destroy (&thread_args.out_file_mutex);
	pthread_mutex_destroy (&thread_args.current_num_primes_mutex);
	pthread_mutex_destroy (&thread_args.rand_state_mutex1);
	pthread_mutex_destroy (&thread_args.rand_state_mutex2);
	pthread_exit (EXIT_SUCCESS);
}
