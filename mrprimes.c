/* This is a program designed by Evan Brown which generates large prime numbers using the Miller-Rabin probabalistic
 * primality test implemented with the GNU Multiple Precision math library and POSIX threads */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <pthread.h>
#include "gmp.h"

#define BASE 10
#define THREAD_STACKSIZE 1024
#define USECS_PER_SEC 1000000
enum boolean {FALSE, TRUE};

/* The following structure contains the offsets from numbers divisible by the first ten odd primes */
struct offsets_t {
	char off3, off5, off7, off11, off13, off17, off19, off23, off27, off29;
};

/* The following structure contains the necessary arguments to allow the threads to perform their function */
struct threadData_t {
	const long numDigits;
	const long precision;
	long currentNumPrimes;
	char *outFileNamePointer;
	/* two randstates are necessary when multithreading to ensure the same primes are found when the same seed is used */
	gmp_randstate_t random1; // for generating starting points
	gmp_randstate_t random2; // for generating random numbers for the test itself
	pthread_mutex_t currentNumPrimesMutex;
	pthread_mutex_t outFileMutex;
	pthread_mutex_t randStateMutex1;
	pthread_mutex_t randStateMutex2;
};

/* This function takes the remainder of a product and assigns the result to the first parameter (mpz_t acts as a reference) */
static void mulMod(mpz_t result, const mpz_t factor0, const mpz_t factor1, const mpz_t mod) {
	mpz_mul(result, factor0, factor1);
	mpz_mod(result, result, mod);
}

/* This function uses the Miller-Rabin method to test the primality of an odd integer n with precision variable k */
static enum boolean millerRabin(const mpz_t n, int k, gmp_randstate_t random, pthread_mutex_t *randStateMutex) {
	/* write n - 1 as 2^s*d with d odd by factoring powers of 2 from n - 1 */
	unsigned long s = 0;
	mpz_t d, a, x, tmp;
	mpz_inits(d, a, x, tmp, NULL);
	mpz_sub_ui(d, n, 1); // d = n - 1
	while (mpz_even_p(d)) {
		mpz_tdiv_q_2exp(d, d, 1); // d /= 2
		++s;
	}
	
	for (int i = 0; i < k; ++i) {
		unsigned long r;
		
		/* generate random number a in range [2, n - 2] */
		mpz_set(tmp, n); // tmp = n
		mpz_sub_ui(tmp, tmp, 4); // tmp -= 4
		pthread_mutex_lock(randStateMutex);
		mpz_urandomm(a, random, tmp); // a in [0, n - 4]
		pthread_mutex_unlock(randStateMutex);
		mpz_add_ui(a, a, 2); // a in [2, n - 2]
		
		/* x = a^d % n */
		mpz_powm(x, a, d, n);
		
		mpz_add_ui(tmp, tmp, 3); // tmp = n - 1
		
		if (!mpz_cmp_ui(x, 1) || !mpz_cmp(x, tmp)) continue; // if (x == 1 || x == n - 1)
		for (r = 1; r < s; ++r) {
			mulMod(x, x, x, n); // x = x * x % n
			if (!mpz_cmp_ui(x, 1)) { // if (x == 1)
				mpz_clears(d, a, x, tmp, NULL);
				return FALSE;
			}
			if (!mpz_cmp(x, tmp)) break; // if (x == n - 1)
		}
		if (r == s) { // if the for loop completed
			mpz_clears(d, a, x, tmp, NULL);
			return FALSE;
		}
	}
	mpz_clears(d, a, x, tmp, NULL);
	return TRUE;
}

/* This method generates a random integer with a certain number of digits in the parameter n to be used as a starting point */
static void genStart(mpz_t n, int numDigits, gmp_randstate_t random, pthread_mutex_t *randStateMutex) {
	char s[numDigits + 1];
	mpz_t tmp;
	mpz_init(tmp);
	
	/* set tmp = 45 followed by (numDigits - 2) zeroes */
	s[0] = '4';
	s[1] = '5';
	int i;
	for (i = 2; i < numDigits; ++i) s[i] = '0';
	s[numDigits] = '\0';
	mpz_set_str(tmp, (const char *)&s, BASE);
	
	pthread_mutex_lock(randStateMutex);
	mpz_urandomm(n, random, tmp); // n in [0, 44999...]
	pthread_mutex_unlock(randStateMutex);
	mpz_mul_2exp(n, n, 1); // n *= 2 -> n in [0, 8999... - 1]
	
	/* set tmp = 1 followed by (numDigits - 2) zeroes */
	s[0] = '1';
	s[1] = '0';
	s[i] = '\0';
	mpz_set_str(tmp, (const char *)&s, BASE);
	
	mpz_add(n, n, tmp); // n in [1000... , 999... - 1] (set of all even integers with specified number of digits)
	mpz_add_ui(n, n, 1); // ++n
	mpz_clear(tmp);
}

/* This function initializes the offsets from the starting point (a random odd integer with the specified number of digits) */
static void offsetInit(const mpz_t startPoint, struct offsets_t *offsets) {
	/* First take the starting point mod each low prime, then use this value to find the offset. See readme for explanation */
	offsets->off3 = (char)mpz_tdiv_ui(startPoint, 3);
	offsets->off3 = (offsets->off3 + (offsets->off3 % 2) * 3) / 2;
	
	offsets->off5 = (char)mpz_tdiv_ui(startPoint, 5);
	offsets->off5 = (offsets->off5 + (offsets->off5 % 2) * 5) / 2;
	
	offsets->off7 = (char)mpz_tdiv_ui(startPoint, 7);
	offsets->off7 = (offsets->off7 + (offsets->off7 % 2) * 7) / 2;
	
	offsets->off11 = (char)mpz_tdiv_ui(startPoint, 11);
	offsets->off11 = (offsets->off11 + (offsets->off11 % 2) * 11) / 2;
	
	offsets->off13 = (char)mpz_tdiv_ui(startPoint, 13);
	offsets->off13 = (offsets->off13 + (offsets->off13 % 2) * 13) / 2;
	
	offsets->off17 = (char)mpz_tdiv_ui(startPoint, 17);
	offsets->off17 = (offsets->off17 + (offsets->off17 % 2) * 17) / 2;
	
	offsets->off19 = (char)mpz_tdiv_ui(startPoint, 19);
	offsets->off19 = (offsets->off19 + (offsets->off19 % 2) * 19) / 2;
	
	offsets->off23 = (char)mpz_tdiv_ui(startPoint, 23);
	offsets->off23 = (offsets->off23 + (offsets->off23 % 2) * 23) / 2;
	
	offsets->off27 = (char)mpz_tdiv_ui(startPoint, 27);
	offsets->off27 = (offsets->off27 + (offsets->off27 % 2) * 27) / 2;
	
	offsets->off29 = (char)mpz_tdiv_ui(startPoint, 29);
	offsets->off29 = (offsets->off29 + (offsets->off29 % 2) * 29) / 2;
}

/* This function updates the offsets after the test value has been incremented */
static void updateOffsets(struct offsets_t *offsets) {
	offsets->off3 = (offsets->off3 + 1) % 3;
	offsets->off5 = (offsets->off5 + 1) % 5;
	offsets->off7 = (offsets->off7 + 1) % 7;
	offsets->off11 = (offsets->off11 + 1) % 11;
	offsets->off13 = (offsets->off13 + 1) % 13;
	offsets->off17 = (offsets->off17 + 1) % 17;
	offsets->off19 = (offsets->off19 + 1) % 19;
	offsets->off23 = (offsets->off23 + 1) % 23;
	offsets->off27 = (offsets->off27 + 1) % 27;
	offsets->off29 = (offsets->off29 + 1) % 29;
}

/* This function finds the next odd number which should be tested */
static void nextTest(mpz_t testValue, struct offsets_t *offsets) {
	while (!offsets->off3 || !offsets->off5 || !offsets->off7 || !offsets->off11 || !offsets->off13 ||
		   !offsets->off17 || !offsets->off19 || !offsets->off23 || !offsets->off27 || !offsets->off29) {
		mpz_add_ui(testValue, testValue, 2); // next odd integer
		updateOffsets(offsets);
	}
}

/* This method defines the behavior of each thread in finding primes and printing them to the output file */
void *findPrime(void *args) {
	struct threadData_t *data = (struct threadData_t *)args;
	enum boolean probablyPrime;
	mpz_t testValue;
	mpz_init(testValue);
	FILE *outFile;
	struct offsets_t offsets;
	
	/* generate random starting position for search from the set of odd integers with the specified number of digits */
	genStart(testValue, data->numDigits, data->random1, &data->randStateMutex1);
	
	/* Keeping track of the offsets from odd integers divisible by low prime numbers allows for skipping the 
	 * testing of odd numbers divisible by these low primes. See readme for explanation of this principle. */
	offsetInit(testValue, &offsets);
	nextTest(testValue, &offsets);
	probablyPrime = millerRabin(testValue, data->precision, data->random2, &data->randStateMutex2);
	
	while (!probablyPrime) {
		mpz_add_ui(testValue, testValue, 2); // next odd integer
		updateOffsets(&offsets);
		nextTest(testValue, &offsets);
		/* Moving the test to the end of the loop (which involves running the test once before the loop begins)
		 * ensures that no part of the loop executes unnecessarily once a result of TRUE has already been returned. */
		probablyPrime = millerRabin(testValue, data->precision, data->random2, &data->randStateMutex2);
	}
	
	/* print and increment current number of primes found */
	pthread_mutex_lock(&data->currentNumPrimesMutex);
	printf("Prime #%d found\n", ++data->currentNumPrimes);
	pthread_mutex_unlock(&data->currentNumPrimesMutex);
	
	/* Open output file, append prime, and close output file. Opening and closing the output file each time
	 * allows the program to be aborted without losing the primes which have already been found. */
	pthread_mutex_lock(&data->outFileMutex);
	outFile = fopen(data->outFileNamePointer, "a");
	if (!outFile) {
		fprintf(stderr, "Error: failure to open output file.\n");
		exit(EXIT_FAILURE);
	}
	gmp_fprintf(outFile, "%Zd\n", testValue);
	fclose(outFile);
	pthread_mutex_unlock(&data->outFileMutex);
	
	mpz_clear(testValue);
	pthread_exit(EXIT_SUCCESS);
}

/* The main method generates prime numbers based on the command line arguments */
int main(int argc, char *argv[]) {
	/* set default argument values */
	char outFileName[11] = "primes.txt";
	char *outFileNamePointer = outFileName;
	long numDigits = 300;
	long numPrimes = 10;
	int precision = 8;
	long seed = time(NULL);
	enum boolean append = FALSE; // whether the program should try to append output to an existing output file
	
	/* set argument values and check for validity */
	if (argc > 1) {
		char **invalidInt = NULL; // for checking whether integer arguments are valid
		
		for (int i = 1; i < argc; ++i) {
			if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
				++i;
				if (i < argc) {
					outFileNamePointer = argv[i];
				} else {
					fprintf(stderr, "Error: %s takes an argument. See readme for usage.\n", argv[i - 1]);
					exit(EXIT_FAILURE);
				}
			}
			else if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--numprimes") == 0) {
				++i;
				if (i < argc) {
					numPrimes = (int)strtol(argv[i], invalidInt, BASE);
					if (numPrimes <= 0 || invalidInt) {
						fprintf(stderr, "Error: number of primes must be a valid integer greater than 0.\n");
						exit(EXIT_FAILURE);
					}
				} else {
					fprintf(stderr, "Error: %s takes an argument. See readme for usage.\n", argv[i - 1]);
					exit(EXIT_FAILURE);
				}
			}
			else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--numdigits") == 0) {
				++i;
				if (i < argc) {
					numDigits = (int)strtol(argv[i], invalidInt, BASE);
					if (numDigits <= 2 || invalidInt) {
						fprintf(stderr, "Error: number of digits must be a valid integer greater than 2.\n");
						exit(EXIT_FAILURE);
					}
				} else {
					fprintf(stderr, "Error: %s takes an argument. See readme for usage.\n", argv[i - 1]);
					exit(EXIT_FAILURE);
				}
			}
			else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--seed") == 0) {
				++i;
				if (i < argc) {
					seed = strtol(argv[i], invalidInt, BASE);
					if (invalidInt) {
						fprintf(stderr, "Error: seed value must be a valid long integer.\n");
						exit(EXIT_FAILURE);
					}
				} else {
					fprintf(stderr, "Error: %s takes an argument. See readme for usage.\n", argv[i - 1]);
					exit(EXIT_FAILURE);
				}
			}
			else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--precision") == 0) {
				++i;
				if (i < argc) {
					precision = strtol(argv[i], invalidInt, BASE);
					if (precision <= 0 || precision >= 200 || invalidInt) {
						fprintf(stderr, "Error: Miller Rabin test precision must be a valid integer greater than 0 and less than 200.\n");
						exit(EXIT_FAILURE);
					}
				} else {
					fprintf(stderr, "Error: %s takes an argument. See readme for usage.\n", argv[i - 1]);
					exit(EXIT_FAILURE);
				}
			}
			else if (strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--append") == 0) {
				append = TRUE;
			}
			else {
				fprintf(stderr, "Error: invalid arguments. See readme for usage.\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	
	/* delete contents of output file unless user specified otherwise */
	if (!append) {
		FILE *outFile = fopen(outFileNamePointer, "w");
		fclose(outFile);
	}
	
	/* initialize and set thread detached and stacksize attributes */
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	pthread_attr_setstacksize(&attr, THREAD_STACKSIZE);
	
	/* initialize thread arguments */
	struct threadData_t args = {numDigits, precision}; // numDigits and precision must be initialized immediately because they are const
	args.outFileNamePointer = outFileNamePointer;
	gmp_randinit_mt(args.random1);
	gmp_randseed_ui(args.random1, seed);
	gmp_randinit_mt(args.random2);
	gmp_randseed_ui(args.random2, seed);
	args.currentNumPrimes = 0;
	pthread_mutex_init(&args.outFileMutex, NULL);
	pthread_mutex_init(&args.currentNumPrimesMutex, NULL);
	pthread_mutex_init(&args.randStateMutex1, NULL);
	pthread_mutex_init(&args.randStateMutex2, NULL);
	
	/* initialize time vars and get time */
	struct timeval startTime, endTime;
	gettimeofday(&startTime, NULL);
	
	/* create one thread to find each prime */
	pthread_t threads[numPrimes];
	int rc;
	for (int i = 0; i < numPrimes; ++i) {
		rc = pthread_create(&threads[i], &attr, findPrime, (void *)&args);
		if (rc) {
			fprintf(stderr, "Error: return code from pthread_create() is %d\n", rc);
			exit(EXIT_FAILURE);
		}
	}
	
	/* free attribute and wait for the threads to finish */
	pthread_attr_destroy(&attr);
	for(int i = 0; i < numPrimes; ++i) {
		rc = pthread_join(threads[i], NULL);
		if (rc) {
			fprintf(stderr, "Error: return code from pthread_join() is %d\n", rc);
			exit(EXIT_FAILURE);
		}
	}
	
	/* get end time and print time taken */
	gettimeofday(&endTime, NULL);
	int usec;
	if (endTime.tv_usec > startTime.tv_usec) usec = endTime.tv_usec - startTime.tv_usec;
	else {
		usec = USECS_PER_SEC + endTime.tv_usec - startTime.tv_usec;
		--endTime.tv_sec;
	}
	printf("Execution time: %ld seconds, %d microseconds.\n", endTime.tv_sec - startTime.tv_sec, usec);
	
	/* clear thread arguments and exit */
	gmp_randclear(args.random1);
	gmp_randclear(args.random2);
	pthread_mutex_destroy(&args.outFileMutex);
	pthread_mutex_destroy(&args.currentNumPrimesMutex);
	pthread_mutex_destroy(&args.randStateMutex1);
	pthread_mutex_destroy(&args.randStateMutex2);
	pthread_exit(EXIT_SUCCESS);
}
