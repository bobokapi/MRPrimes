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

#include <time.h>
#include <sys/time.h>

/* Initialize time structs and precision. */
#ifdef CLOCK_MONOTONIC
	#define CLOCK_PRECISION 9
	static struct timespec ts[2];
#else
	#define CLOCK_PRECISION 6
	static struct timeval tv[2];
#endif

/* Keep track of which of the two structs is from the last function call. */
char state = -1; // -1 = first run state

/* This function gets the current time, and returns the difference from the last call. */
static double
timer ()
{
	double dtime = 0.0;
	#ifdef CLOCK_MONOTONIC
		clock_gettime (CLOCK_MONOTONIC, &ts[!state]);
		if (state != -1)
		{
			dtime = ts[!state].tv_sec + ts[!state].tv_nsec * 1e-9;
			dtime -= ts[state].tv_sec + ts[state].tv_nsec * 1e-9;
		}
	#else
		gettimeofday (&tv[!state], NULL);
		if (state != -1)
		{
			dtime = tv[!state].tv_sec + tv[!state].tv_usec * 1e-6;
			dtime -= tv[state].tv_sec + tv[state].tv_usec * 1e-6;
		}
	#endif
	state = !state; // -1,+1 -> 0; 0 -> +1
	return dtime;
}
