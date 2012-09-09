MRPrimes
Copyright 2012 Evan Brown
Released under GNU General Public License (see license.txt)
Contact: Evan Brown (ebrown255@gmail.com)

Description
-----------

MRPrimes is designed to allow for the generation of large prime numbers using
the probabalistic Miller-Rabin primality test. This means that the results are
not guaranteed to be prime, although using the default precision of 8 rounds
gives a certainty of at least 1 - (1/4) ^ 8 that the output is prime. See
http://en.wikipedia.org/wiki/Miller-Rabin. MRPrimes generates the probable
primes and prints them to a file ("primes.txt" by default).

Compilation
-----------

MRPrimes is dependent on the GNU Multiple Precision math library. If you do not
have this library installed, then you can find it at http://gmplib.org/.
MRPrimes also relies on POSIX multithreading.

It is recommended to compile MRPrimes with the GNU Compiler Collection. It is
possible to compile MRPrimes using:
  gcc -std=c99 -O3 -lgmp -lpthread -o mrprimes mrprimes.c

This will produce the executable file mrprimes in the current directory.

Usage
-----

MRPrimes can be used with a number of optional arguments.

[-v] or [--version] can be used to print the version of the software along
with licensing information.
example: ./mrprimes -v

[-h] or [--help] can be used to print a brief summary of the command options.
example: ./mrprimes -h

[-o] or [--output] can be used to specify the name of the output file.
example: ./mrprimes -o output.txt
This would result in the output being printed to the file "output.txt". The
default output filename is "primes.txt".

[-n] or [--numprimes] can be used to specify the number of primes to be
generated.
example: ./mrprimes -n 50
This would result in 50 primes being generated. The default number of primes
generated is 10.

[-d] or [--numdigits] can be used to specify the length of each prime to be
generated in digits.
example: ./mrprimes -d 500
This would result in the program generating 500 digit prime numbers. The
default number of digits used is 300.

[-s] or [--seed] can be used to specify a value for the random seed used by
the program to generate random starting points and to generate random numbers
for the Miller-Rabin test itself.
example: ./mrprimes -s 1
This would use the seed value 1. The default seed value calls the time(NULL)
function defined in <time.h>.

[-p] or [--precision] can be used to specify how many rounds of the Miller-
Rabin test to perform.
example: ./mrprimes -p 10
This would make the program perform 10 rounds of the test each time. The
default precision value is 8.

[-a] or [--append] can be used to tell the program to append to an output file
rather than overwriting it.
example: ./mrprimes -o out.txt -a
This would result in the program appending its results to the file "out.txt"
instead of overwriting it if it already exists. By default, the program
overwrites the output file rather than appending to it.

Explanation of Offsets
----------------------

The Miller-Rabin test operates on odd integers, but if it can be determined
that a given odd integer is divisible by some low prime number then performing
the test is superfluous. Determining whether a given odd integer, n, is
divisible by a given prime number, p, is simply a matter of taking the modulus
of n with respect to p and seeing whether or not it is equal to zero. The
result of this modulus, however, can be used not only to determine whather n is
divisible by p, but also how far n is from the next odd integer divisible by p.
This knowledge can be used to dramatically speed up the process of sequentially
testing the primality of odd integers.

Let us consider the case of the third odd prime number, 7. Every 7th integer
beginning with 7 is divisible by 7. Every 14th integer beginning with 7 is an
odd integer divisible by 7. Therefore, every 7th odd integer is divisible by 7.
This means that any random odd integer, n, can be classified as being of some
"offset" from the preceding odd integer divisible by 7, div, and this offset
can be related to the value of n mod 7.

n           n mod 7     offset
------------------------------
div         0           0
div + 2     2           1
div + 4     4           2
div + 6     6           3
div + 8     1           4
div + 10    3           5
div + 12    5           6
div + 14    0           0

For numbers more than 7 away from div, n mod 7 is 7 less than the difference
between n and div. Relating the value n mod 7 to the offset is simply a matter
of adding back in this 7 when appropriate and dividing by 2. The formula is
offset = (n mod 7 + ((n mod 7) mod 2) * 7) / 2
or generally:
offset = (n mod p + ((n mod p) mod 2) * p) / 2

MRPrimes finds and employs the offsets for the first 100 odd primes.
