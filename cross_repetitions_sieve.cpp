//============================================================================
// Name        : cross_repetitions_sieve.cpp
// Author      : Felix Baessler
// Version     : 01.04.2022
// Copyright   : Felix Baessler, felix.baessler@gmail.com
// SEE TLDR; VERSION OF THE LICENSE: https://creativecommons.org/licenses/by-nc/4.0/legalcode
// SEE FULL LICENSE DETAILS HERE   : https://creativecommons.org/licenses/by-nc/4.0/
// Description : Repeated Substrings / cross_repetitions_sieve
// Compilation flags:
// -O3 -g3 -Wall         : optimization
// -Wl,--stack,0xFFFFFF  : long arrays
// if required: upgrade minGW to x86_64 : 64 bit executable
//
// introduction available on: https://sites.google.com/view/repsieve
//============================================================================

/*
Repeated Substrings : Cross-Repetitions Sieve
=============================================

This program illustrates the principles of a sieve, that solves the following problem:

Given two strings, find all comparatively short substrings that have at least one
repetition in both, string S1 and string S2.

For the purpose of the experiment, the long input strings S1 and S2 (~ Giga Bytes) are generated
by means of the standard MT19937 pseudo-random number generator. Also to make sure that
at least one repetition can be found, a test duplicate (~ 7 Bytes) is copied from S1 to S2.

The conjecture is: as long as the (minimal) string length of the repetitions we are looking for,
are longer than the length of the majority of the inherent random repetitions, the sieve can be
competitive:
 - in terms of space, as only three bit-vectors of the lengths S1, S2 are required
 - in terms of time, as experimental results indicate a linear behavior
The question, how well an execution environment can cope with an unpredictably accessed,
very long bit-vector (exceeding the processor cache), requires hardware dependent analysis.

The method presented below, relies essentially on:
 -	string shingling
	see for example "Mining of Massive Datasets", by Leskovec, Rajaraman and Ullman
 -	the Karp-Rabin signature algorithm (a.k.a. the "fingerprint" method)
	see "Efficient randomized pattern-matching algorithms", by Karp and Rabin
 -	the Sattolo shuffle algorithm (to generate uniformly distributed cycles)

Shingling of Strings
Shingling of strings is a simple technique to represent strings as sets of substrings,
in our case, for the purpose of identifying repeated substrings.
Define a shingle for a string to be any substring of length L found within a given string.
Instead of using shingles directly, a well-known approach is to use a hash function
that maps shingles of several bytes length to a number represented by only a few bytes.
This short representative of a shingle can be viewed as its fingerprint or signature and
L is the number of successive bytes used to compute it.
Collisions occur when shingles, identical or not, have the same signature.
Colliding shingles that are not identical are called false positives.

Shingle Sieve
=============
The shingle sieve operates as a two-way signature filter:
 -	forward  : the shingles are processed from left to right
 -	backward : the shingles are processed from right to left
Shingles that produce unique signatures fall through the sieve and are eliminated.
Pre-shuffling at the beginning of each sieving round boosts performance.

forward shingling: the hash is computed from left to right
----------------------------------------------------------
N: length of input string S

	01234567890 ...
1	xxxxxxx
2	 xxxxxxx      L = 7  :  shingle length
3	  xxxxxxx
4	   ...

j         : shingle begin  (on the left)      ==> shingle start_index (left adjusted)
j + L - 1 : shingle end    (on the right)
N - L     : last, rightmost shingle begin index

backward shingling: the hash is computed from right to left
-----------------------------------------------------------
N: length of input string S

... 012345678901234567890
1                 xxxxxxx
2                xxxxxxx      L = 7  :  shingle length
3               xxxxxxx
4                  ...

i           : shingle begin (on the right)
i - L + 1   : shingle end   (on the left)     ==> shingle start_index (left adjusted!)
L -	1       : last, leftmost shingle begin index

Note that the start index of a shingle is defined by the index of its leftmost byte
 -	on forward  processing: the begin location
 -	on backward processing: the end   location
This allows simple switching between the state vectors when the processing direction is reversed.

Start Shingles
--------------
For repetitions longer than L, the sieve will retain a whole sequence of contiguous shingles.
In this case, it is more economic to record only the start shingle index, instead of storing individually
the indices of each of the involved shingles.

Looking more closely into the sieving process, three types of shingles can be distinguished:

T-Shingles (T: trigger)
----------
 -	a T-shingle is the first shingle to produce a beforehand unknown signature.
 -	a T-shingle may trigger one or more collisions.
 -	T-shingles will be tracked by the state vector t.

X-Shingles (X: collision)
----------
 -	an X-shingle is a later processed shingle that collides with a T-shingle
 	and possibly other X-shingles.
 -	an X-shingle may be a repetition of a T-shingle (mind false positives).
 -	X-shingles will be tracked by the state vector x.

Collision free shingles can be eliminated from further processing, as their signature
cannot be identical with the signature of any other shingle in the sieve.

E-Shingles (E: eliminated)
----------
 - an E-shingle is a shingle that is eliminated from the sieve
 - E-shingles are unique, i.e. free of repetition

State Vectors
-------------
Both state vectors, t and x, are bit-vectors of the same lengths as S1, S2.
Along with the input string, they are sequentially accessed during the sieving process.
Trigger Vector:
 -	t[j] = true : the shingle j is a  T-shingle
 -	t[j] = false: the shingle j is an E-shingle
Collision Vector:
 -	x[j] = true : the shingle j is an X-shingle, it MAY repeat with a T-shingle
 -	x[j] = false: the shingle j is an E-shingle
At the beginning of the sieving process, t and x are true for all shingles.
Moving forward  (t1,x2): the T-shingles are in S1 and the X-shingles are in S2.
Reversing the processing direction: t1 becomes x1 and x2 becomes t2.
Moving backward (t2,x1): the T-shingles are now in S2 and the X-shingles in S1.

The shingle signatures, generated during the sieving process, will be recorded
by means of the signature incident vector.

Signature Incident Vector
-------------------------
The signature incident bit-vector z reflects the signature production history
during the hashing process:
 -	z[h] = true : the signature h has already been computed for one or more shingles
 -	z[h] = false: the signature h has never been produced beforehand
At each beginning of forward and backward sieving, z is initialized to false.
Note that the objective is to keep z small; a reasonable choice is to use a z of
the same length as S1, assuming S1 smaller than S2.

State Check
-----------
During the sieving process, it is always possible to check whether a test pair (j1,j2)
is still in the sieve:
if (t1[j1] && x2[j2]) printf("The test pair is still in the sieve! \n");
where j1, j2 denote the respective indices of the test pair in S1, S2.

Sieve Filter
------------
	if (TestBit(x, j)) {
		// shingle j is an X-shingle
		if (!TestBit(z, hash))  {
			// shingle j has a unique signature
			uniques++;
			ClearBit(x, j);
		} else {
			// collisions++;
		}
	} else {
		// shingle j is an E-shingle (elimination)
	}


Program Usage
=============
The program parameters N1, N2, M and L are set by define statements:

 -  N1, N2 are the lengths of S1, S2, i.e. the number of random bytes
    generated at the beginning of the program.

 -	M is the modulus of the Rabin fingerprint;
	it defines the length of the signature incident vector

 -	L is the length of the shingles;
	it is also used to set the string length of the test pairs.

Two additional define statements concern:

 -  B = 257, the base of the Rabin fingerprint (first prime > 256)
 -  C = (B ^ L) % M, a derived constant to help compiler optimization

For given N1, N2, the performance of the sieve depends on the choice of M and L.

Prime numbers can be obtained from: https://www.walter-fendt.de/html5/mde/primenumbers_de.htm

Sample OUTPUT ( 1/10 and 1 Giga Byte)
=====================================
#define B   257
#define N1   100000000ULL	// 100 MB
#define N2  1000000000ULL	//   1 GB
#define M     99999989ULL
#define L   7
#define C     77150229ULL

parameter check
storage allocation

Repeated Substrings Cross-Sieve S1 x S2
=======================================
base               B: 257
modulus            M: 99999989
substrings length  L: 7
(B ^ L) % M        C: 77150229
string S1 length  N1: 100000000
string S2 length  N2: 1000000000
E(S1 x S2)        Ec:  1.39  expected cross-repetitions between S1 and S2
E(movements)      Em: 13.25  expected number of required sieve movements

start generation of input string
first, last bytes of input string: 16 147
overlapping end   of S1: 67 ... 212
overlapping begin of S2: 67 ... 212
test repetitions indices: 50000000 500000000
prime_performance_test relative hash time: 0.004727 microseconds 	 99999994 	 160121555
prime_performance_test relative hash time: 0.004699 microseconds 	 999999994 	 1012018771

ROUND 1
  > S1 forward  :     2010.85 milliseconds
  > S2 forward  :    22938.88 milliseconds 	 residue  : 632150014 (36.784998 %)
elapsed time   1:    24999.68 milliseconds
ROUND 2
  > S1 forward  :     1961.85 milliseconds
  > S2 forward  :    20648.19 milliseconds 	 residue  : 399603347 (36.786627 %)
elapsed time   2:    47611.73 milliseconds
ROUND 3
  > S1 forward  :     1947.85 milliseconds
  > S2 forward  :    17061.25 milliseconds 	 residue  : 252580441 (36.792211 %)
elapsed time   3:    66622.84 milliseconds
ROUND 4
  > S1 forward  :     1947.87 milliseconds
  > S2 forward  :    13512.26 milliseconds 	 residue  : 159667435 (36.785511 %)
elapsed time   4:    82084.98 milliseconds
ROUND 5
  > S1 forward  :     2004.82 milliseconds
  > S2 forward  :    11149.63 milliseconds 	 residue  : 100929068 (36.787944 %)
elapsed time   5:    95241.45 milliseconds
  < S2 backward :     8388.19 milliseconds
  < S1 backward :     2288.69 milliseconds 	 residue  : 63558166 (36.441830 %)
elapsed time   6:   105920.34 milliseconds
ROUND 6
  > S1 forward  :     1642.06 milliseconds
  > S2 forward  :     9971.29 milliseconds 	 residue  : 47468821 (52.968137 %)
elapsed time   7:   117535.68 milliseconds
  < S2 backward :     6878.06 milliseconds
  < S1 backward :     2149.77 milliseconds 	 residue  : 24017430 (62.211889 %)
elapsed time   8:   126565.52 milliseconds
ROUND 7
  > S1 forward  :     1100.36 milliseconds
  > S2 forward  :     7291.80 milliseconds 	 residue  : 10134684 (78.649809 %)
elapsed time   9:   134959.69 milliseconds
  < S2 backward :     5360.95 milliseconds
  < S1 backward :     1234.29 milliseconds 	 residue  : 2314713 (90.362362 %)
elapsed time  10:   141556.92 milliseconds
ROUND 8
  > S1 forward  :      595.63 milliseconds
  > S2 forward  :     5420.89 milliseconds 	 residue  : 231850 (97.712312 %)
elapsed time  11:   147575.45 milliseconds
  < S2 backward :     4702.33 milliseconds
  < S1 backward :      613.65 milliseconds 	 residue  : 5272 (99.772240 %)
elapsed time  12:   152893.44 milliseconds
ROUND 9
  > S1 forward  :      473.73 milliseconds
  > S2 forward  :     4695.31 milliseconds 	 residue  : 15 (99.993530 %)
elapsed time  13:   158064.45 milliseconds
  < S2 backward :     4698.32 milliseconds
  < S1 backward :      474.73 milliseconds 	 residue  : 3 (99.943096 %)
elapsed time  14:   163239.50 milliseconds
ROUND 10
  > S1 forward  :      475.72 milliseconds
  > S2 forward  :     4716.30 milliseconds 	 residue  : 3 (80.000000 %)
elapsed time  15:   168433.53 milliseconds
  < S2 backward :     4661.32 milliseconds
  < S1 backward :      468.71 milliseconds 	 residue  : 3 (0.000000 %)
elapsed time  16:   173565.56 milliseconds
END

total elapsed time: 173565.562500 milliseconds
check number of shingles in S1: 3 == 3
check number of shingles in S2: 3 == 3

residue cross-pairing S1 x S2
=============================
the trigger substrings T are in S2
note: S1 != S2

number of start substrings:
- T trigger  : 3
- X collision: 3

- T start indices:
  0 330038357
  1 500000000
  2 673390940
- X start indices:
  0 11552071
  1 50000000
  2 51941278

repeated substring pairs (T,X) of length 7
-------------------------------------------
=> 0 2 330038357 51941278
=> 1 1 500000000 50000000
=> 2 0 673390940 11552071

number of substring repetitions: 3

The test pair is still in the sieve!

*/

#include <time.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <random>
#include <cstring>
#include <chrono>
using namespace std;
using Time = std::chrono::time_point<std::chrono::high_resolution_clock>;

#define SetBit(A,k)     ( A[(k/64)] |=  (1ULL << (k%64)) )
#define ClearBit(A,k)   ( A[(k/64)] &= ~(1ULL << (k%64)) )
#define TestBit(A,k)    ( A[(k/64)] &   (1ULL << (k%64)) )

uint64_t blq();			// returns (B^L) % M
void rcp_generator(std::mt19937& mt_rand, uint8_t p[]);
void string_generator(std::mt19937& mt_rand, uint8_t s[], uint64_t n);
uint64_t prime_performance_test(const uint8_t s[], const uint64_t ns, const uint8_t p[]);
float cross_sieve_movements(const float xr);

void forward_S1 (const uint8_t s[], const uint64_t ns, const uint64_t t[], uint64_t z[], const uint8_t p[]);
bool forward_S2 (const uint8_t s[], const uint64_t ns, uint64_t &nr, uint64_t x[], const uint64_t z[], const uint8_t p[]);
void backward_S2(const uint8_t s[], const uint64_t ns, const uint64_t t[], uint64_t z[], const uint8_t p[]);
bool backward_S1(const uint8_t s[], const uint64_t ns, uint64_t &nr, uint64_t x[], const uint64_t z[], const uint8_t p[]);

uint64_t start_extractor(uint64_t w[], const uint64_t ns, uint64_t ss_ind[]);
void residue_pairing(const uint8_t s2t[], const uint64_t n2t, uint64_t t2[], const uint8_t s1x[], const uint64_t n1x, uint64_t x1[]);

#define SS_DIM  300			// max. number of start shingles

#define B 257				// first prime  >  256

/*
//  20MB
#define N1  20000000ULL
#define N2  20000000ULL
#define M   20000003ULL
#define L   6
#define C   4295135ULL

//  30MB
#define N1  30000000ULL
#define N2  30000000ULL
#define M   30000001ULL
#define L   6
#define C   27911090ULL

//  40MB
#define N1  40000000ULL
#define N2  40000000ULL
#define M   40000003ULL
#define L   6
#define C   25905392ULL

//  50MB
#define N1  50000000ULL
#define N2  50000000ULL
#define M   50000017ULL
#define L   6
#define C   9549171ULL

//  60MB
#define N1  60000000ULL
#define N2  60000000ULL
#define M   60000023ULL
#define L   6
#define C   17063255ULL

//  70MB
#define N1  70000000ULL
#define N2  70000000ULL
#define M   70000027ULL
#define L   6
#define C   36377223ULL

//  80MB
#define N1  80000000ULL
#define N2  80000000ULL
#define M   80000027ULL
#define L   6
#define C   70269533ULL

//  90MB
#define N1  90000000ULL
#define N2  90000000ULL
#define M   90000049ULL
#define L   6
#define C   30641267ULL

// 100MB
#define N1  100000000ULL
#define N2  100000000ULL
#define M    99999989ULL
#define L   6
#define C   39210697ULL

// 110MB
#define N1  110000000ULL
#define N2  110000000ULL
#define M   110000017ULL
#define L   6
#define C   12985424ULL

// 120MB
#define N1  120000000ULL
#define N2  120000000ULL
#define M   120000007ULL
#define L   6
#define C   110707676ULL

// 130MB
#define N1  130000000ULL
#define N2  130000000ULL
#define M   130000027ULL
#define L   6
#define C   67671877ULL

// 140MB
#define N1  140000000ULL
#define N2  140000000ULL
#define M   139999991ULL
#define L   6
#define C   26038729ULL

// 150MB
#define N1  150000000ULL
#define N2  150000000ULL
#define M   150000029ULL
#define L   6
#define C   101809230ULL

*/

/*
Sample ( 1/10 and 1 Giga Byte)
==============================
#define B   257
*/
#define N1   100000000ULL	// 100 MB
#define N2  1000000000ULL	//   1 GB
#define M     99999989ULL
#define L   7
#define C     77150229ULL

#define M_64  ((M +1)/64 + 1)
#define N1_64 ((N1+1)/64 + 1)
#define N2_64 ((N2+1)/64 + 1)

Time start_timer() {
    return std::chrono::high_resolution_clock::now();
}

double get_elapsed_time(Time start) {
    Time end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d = end - start;
    std::chrono::microseconds us = std::chrono::duration_cast<std::chrono::microseconds>(d);
    return us.count() / 1000.0f;
}

//	*********************************************************************************************************************************************
//	*********************************************************************************************************************************************

int main() {
	uint64_t ss_ind1;	// start indices of the test pair
	uint64_t ss_ind2;
	bool forward_flag;	// ´true´ if end is reached (no unique signatures)
	bool backward_flag;
	uint64_t nr1;		// number of remaining shingles in the sieve (residue)
	uint64_t nr2;
	uint8_t  p[256];	// random cyclic permutations
	Time     start_time;			// start   time of timer
    uint64_t bit_count; 			// counter of bits that are set
    // concatenated length with L-1 overlaps (at end S1 and begin S2):
    // NC = (N1 + N2 - L + 1)

    const float pb= 1.0 / 256.0;	// byte probability
	const float ps= pow(pb, L);		// substring probability: p= pb^L (product over the substring)
    const float xr					// expected cross-repetitions S1 x S2, (without the artificially introduced repetitions)
			= (0.5 * (((N1+N2-L+1)-L)*((N1+N2-L+1)-L+1) - ((N1-L)*(N1-L+1)) - ((N2-L)*(N2-L+1))) * ps);
	const float xm					// expected number of required sieve movements
			= cross_sieve_movements(xr);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// M: Modulus; B: Base
// MAX        2MB                   255M           255			max= 2MB + 255   <=   2^64 - 1    ->   MB <= 2^63 - 2^8
// hash =  (hash + M) * B   -  C * p[s[j]]   +   p[s[j+L]];
// MIN         MB                     0              0          min=  MB - 255M  >=   0   ->   B >= 255
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	printf("\nparameter check \n");
	if (B < 255)  exit (1);
	if ((log2(M) + log2(B)) > 63) exit (2);
	if ((M * B) > (1ULL<<63) - (1ULL<<8)) exit (3);
	if (C != blq()) {
		printf("#define C= (B ^ L) %% M : %llu !\n", blq());
		exit (4);
	}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	// allocate s [bytes]; u1, v2 and z [bits]
	// =======================================
	printf("\nstorage allocation \n");
	// random byte string
	uint8_t *s;
	s= (uint8_t *)malloc(N1+N2+1);
	if (s == NULL) exit (99);
	// current shingle state (trigger / collision)
	uint64_t *u1;
	u1= (uint64_t *)malloc(8*N1_64);
	if (u1 == NULL) exit (98);
	uint64_t *v2;
	v2= (uint64_t *)malloc(8*N2_64);
	if (v2 == NULL) exit (97);
	// current incident vector (unique / colliding hash value)
	uint64_t *z;
	z= (uint64_t *)malloc(8*M_64);
	if (z == NULL) exit (96);

	// concatenate s1,s2 with L-1 overlaps (see Staged Sieve)
	// ===================================
	uint8_t *s1= s;
	uint8_t *s2= s+N1-(L-1);

	printf("\n");
	printf("Repeated Substrings Cross-Sieve S1 x S2 \n");
	printf("======================================= \n");
	printf("base               B: %d   \n", B);
	printf("modulus            M: %llu \n", M);
	printf("substrings length  L: %d   \n", L);
	printf("(B ^ L) %% M        C: %llu \n", C);
	printf("string S1 length  N1: %llu  \n", N1);
	printf("string S2 length  N2: %llu  \n", N2);
    printf("E(S1 x S2)        Ec: %5.2f  expected cross-repetitions between S1 and S2 \n", xr);
	printf("E(movements)      Em: %5.2f  expected number of required sieve movements \n",  xm );
	fflush(stdout);

	printf("\n");
	printf("start generation of input string \n");
	fflush(stdout);
	time_t  cur_time = time(NULL);		// current time
	mt19937 mt_rand(time(&cur_time));	// random number initialization

	// generate input string
	// ---------------------
	string_generator(mt_rand, s, N1+N2);
	// print the first and last element of the generated input string
	printf("first, last bytes of input string: %d %d \n", s[0], s[N1+N2-1]);
	// check overlap between S1 and S2
	printf("overlapping end   of S1: %d ... %d \n", s1[N1-(L-1)], s1[N1-1]);
	printf("overlapping begin of S2: %d ... %d \n", s2[0], s2[L-2]);
	fflush(stdout);

	// generate test repetition: copy the substring from position ss_ind1 to positions: ss_ind2
	ss_ind1= N1 / 2;	// arbitrary choice of positions ss_ind1 and ss_ind2
	ss_ind2= N2 / 2;
	// copy/paste the test duplicate from S1[ss_ind1] to S2[ss_ind2]
	for (int k= 0; k < L; k++ ) {
		s2[k + ss_ind2]=  s1[k + ss_ind1];
	}
	printf("test repetitions indices: %llu %llu \n", ss_ind1, ss_ind2);
	fflush(stdout);

	// TEST: 0.0047 microseconds is a good result (hardware dependent)
    rcp_generator(mt_rand, p);
	prime_performance_test(s1, N1, p);
	prime_performance_test(s2, N2, p);

	// start total elapsed time
    start_time = start_timer();

    // initialize shingle bit-vectors
    // ==============================
    // first forward: all S1 shingles are T-shingles
	memset(u1, -1, 8*N1_64);
    // first forward: all S2 shingles are X-shingles
	memset(v2, -1, 8*N2_64);

    // initialize residue
    // ==================
    // start with all shingles in the sieve
	nr1= N1 - L + 1;
	nr2= N2 - L + 1;

	// cross-sieve the shingles
	// ========================
	printf("\n");
	for (int kk= 1, i= 0; kk < 100; kk++) {

		printf("ROUND %d \n", kk);
		fflush(stdout);
		// shuffle
	    rcp_generator(mt_rand, p);

		// move forward from left to right
		// -------------------------------
		// reset incidents
		memset(z, 0, 8*M_64);
		forward_S1 (s1, N1, u1, z, p);
		forward_flag= forward_S2 (s2, N2, nr2, v2, z, p);
	    printf("elapsed time %3d:  %10.2f milliseconds \n", ++i, get_elapsed_time(start_time));
	    if (forward_flag) break;

		// is it worth to reverse the flow?
	    if (nr2 > 1.1*M) continue;

		// move backward from right to left
		// --------------------------------
		// reset incidents
		memset(z, 0, 8*M_64);
		backward_S2 (s2, N2, v2, z, p);
		backward_flag= backward_S1 (s1, N1, nr1, u1, z, p);
	    printf("elapsed time %3d:  %10.2f milliseconds \n", ++i, get_elapsed_time(start_time));
	    if (backward_flag) break;
	}
	printf("END \n\n");

	// get total elapsed time
    printf("total elapsed time: %f milliseconds \n", get_elapsed_time(start_time));

    // check resulting number of shingles in the sieve
    // ===============================================
    bit_count= 0;
	for (uint64_t j= 0; j < N1-L+1; j++) {
		if (TestBit(u1, j)) bit_count++;
	}
	printf("check number of shingles in S1: %llu == %llu \n", bit_count, nr1);
    bit_count= 0;
	for (uint64_t j= 0; j < N2-L+1; j++) {
		if (TestBit(v2, j)) bit_count++;
	}
	printf("check number of shingles in S2: %llu == %llu \n", bit_count, nr2);

	// residue cross-pairing S1 x S2: match substrings marked by u1 / v2
	// =============================
	printf("\n");
	printf("residue cross-pairing S1 x S2 \n");
	printf("============================= \n");
	if (forward_flag)  {
		printf("the trigger substrings T are in S1\n");
		residue_pairing(s1, N1, u1, s2, N2, v2);
	}
	if (backward_flag) {
		printf("the trigger substrings T are in S2\n");
		residue_pairing(s2, N2, v2, s1, N1, u1);
	}

    // check test pair
    // ===============
	printf("\n");
	if (TestBit(u1, ss_ind1) && TestBit(v2, ss_ind2)) printf("The test pair is still in the sieve! \n");
	else  printf("The test pair is NOT in the sieve?! \n");
}

//	*********************************************************************************************************************************************
//	*********************************************************************************************************************************************

uint64_t blq() {
	// returns C= (B ^ L) % M
	uint64_t result = 1ULL;
	for (int k = 0; k < L; k++) {
		result *= B;
		result %= M;
	}
	return result;
}

void rcp_generator(std::mt19937& mt_rand, uint8_t p[]) {
    std::uniform_int_distribution<uint8_t> dist(0, 255);
	// p: random cyclic permutations (see Sattolo / Fisher–Yates)
	int i, j;
	uint8_t rand;	// random number
	bool z[256];	// assigned random number
	for (j= 0; j < 256; j++) z[j]= false;
	for (j= 0; j < 256; j++) {
		rand= dist(mt_rand);
		while (z[rand]) {rand++;}
		z[rand]= true;
		p[j]= rand;
	}
	// test only
	for (j= 0; j < 256; j++ ) {
		for (i= j+1; i < 256; i++ ) {
			if (p[j] == p[i]) {
				printf("ERROR\n");
				exit (95);
			}
		}
	}
}

void string_generator(std::mt19937& mt_rand, uint8_t s[], const uint64_t ns) {
    std::uniform_int_distribution<uint8_t> dist(0, 255);
	// s : random string (uniform_int_distribution)
	// ns: length of string s
	for (uint64_t j= 0; j < ns; j++ ) {
		s[j]= dist(mt_rand);
	}
}

//	*********************************************************************************************************************************************

uint64_t prime_performance_test(const uint8_t s[], const uint64_t ns, const uint8_t p[]) {
	// the time it takes to roll through a string
	// depends on the selected prime ?!
	// try neighbor-primes of roughly the same size.

	// forward rolling on S1:
	// For each shingle compute its signature advancing from left to right
	// and record incidents in z
	// returns the last computed hash value
	// s : string S
	// ns: length of string S
	// L : shingle length
	// p : random cyclic permutations

	uint64_t hash;			// hash value of the current shingle
	uint64_t j;				// shingle begin location on the left	==> shingle index
	uint64_t uniques;		// number of unique shingles
	uint64_t avoid;			// avoid compiler optimization

	// start elapsed time
    Time start_time = start_timer();

	// initiate hash of the first, leftmost shingle
	j= 0;
	hash= 0;
	for (uint64_t k= j; k < j+L; k++) {
		hash= (hash * B + p[s[k]]) % M;
	}

	// roll hash
	// ---------
	avoid= uniques= 0;
	while (true) {
		// j:      shingle begin location on the left  ==> shingle index (left adjusted)
		// j+L-1:  shingle end   location on the right

		// avoid compiler optimization (elimination of the loop)
		// ---------------------------
		uniques++;
		avoid= hash + uniques;

		// last shingle start index: ns - L
		if (j >= ns - L) break;

		// update hash for the next shingle, rolling forward
		hash= ((hash + M) * B   -  C * p[s[j]]   +   p[s[j+L]]) % M;
		j++;
	}

    printf("prime_performance_test relative hash time: %f microseconds \t %llu \t %llu \n", 1000.0 * get_elapsed_time(start_time) / (ns - L), uniques, avoid);
	fflush(stdout);
	return (hash);
}

float cross_sieve_movements(const float xr) {
	// expected number of required cross sieve movements
	// in function of the modulus M (parameter)
	// note: a sieving round consists of two movements: forward followed by backward
	// N1    string S1 length [bytes]
	// N2    string S2 length [bytes]
	// L     length of the repeated substrings [bytes]
	// xr    expected cross-repetitions S1 x S2
	// nr    current number of shingles in the sieve (residue)
	float nr1= N1-L+1;	// number of shingles in S1
	float nr2= N2-L+1;	// number of shingles in S2
	int   i;			// current number of sieve movements
	float pe;			// probability that a slot remains empty after mapping nr1 (forward) or nr2 (backward)
	float ne;			// number of unique shingles (those that will be eliminated)
	float xm;			// expected number movements

	// expected number of initial forward only movements
	pe= pow ((1.0 - 1.0/M), nr1);
	xm= (log(nr1) - log(nr2)) / log(1.0 - pe);
/*
	printf("\n");
	printf("Expected Number of Required Cross-Sieve Movements \n");
	printf("================================================= \n");
	printf ("L   = %15d    [bytes] \n", L);
	printf ("N1  = %15llu    [bytes] \n", N1);
	printf ("N2  = %15llu    [bytes] \n", N2);
	printf ("M   = %18.2f modulus \n", (float) M);
	printf ("    =>%18.2f expected cross-repetitions S1 x S2 \n", xr);
	printf ("    =>%18.2f expected initial forward only sieve movements \n", xm);
	printf ("\n");
*/
	i= 0;
	while (true) {

		// forward movement
		// ----------------
		// > S1 : probability that a slot remains empty after mapping nr1
		pe= pow ((1.0 - 1.0/M), nr1);
		// > S2 : number of shingles nr2 that find an empty slot
		ne= pe * nr2;

		if (nr2 < xr + ne) {
			// expected final movement
			xm= i + ((log(xr) - log(nr2)) / log(1.0 - pe));
			i++;
//			printf ("  > forward  : i= %3d pe= %f ne= %.2f nr2= %.2f\n", i, pe, ne, nr2-ne);
			break;
		}

		// eliminate unique shingles nr2 that do not occur in nr1
		nr2-= ne;
		i++;
//		printf ("  > forward  : i= %3d pe= %f ne= %.2f nr2= %.2f\n", i, pe, ne, nr2);

		// is it worth to reverse the flow?
		if (nr2 > 1.1*N1) continue;

		// backward movement
		// -----------------
		// < S2 : probability that a slot remains empty after mapping nr2
		pe= pow ((1.0 - 1.0/M), nr2);
		// < S1 : number of shingles nr1 that find an empty slot
		ne= pe * nr1;

		if (nr1 < xr + ne) {
			// expected final movements
			xm= i + ((log(xr) - log(nr1)) / log(1.0 - pe));
			i++;
//			printf ("  < backward : i= %3d pe= %f ne= %.2f nr1= %.2f\n", i, pe, ne, nr1-ne);
			break;
		}

		// eliminate unique shingles nr1 that do not occur in nr2
		nr1-= ne;
		i++;
//		printf ("  < backward : i= %3d pe= %f ne= %.2f nr1= %.2f\n", i, pe, ne, nr1);

	}
/*
	printf ("    =>%18.2f relative cost (RAM(z) * CPU / N1) \n", ((float)i * M) / N1 );
	printf ("    =>%18.2f expected number of required sieve movements \n", xm);
	printf ("\n");
*/
	return (xm);
}

//	*********************************************************************************************************************************************
//	*********************************************************************************************************************************************

void forward_S1(const uint8_t s[], const uint64_t ns, const uint64_t t[], uint64_t z[], const uint8_t p[]) {
	// forward rolling on S1:
	// For each shingle compute its signature advancing from left to right
	// and record incidents in z
	// s : string S1
	// ns: length of string S1
	// L : shingle length
	// t : T-shingles (trigger in S1)
	// z : incident vector
	// p : random cyclic permutations

	uint64_t hash;		// hash value of the current shingle
	uint64_t j;			// shingle begin location on the left	==> shingle index

	// start elapsed time
    Time start_time = start_timer();

	// compute hash of the first, leftmost shingle
	j= 0;
	hash= 0;
	for (uint64_t k= j; k < j+L; k++) {
		hash= (hash * B + p[s[k]]) % M;
	}

	// incident loop
	while (true) {
		// j:      shingle begin location on the left  ==> shingle index (left adjusted)
		// j+L-1:  shingle end   location on the right

		// skip if shingle is eliminated
		if (TestBit(t, j)) {
			// record incident
			SetBit(z, hash);
		}
		// last shingle start index: ns - L
		if (j >= ns - L) break;

		// update hash for the next shingle, rolling forward
		hash= ((hash + M) * B   -  C * p[s[j]]   +   p[s[j+L]]) % M;
		j++;
	}

    printf("  > S1 forward  :  %10.2f milliseconds \n", get_elapsed_time(start_time));
	fflush(stdout);
}

bool forward_S2(const uint8_t s[], const uint64_t ns, uint64_t &nr, uint64_t x[], const uint64_t z[], const uint8_t p[]) {
	// forward rolling on S2:
	// For each shingle compute its signature advancing from left to right
	// and eliminate shingles with unique signatures;
	// returns ´true´ if no unique shingles
	// s : string S2
	// ns: length of string S2
	// nr: number of shingles currently in S2 (residue)
	// L : shingle length
	// x : X-shingles (collisions in S2)
	// z : incident vector
	// p : random cyclic permutations

	uint64_t hash;		// hash value of the current shingle
	uint64_t j;			// shingle begin location on the left	==> shingle index
	uint64_t uniques;	// number of unique shingles
	if (nr == 0) return(true);

	// start elapsed time
    Time start_time = start_timer();

    // compute hash of the first, leftmost shingle
	j= 0;
	hash= 0;
	for (uint64_t k= j; k < j+L; k++) {
		hash= (hash * B + p[s[k]]) % M;
	}

	// incident loop
	uniques= 0;
	while (true) {
		// j:      shingle begin location on the left  ==> shingle index (left adjusted)
		// j+L-1:  shingle end   location on the right

		// skip if shingle is eliminated
		if (TestBit(x, j) && !TestBit(z, hash))  {
			// shingle j of S2 has a unique signature relative to S1
			ClearBit(x, j);
			uniques++;
		}

		// last shingle start index: ns - L
		if (j >= ns - L) break;

		// update hash for the next shingle, rolling forward
		hash= ((hash + M) * B   -  C * p[s[j]]   +   p[s[j+L]]) % M;
		j++;
	}

    printf("  > S2 forward  :  %10.2f milliseconds \t residue  : %llu (%f %%) \n", get_elapsed_time(start_time), nr - uniques, (100.0 * uniques) / nr);
	fflush(stdout);
	nr= nr - uniques;
	return (uniques == 0);
}

//	*********************************************************************************************************************************************

void backward_S2(const uint8_t s[], const uint64_t ns, const uint64_t t[], uint64_t z[], const uint8_t p[]) {
	// backward rolling on S2:
	// For each shingle compute its signature advancing from right to left
	// and record incidents in z
	// s : string S2
	// ns: length of string S2
	// L : shingle length
	// t : T-shingles (trigger)
	// z : incident vector
	// p : random cyclic permutations

	uint64_t hash;		// hash value of the current shingle
	uint64_t i;			// shingle begin location on the right
	uint64_t j;			// shingle end   location on the left	==> shingle index

	// start elapsed time
    Time start_time = start_timer();

	// compute hash of the first, rightmost shingle
	i= ns - 1;
	hash= 0;
	for (uint64_t k= i; k > i-L; k--) {
		hash= (hash * B + p[s[k]]) % M;
	}

	// incident loop
	while (true) {
		// i:		   shingle begin location on the right
		j= i-L+1;	// shingle end   location on the left      ==> shingle index (left adjusted)

		// skip if shingle is eliminated
		if (TestBit(t, j)) {
			// record incident
			SetBit(z, hash);
		}

		// last shingle start index: L - 1
		if (i <= L-1) break;

		// update hash for the next shingle, rolling backward
		hash= ((hash + M) * B   -  C * p[s[i]]   +   p[s[i-L]]) % M;
		i--;
	}

    printf("  < S2 backward :  %10.2f milliseconds \n", get_elapsed_time(start_time));
	fflush(stdout);
}

bool backward_S1(const uint8_t s[], const uint64_t ns, uint64_t &nr, uint64_t x[], const uint64_t z[], const uint8_t p[]) {
	// backward rolling on S1:
	// For each shingle compute its signature advancing from right to left
	// and eliminate shingles with unique signatures;
	// returns ´true´ if no unique shingles
	// s : string S1
	// ns: length of string S1
	// nr: number of shingles currently in S1 (residue)
	// L : shingle length
	// x : X-shingles (collisions in S1)
	// z : incident vector
	// p : random cyclic permutations

	uint64_t hash;		// hash value of the current shingle
	uint64_t i;			// shingle begin location on the right
	uint64_t j;			// shingle end   location on the left	==> shingle index
	uint64_t uniques;	// number of unique shingles
	if (nr == 0) return(true);

	// start elapsed time
    Time start_time = start_timer();

	// compute hash of the first, rightmost shingle
	i= ns - 1;
	hash= 0;
	for (uint64_t k= i; k > i-L; k--) {
		hash= (hash * B + p[s[k]]) % M;
	}

	// incident loop
	uniques= 0;
	while (true) {
		// i:		   shingle begin location on the right
		j= i-L+1;	// shingle end   location on the left      ==> shingle index (left adjusted)

		// skip if shingle is eliminated
		if (TestBit(x, j) && !TestBit(z, hash))  {
			// shingle j of S1 has a unique signature relative to S2
			ClearBit(x, j);
			uniques++;
		}

		// last shingle start index: L - 1
		if (i <= L-1) break;

		// update hash for the next shingle, rolling backward
		hash= ((hash + M) * B   -  C * p[s[i]]   +   p[s[i-L]]) % M;
		i--;
	}

    printf("  < S1 backward :  %10.2f milliseconds \t residue  : %llu (%f %%) \n", get_elapsed_time(start_time), nr - uniques, (100.0 * uniques) / nr);
	fflush(stdout);
	nr= nr - uniques;
	return (uniques == 0);
}

//	*********************************************************************************************************************************************

uint64_t start_extractor(uint64_t w[], const uint64_t ns, uint64_t ss_ind[]) {
	// extract the start shingles of all contiguous shingle sequences
	// returns the number of detected start shingles
	// w : shingles (X / T)
	// ns: length of associated string (S1 / S2)
	// L : shingle length
	// ss_ind[i]:  start shingle indices
	const uint64_t w_len= ns - L + 1;

	uint64_t ss_count;	// number of start shingle
	uint64_t start;		// start  of adjacent shingles
	uint64_t stop;		// stop   of adjacent shingles
	uint64_t j;			// shingle index

	j= 0;
	ss_count= 0;
	while (true) {

		// get start/stop of adjacent X-shingles
		start= stop= 0;
		// w[w_len]= true;
		SetBit(w, w_len);
		// while (!w[j++]);
		while (!TestBit(w, j)) {j++;}
		// if(--j >= w_len) return(ss_count);
		if(j >= w_len) return(ss_count);
		start= j;
		// w[w_len]= false;
		ClearBit(w, w_len);
		// while (w[j++]);
		while (TestBit(w, j)) {j++;}
		// stop= j - 1;
		stop= j;
		// number of adjacent shingles = stop - start;

		// record first SS_DIM start shingles
		if (ss_count < SS_DIM) ss_ind[ss_count]= start;
		else return(ss_count);
		ss_count++;

		if (j >= w_len) return(ss_count);
	}
}

void residue_pairing(const uint8_t st[], const uint64_t nt, uint64_t t[], const uint8_t sx[], const uint64_t nx, uint64_t x[]) {
	// st: string T
	// nt: length of string T
	// t : T-shingles (trigger substrings are in T)
	// sx: string X
	// nx: length of string X
	// x : X-shingles (collision substrings are in X)

	uint64_t x_ind[SS_DIM];		// indices of the X start shingles
	uint64_t t_ind[SS_DIM];		// indices of the T start shingles
	int x_count;				// number  of start shingles in X
	int t_count;				// number  of start shingles in T
	int rep_count;				// number  of repeated shingles

	// note:
	if (sx != st) printf("note: S1 != S2 \n\n");
	else printf("note: S1 == S2 \n\n");

	// start shingle extractor
	// -----------------------
	printf("number of start substrings: \n");
	t_count= start_extractor(t, nt, t_ind);
	printf("- T trigger  : %d \n", t_count);
	x_count= start_extractor(x, nx, x_ind);
	printf("- X collision: %d \n", x_count);
	printf("\n");
	fflush(stdout);

	// shingle pairing
	// ---------------
	if (t_count >= SS_DIM || x_count >= SS_DIM) {
		printf("\n");
		printf("max. number of start substrings exceeded: SS_DIM= %d\n", SS_DIM);
		printf("too many substring repetitions \n");
	} else {
		printf("- T start indices: \n");
		for (int k= 0; k < t_count; k++) printf("  %d %llu \n", k, t_ind[k]);
		printf("- X start indices: \n");
		for (int k= 0; k < x_count; k++) printf("  %d %llu \n", k, x_ind[k]);

		printf("\n");
		printf("repeated substring pairs (T,X) of length %d \n", L);
		printf("------------------------------------------- \n");
		rep_count= 0;
		// T-shingle sequence
		for (int i= 0; i < t_count; i++) {
			// begin and end indices of the shingle sequence i
			uint64_t cur_t_ind= t_ind[i];
			// loop over all shingles in sequence i
			while (TestBit(t, cur_t_ind)) {
				// X-shingle sequence
				for (int j= 0; j < x_count; j++) {
					// begin and end indices of the shingle sequence j
					uint64_t cur_x_ind= x_ind[j];
					// loop over all shingles in sequence j
					while (TestBit(x, cur_x_ind)) {
						if (sx != st || cur_t_ind != cur_x_ind) {
							// check the substring pair byte by byte
							int k;
							for (k= 0; k < L ; k++) {
								if (st[cur_t_ind + k] != sx[cur_x_ind + k]) break;
							}
							if (k == L) {
								rep_count++;
								printf("=> %d %d %llu %llu \n", i, j, cur_t_ind, cur_x_ind);
							}
						}
						// X: begin and end index of the next prefix in cluster j
						if (++cur_x_ind + L - 1 >= nx) break;
					}
				}
				// T: begin and end index of the next prefix in cluster i
				if (++cur_t_ind + L - 1 >= nt) break;
			}
		}
		printf("\n");
		printf("number of substring repetitions: %d \n", rep_count);
	}
}

//	*********************************************************************************************************************************************
