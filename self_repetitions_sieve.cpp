//============================================================================
// Name        : self_repetitions_sieve.cpp
// Author      : Felix Baessler
// Version     : 01.04.2022
// Copyright   : Felix Baessler, felix.baessler@gmail.com
// SEE TLDR; VERSION OF THE LICENSE: https://creativecommons.org/licenses/by-nc/4.0/legalcode
// SEE FULL LICENSE DETAILS HERE   : https://creativecommons.org/licenses/by-nc/4.0/
// Description : Repeated Substrings / self_repetitions_sieve
// Compilation flags:
// -O3 -g3 -Wall         : optimization
// -Wl,--stack,0xFFFFFF  : long arrays
// if required: upgrade minGW to x86_64 : 64 bit executable
//
// introduction available on: https://sites.google.com/view/repsieve
//============================================================================

/*
Repeated Substrings : Self-Repetitions Sieve
============================================

This program illustrates the principles of a sieve, that solves the following problem:

Find all comparatively short repetitions hidden within a long string of bytes.

For the purpose of the experiment, the long input string S (~ Giga Bytes) is generated
by means of the standard MT19937 pseudo-random number generator. Also to make sure that
at least one repetition can be found, non-overlapping and non-adjacent test duplicates
(~ 7 Bytes) are copy/pasted within S.

The conjecture is: as long as the (minimal) string length of the repetitions we are looking for,
are longer than the length of the majority of the inherent random repetitions, the sieve can be
competitive:
 - in terms of space, as only three bit-vectors of the length of S are required
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
Both state vectors, t and x, are bit-vectors of the same length as S.
Along with the input string, they are sequentially accessed during the sieving process.
Trigger Vector:
 -	t[j] = true : the shingle j is a T-shingle
 -	t[j] = false: the shingle j is NOT a T-shingle
At the beginning of the sieving process, t is true for all shingles.
Collision Vector:
 -	x[j] = true : the shingle j is an X-shingle, it MAY repeat with a T-shingle
 -	x[j] = false: the shingle j is NOT an X-shingle
At the beginning of the sieving process, x is false for all shingles.

Reversing the processing direction is achieved by switching between t and x.

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
the same length as S.

State Space
-----------
The global state space of the sieve presents itself as follows:
 -	t[j] = false &&  x[j] = false :  shingle j is identified as an E-shingle
 -	t[j] = true  &&  x[j] = false :  shingle j is identified as a  T-shingle
 -	t[j] = false &&  x[j] = true  :  shingle j is identified as an X-shingle
 -	t[j] = true  &&  x[j] = true  :  will never occur

State Check
-----------
During the sieving process, it is always possible to check whether a test pair is still in the sieve:
if ((t[j1] && x[j2]) || (x[j1] && t[j2])) printf("The test pair is still in the sieve! \n");
where j1, j2 denote the respective indices in S.

Sieve Filter
------------
	if (TestBit(t, j)) {
		// shingle j is a T-shingle
		if (!TestBit(z, hash))  {
			SetBit(z, hash);
		} else {
			// collisions++;
			// state change: T -> X
			SetBit(x, j);
			ClearBit(t, j);
		}
	} else {
		if (TestBit(x, j)) {
			// shingle j is an X-shingle
			if (!TestBit(z, hash))  {
				// shingle j has a unique signature
				uniques++;
				ClearBit(x, j);
			} else {
				// collisions++;
			}
		}
		// shingle j is an E-shingle (elimination)
	}


Program Usage
=============
The program parameters N, M and L are set by define statements:

 -  N is the length of S, i.e. the number of random bytes
    generated at the beginning of the program.

 -	M is the modulus of the Rabin fingerprint;
	it defines the length of the signature incident vector

 -	L is the length of the shingles;
	it is also used to set the string length of the test pairs.

Two additional define statements concern:

 -  B = 257, the base of the Rabin fingerprint (first prime > 256)
 -  C = (B ^ L) % M, a derived constant intended to help compiler optimization

For a given N, the performance of the sieve depends on the choice of M and L.

Prime numbers can be obtained from: https://www.walter-fendt.de/html5/mde/primenumbers_de.htm

Sample OUTPUT ( 1 Giga Byte)
============================
#define B   257
#define N   1000000000ULL
#define M   1000000007ULL
#define L   7
#define C   13163680ULL

parameter check
storage allocation

Repeated Substrings Sieve
=========================
base              B: 257
modulus           M: 1000000007
substrings length L: 7
(B ^ L) % M       C: 13163680
string S length   N: 1000000000
E(S x S)         Es:  6.94  expected repetitions in S
E(movements)     Em:  7.00  expected number of required sieve movements

start generation of input string
first, last bytes of input string: 17 145
test repetitions indices: 0 333333333 666666666 999999993
prime_performance_test relative hash time: 0.004725 microseconds 	 999999994 	 1114245839

ROUND 1
  > S  forward  :    30459.57 milliseconds 	 residue  : 999999994 (0.000000 %)
elapsed time   1:    30608.47 milliseconds
  < S  backward :    31038.22 milliseconds 	 residue  : 512846627 (48.715337 %)
elapsed time   2:    61669.66 milliseconds
ROUND 2
  > S  forward  :    23076.78 milliseconds 	 residue  : 251211615 (51.016229 %)
elapsed time   3:    84769.43 milliseconds
  < S  backward :    15981.81 milliseconds 	 residue  : 75427715 (69.974432 %)
elapsed time   4:   100775.27 milliseconds
ROUND 3
  > S  forward  :     9274.71 milliseconds 	 residue  : 12344141 (83.634476 %)
elapsed time   5:   110071.96 milliseconds
  < S  backward :     5835.67 milliseconds 	 residue  : 694673 (94.372448 %)
elapsed time   6:   115929.60 milliseconds
ROUND 4
  > S  forward  :     4806.21 milliseconds 	 residue  : 7161 (98.969155 %)
elapsed time   7:   120757.81 milliseconds
  < S  backward :     4712.31 milliseconds 	 residue  : 19 (99.734674 %)
elapsed time   8:   125492.12 milliseconds
ROUND 5
  > S  forward  :     4720.26 milliseconds 	 residue  : 16 (15.789474 %)
elapsed time   9:   130234.39 milliseconds
  < S  backward :     4791.22 milliseconds 	 residue  : 16 (0.000000 %)
elapsed time  10:   135049.62 milliseconds
END
total elapsed time: 135049.625000 milliseconds
check number of shingles: 16 == 16

residue self-pairing S x S
==========================
the trigger substrings T are in v
note: S1 == S2

number of start substrings:
- T trigger  : 7
- X collision: 9

- T start indices:
  0 232219364
  1 666787182
  2 763309560
  3 857807308
  4 881113915
  5 930664034
  6 999999993
- X start indices:
  0 0
  1 192235252
  2 228714536
  3 333333333
  4 491005087
  5 546648114
  6 575380659
  7 666666666
  8 840796034

repeated substring pairs (T,X) of length 7
-------------------------------------------
=> 0 2 232219364 228714536
=> 1 6 666787182 575380659
=> 2 5 763309560 546648114
=> 3 8 857807308 840796034
=> 4 1 881113915 192235252
=> 5 4 930664034 491005087
=> 6 0 999999993 0
=> 6 3 999999993 333333333
=> 6 7 999999993 666666666

number of substring repetitions: 9

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
void string_generator(std::mt19937& mt_rand, uint8_t s[], const uint64_t ns);
uint64_t prime_performance_test(const uint8_t s[], const uint64_t ns, const uint8_t p[]);
int sieve_movements(const float xr);

bool forward (const uint8_t s[], const uint64_t ns, uint64_t &nr, uint64_t t[], uint64_t x[], uint64_t z[], const uint8_t p[]);
bool backward(const uint8_t s[], const uint64_t ns, uint64_t &nr, uint64_t t[], uint64_t x[], uint64_t z[], const uint8_t p[]);

uint64_t start_extractor(uint64_t w[], const uint64_t ns, uint64_t ss_ind[]);
void residue_pairing(const uint8_t s1x[], const uint64_t n1x, uint64_t x[], const uint8_t s2t[], const uint64_t n2t, uint64_t t[]);

#define SS_DIM  300			// max. number of start shingles

#define B 257				// first prime  >  256
//#define B 8191			// first Mersenne prime  >  256

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

//#define N    N1

/*
Sample: 1 Giga Byte String
==========================
#define B   257
*/
#define N   1000000000ULL
#define M   1000000007ULL
#define L   7
#define C   13163680ULL

#define M_64 ((M+1)/64 + 1)
#define N_64 ((N+1)/64 + 1)

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
	uint64_t nr;		// number of remaining shingles in the sieve (residue)
	uint8_t  p[256];	// random cyclic permutations

	Time     start_time;			// start   time of timer
    uint64_t bit_count; 			// counter of bits that are set

    const float pb= 1.0 / 256.0;	// byte probability
	const float ps= pow(pb, L);		// substring probability: p= pb^L (product over the substring)
    const float xr					// expected repetitions S x S, (without the artificially introduced repetitions)
			= (0.5 * ((N-L)*(N-L+1)) * ps);
	const float xm					// expected number of required sieve movements
			= sieve_movements(xr);

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

	// allocate s [bytes]; u, v and z [bits]
	// =====================================
	printf("\nstorage allocation \n");
	// random byte string
	uint8_t *s;
	s= (uint8_t *)malloc(N+1);
	if (s == NULL) exit (99);
	// current shingle state (trigger / collision)
	uint64_t *u;
	u= (uint64_t *)malloc(8*N_64);
	if (u == NULL) exit (98);
	uint64_t *v;
	v= (uint64_t *)malloc(8*N_64);
	if (v == NULL) exit (97);
	// current incident vector (unique / colliding hash value)
	uint64_t *z;
	z= (uint64_t *)malloc(8*M_64);
	if (z == NULL) exit (96);

	printf("\n");
	printf("Repeated Substrings Sieve \n");
	printf("========================= \n");
	printf("base              B: %d   \n", B);
	printf("modulus           M: %llu \n", M);
	printf("substrings length L: %d   \n", L);
	printf("(B ^ L) %% M       C: %llu \n", C);
	printf("string S length   N: %llu  \n", N);
    printf("E(S x S)         Es: %5.2f  expected repetitions in S \n", xr);
	printf("E(movements)     Em: %5.2f  expected number of required sieve movements \n",  xm );
	fflush(stdout);

	printf("\n");
	printf("start generation of input string \n");
	fflush(stdout);
	time_t  cur_time = time(NULL);		// current time
	mt19937 mt_rand(time(&cur_time));	// random number initialization

	// generate input string
	// ---------------------
	string_generator(mt_rand, s, N);
	// print the first and last element of the generated input string
	printf("first, last bytes of input string: %d %d \n", s[0], s[N-1]);
	fflush(stdout);

	// generate test repetitions: copy the substring from position ss_ind1
	// to positions: 0, ss_ind2, N-L
	ss_ind1= N / 3;	// arbitrary choice of positions ss_ind1 and ss_ind2
	ss_ind2= 2 * ss_ind1;
	// copy/paste the test duplicate from ss_ind1 to ss_ind2
	for (int k= 0; k < L; k++ ) {
		s[k + ss_ind2]=  s[k + ss_ind1];
		//		/*
		// insert also a duplicate at the beginning and the end
		s[k]=  s[k + ss_ind1];
		s[k + N - L]=  s[k + ss_ind1];
		//		*/
	}
	printf("test repetitions indices: 0 %llu %llu %llu \n", ss_ind1, ss_ind2, N - L);
	fflush(stdout);

	// TEST: 0.0047 microseconds is a good result (hardware dependent)
    rcp_generator(mt_rand, p);
	prime_performance_test(s, N, p);

	// start total elapsed time
    start_time = start_timer();

    // initialize shingle bit-vectors
    // ==============================
    // first forward: all shingles are T-shingles
    memset(u, -1, 8*N_64);
    // first forward: there are no X-shingles
    memset(v,  0, 8*N_64);

    // initialize residue
    // ==================
    // start with all shingles in the sieve
	nr= N - L + 1;

	// sieve the shingles
	// ==================
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
		forward_flag= forward(s, N, nr, u, v, z, p);
	    printf("elapsed time %3d:  %10.2f milliseconds \n", ++i, get_elapsed_time(start_time));
	    if (kk > 1 && forward_flag) break;

		// move backward from right to left (switch u <-> v)
		// -------------------------------------------------
		// reset incidents
		memset(z, 0, 8*M_64);
		backward_flag= backward(s, N, nr, v, u, z, p);
	    printf("elapsed time %3d:  %10.2f milliseconds \n", ++i, get_elapsed_time(start_time));
		if (backward_flag) break;
	}
	printf("END \n");

	// get total elapsed time
    printf("total elapsed time: %f milliseconds \n", get_elapsed_time(start_time));

    // check resulting number of shingles in the sieve
    // ===============================================
    bit_count= 0;
	for (uint64_t j= 0; j < N-L+1; j++) {
		if (TestBit(u, j) || TestBit(v, j)) bit_count++;
	}
	printf("check number of shingles: %llu == %llu \n", bit_count, nr);

	// residue self-pairing S x S: match substrings marked by u / v
	// ==========================
	printf("\n");
	printf("residue self-pairing S x S \n");
	printf("========================== \n");
	if (forward_flag)  {
		printf("the trigger substrings T are in u\n");
		residue_pairing(s, N, u, s, N, v);
	}
	if (backward_flag) {
		printf("the trigger substrings T are in v\n");
		residue_pairing(s, N, v, s, N, u);
	}

    // check test pair
    // ===============
	printf("\n");
	if ((TestBit(u, ss_ind1) || TestBit(v, ss_ind1)) && (TestBit(u, ss_ind2) || TestBit(v, ss_ind2))) printf("The test pair is still in the sieve! \n");
	else  printf("The test pair is NOT in the sieve?! \n");
}

//	*********************************************************************************************************************************************

uint64_t blq() {
	// returns (B^L) % M
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

int sieve_movements(const float xr) {
	// expected number of required sieve movements
	// in function of the modulus M (parameter)
	// note: a sieving round consists of two movements: forward followed by backward
	// N     string length [bytes]
	// L     length of the repeated substrings [bytes]
	// xr    expected number of repetitions S x S

	float nu= N-L+1;	// number of trigger   T-shingles
	float nv= 0;		// number of collision X-shingles
	int   i;			// current number of sieve movements
	float pe;			// probability that a slot remains empty after mapping u
	float ne;			// number of unique shingles (those that are eliminated)
	float nt;			// tmp: number of trigger   T-shingles
	float nc;			// tmp: number of collision X-shingles
/*
	printf("\n");
	printf("Expected Number of Required Sieve Movements \n");
	printf("=========================================== \n");
	printf ("L   = %15d    [bytes] \n", L);
	printf ("N   = %15llu    [bytes] \n", N);
	printf ("M   = %15llu modulus \n", M);
	printf ("    =>%18.2f expected repetitions S x S \n", xr);
	printf ("\n");
*/
	i= 0;
	while (true) {

		// After mapping the fingerprints of u
		// -----------------------------------
		// probability that a slot remains empty
		pe= pow ((1.0 - 1.0/M), nu);
		// expected number of collision free mappings
		nt= M * (1.0 - pe);
		// expected number of colliding mappings
		nc= nu - nt;
		if (nc <= 0) break;

		// Probing with v
		// --------------
		ne= nv * pe;

		// Balance
		// -------
		nu= nt;  //	= nu - nc;
		nv= nv + nc - ne;

		i++;
//		printf (" i= %3d nt= %f nc= %f ne= %f nu= %f nv= %f r= %f \n", i, nt, nc, ne, nu, nv, nu+nv);
		if (nv < xr) break;

		// Switch between u and v
		nu= nv;
		nv= nt;

	}
/*
	printf ("    =>%15llu  relative cost (RAM(z) * CPU / N) \n", (i * M) / N );
	printf ("    =>%15d  expected number of required sieve movements \n", i);
	printf ("\n");
*/
	return (i);
}

//	*********************************************************************************************************************************************
//	*********************************************************************************************************************************************

bool forward(const uint8_t s[], const uint64_t ns, uint64_t &nr, uint64_t t[], uint64_t x[], uint64_t z[], const uint8_t p[]) {
	// forward rolling shingle signature filter
	// For each shingle compute its signature advancing from left to right
	// and eliminate the shingles with unique signature
	// returns ´true´ if end is reached (no unique signatures)
	// s : string S
	// ns: string length
	// nr: current number of shingles in the sieve (residue)
	// L : shingle length
	// t : T-shingles (trigger)    can only change to false
	// x : X-shingles (collisions) can change to true and false
	// z : incident vector         can only change to true
	// p : random cyclic permutations

	uint64_t hash;			// hash value of the current shingle
	uint64_t j;				// shingle begin location on the left	==> shingle index
	// uint64_t collisions;	// number of colliding shingles
	uint64_t uniques;		// number of unique shingles
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
	// -------------
	uniques= 0;
	// collisions= 0;
	while (true) {
		// j:      shingle begin location on the left  ==> shingle index (left adjusted)
		// j+L-1:  shingle end   location on the right

	    // filter (same code for forward and backward)
	    // ======
	    // u and v must be initialized as follows:
	    // memset(u, -1, 8*N_64);
	    // memset(v,  0, 8*N_64);
		if (TestBit(t, j)) {
			if (!TestBit(z, hash))  {
				SetBit(z, hash);
			} else {
				// collisions++;
				SetBit(x, j);
				ClearBit(t, j);
			}
		} else {
			if (TestBit(x, j)) {
				if (!TestBit(z, hash))  {
					// shingle j has a unique signature  -> to be eliminated from the sieve
					uniques++;
					ClearBit(x, j);
				} else {
					// collisions++;
				}
			}
		}
		// test: if (TestBit(t, j) && TestBit(x, j)) exit (12345);

		// last shingle start index: ns - L
		if (j >= ns - L) break;

		// update hash for the next shingle, rolling forward
		hash= ((hash + M) * B   -  C * p[s[j]]   +   p[s[j+L]]) % M;
		j++;
	}

    printf("  > S  forward  :  %10.2f milliseconds \t residue  : %llu (%f %%) \n", get_elapsed_time(start_time), nr - uniques, (100.0 * uniques) / nr);
    //	printf(" > collisions: %llu (%f %%) \n", collisions, (100.0*collisions)/M);
	fflush(stdout);
	nr= nr - uniques;
	return (uniques == 0);
}

bool backward(const uint8_t s[], const uint64_t ns, uint64_t &nr, uint64_t t[], uint64_t x[], uint64_t z[], const uint8_t p[]) {
	// backward rolling shingle signature filter
	// For each shingle compute its signature advancing from right to left
	// and eliminate the shingles with unique signature
	// returns ´true´ if end is reached (no unique signatures detected)
	// s : string S
	// ns: string length
	// nr: current number of shingles in the sieve (residue)
	// L : shingle length
	// t : T-shingles (trigger)    can only change to false
	// x : X-shingles (collisions) can change to true and false
	// z : incident vector         can only change to true
	// p : random cyclic permutations

	uint64_t hash;			// hash value of the current shingle
	uint64_t i;				// shingle begin location on the right
	uint64_t j;				// shingle end   location on the left	==> shingle index
	uint64_t uniques;		// number of unique shingles
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
	// -------------
	uniques= 0;
	while (true) {
		// i:		   shingle begin location on the right
		j= i-L+1;	// shingle end   location on the left      ==> shingle index (left adjusted)

	    // filter (same code for forward and backward)
	    // ======
		if (TestBit(t, j)) {
			if (!TestBit(z, hash))  {
				SetBit(z, hash);
			} else {
				SetBit(x, j);
				ClearBit(t, j);
			}
		} else {
			if (TestBit(x, j)) {
				if (!TestBit(z, hash))  {
					// shingle j has a unique signature  -> to be eliminated from the sieve
					uniques++;
					ClearBit(x, j);
				}
			}
		}
		// test: if (TestBit(t, j) && TestBit(x, j)) exit (12345);

		// last shingle start index: L - 1
		if (i <= L - 1) break;

		// update hash for the next shingle, rolling backward
		hash= ((hash + M) * B   -  C * p[s[i]]   +   p[s[i-L]]) % M;
		i--;
	}

    printf("  < S  backward :  %10.2f milliseconds \t residue  : %llu (%f %%) \n", get_elapsed_time(start_time), nr - uniques, (100.0 * uniques) / nr);
	fflush(stdout);
	nr= nr - uniques;
	return (uniques == 0);
}

//	*********************************************************************************************************************************************

uint64_t start_extractor(uint64_t w[], const uint64_t ns, uint64_t ss_ind[]) {
	// extract the start shingle from a sequence of contiguous shingles
	// returns number of detected start shingles
	// w : shingles (X or T)
	// ns: length of string (S)
	// L : shingle length
	// ss_ind[i]:  start shingle indices
	const uint64_t w_len= ns - L + 1;  // number of shingles

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
