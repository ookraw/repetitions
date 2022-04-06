
Repeated Substrings Sieve
--------------------------

Designed and coded by Felix Baessler, felix.baessler@gmail.com

The sieves are inspired by the Karp-Rabin algorithm; in particular they use the Rabin fingerprint modulo M as rolling hash. A minimalistic hash table reduced to M one-bit slots is used to keep track of the fingerprints. <br/>

Sieving a given string S of length N, consists in successively reducing the residue r, defined by the set of n substrings of length L. <br/>

As illustrated by the example, the initial residue contains n = N - L + 1 items, and ideally in the end, the final residue will contain repeating items only. <br/> 

L is critical as it controls the resulting amount of residue: if it is too large, the sieve will be empty, if L is too narrow, an unmanageable quantity of repetitions will survive. <br/> 

In telecommunications L may be interpreted as a noise level, in the sense that the shorter and more frequent the repetitions the less information they provide. However, by suppressing too much noise, one risks to miss valuable information. <br/> 

**Example:** Consider the string S of length N= (4*3) + (3*7) = 33 : <br/>
&nbsp; &nbsp; &nbsp;  	S = “abc1234567def1234567ghi1234567jkl”. <br/>
With L = 5 and moving forward from left to right, the initial forward residue is <br/>
&nbsp; &nbsp; &nbsp;  	r = { “abc12”, “bc123”, “c1234”,  … , “567jk”, “67jkl” }. <br/>
After sieving, the resulting residue will (implicitly) contain for each repetition the substrings: <br/>
&nbsp; &nbsp; &nbsp;    “12345”, “23456” and “34567”, <br/>
with the leftmost items representing the prefixes of the repeated substrings. <br/>
Similarly, the backward residue is obtained from the reversed S. <br/>

In the included documentation, the overlapping residue items, illustrated by the above example, are also called "shingles". Using this terminology the method can be stated as follows: <br/>

1.&nbsp;	**sieve** unique shingles: eliminate those shingles that produce unique fingerprints <br/>
2.&nbsp;	**extract** start shingles: identify those shingles in the resulting residue that mark beginnings of repetitions (prefixes) <br/>
3.&nbsp;	**pair** the start shingles: find the matching start shingles in the resulting residue <br/>


### Features:

The project demonstrates the principles of substring sieves operating on three types of repetitions: <br/>
- (1)&nbsp;	**Cross-Repetitions	S<sub>1</sub> x S<sub>2</sub> :** 	repetitions across two strings <br/>
at least one instance of a repetition appears in each string S<sub>1</sub> and S<sub>2</sub>. <br/>
- (2)&nbsp;	**Self-Repetitions 	S<sup>2</sup> :** 		repetitions within one and the same string S <br/>
- (3)&nbsp;	**Staged Repetitions	S<sub>1</sub><sup>2</sup> x S<sub>2</sub> :** 	repetitions within the 1st and across the 2nd string <br/>
repetitions are either within S<sub>1</sub> or across S<sub>1</sub> , S<sub>2</sub> ; whereas repetitions within S<sub>2</sub> are not taken into account. <br/>

Throughout, repetitions are exact matches of substrings (byte-by-byte); whether they repeat one or multiple times is irrelevant. <br/>

### Performance:
To sieve a 1 Giga Byte string for self-repetitions requires roughly a hash table of 1 Giga Bit and 135 seconds processing time on a MS Windows laptop (Inspiron 5748, Intel(R) Core(TM) i7-4510U CPU @ 2.00GHz; 8,00 GB RAM). Each time the problem size redoubles, both the memory for the hash table and the processing time are set to double as well. Thus, to upgrade our example to 2 Giga Bytes, would require 270 seconds and a hash table of 2 Giga Bits.

### Context
Intuitively, identifying repetitions and establishing correlations between repetitions is fundamental in any learning process.  <br/>

In telecommunications we are interested in both, learning about:  <br/>

 - **self-repeated information** such as syncs, preambles, delimiters, headers or even entire message blocks, as for example bursts that repeat in the same data stream, and <br/> 
 - **cross-repeated information** of identical subsequences originating from correlated sources that are for example candidates for sharing the same protocol.<br/>

For a specific application it is always interesting, and for large data certainly advisable, to analyze the distribution of the input before starting the sieving process. On the other hand, in the case of small problems, a more pragmatic approach could do the job:
Start with a rather low noise threshold (L) and increase it on the fly until the amount of resulting residue becomes manageable for pairing.<br/> 

### Setup Guide
This repository contains three self-contained programs:
- cross_repetitions_sieve.cpp
- self_repetitions_sieve.cpp
- staged_sieve.cpp  (follows later) <br/>

No attempt has been made to “tune” the code. On the other hand some effort was made to help the compiler optimize the sieving time. For this reason it will be necessary to recompile the programs with new define statements for each problem. <br/>

The following compiler flags have been used:
-	-O3 -g3 -Wall         		: for optimization
-	-Wl,--stack,0xFFFFFF  	: for long arrays  <br/>

Note that you might have to upgrade minGW to x86_64 in order to obtain 64 bit executables.  <br/>

### Project Presentation
An analysis of the algorithms is available on https://sites.google.com/view/repsieve

### LICENSE
This project is released under [CC-BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).<br/>
The licensing TLDR is: You are free to use, copy, distribute and transmit this Software for personal, non-commercial purposes, as long as you give attribution and share any modifications under the same license. Commercial or for-profit use requires a license. <br/>
For more details see the [LICENSE](https://github.com/ookraw/OOK-Raw-Data-Receiver/blob/master/LICENSE)

