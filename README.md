
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
- (1)&nbsp;	**Cross-Repetitions	S1 x S2 :** 	repetitions across two strings <br/>
at least one instance of a repetition appears in each string S1 and S2. <br/>
- (2)&nbsp;	**Self-Repetitions 	S2 :** 		repetitions within one and the same string S <br/>
- (3)&nbsp;	**Staged Repetitions	S12 x S2 :** 	repetitions within the 1st and across the 2nd string <br/>
repetitions are either within S1 or across S1 , S2 ; whereas repetitions within S2 are not taken into account. <br/>

Throughout, repetitions are exact matches of substrings (byte-by-byte); whether they repeat one or multiple times is irrelevant. <br/>

### Performance:
To sieve a 1 Giga Byte string for self-repetitions requires roughly a hash table of 1 Giga Bit and 135 seconds processing time on a MS Windows laptop (Inspiron 5748, Intel(R) Core(TM) i7-4510U CPU @ 2.00GHz; 8,00 GB RAM). Each time the problem size redoubles, both the memory for the hash table and the processing time are set to double as well. Thus, to upgrade our example to 2 Giga Bytes, would require 270 seconds and a hash table of 2 Giga Bits.

### Context
- **ISM** are the preferred radio frequency bands (434 / 868 / 912 MHz) used in smart home automation for remote-control and sensor data acquisition over the air
- **OOK**, On-Off Keying, is the modulation technique most widely found in low cost equipment. Information is transmitted by varying the duration of alternating HIGH- and LOW-signals. In general, these durations are restricted to a limited number of duration levels / categories (-> clusters)
- **Raw Data**, in form of signal duration sequences, is the common base level protocol of any OOK sender / receiver



### LICENSE
This project is released under [CC-BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).<br/>
The licensing TLDR is: You are free to use, copy, distribute and transmit this Software for personal, non-commercial purposes, as long as you give attribution and share any modifications under the same license. Commercial or for-profit use requires a license. <br/>
For more details see the [LICENSE](https://github.com/ookraw/OOK-Raw-Data-Receiver/blob/master/LICENSE)

### Setup Guide
This repository contains all that is needed to setup your Arduino workspace:
- receiver.ino
- categorizer.cpp
- categorizer.h
- categorizer_lib.cpp
- recorder.cpp
- radio_lib.cpp
- radio_lib.h
- RFM69_lib.cpp
- RFM69_registers.h

### Required Hardware
- MCU:   Arduino compatible mini pro 3.3V (MEGA328P)
- Radio: Hope RFM69w
- recommended: USB extension cable (5m) with snap-on Ferrite cores 

The Radio - MCU connections are defined in radio_lib.cpp <br/>
The prototype used during the development of the project can be found on  https://sites.google.com/site/rfm69arduino

### Project Presentation
An introduction to the project is available on https://sites.google.com/view/ookraw
