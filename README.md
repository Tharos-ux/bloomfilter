# C++ implementation of a Bloom filter
## Compilation instructions
In order to compile the code, you may use `g++ -o bf bloomfilter.cpp` with GNU C++ compiler.
Then, program may be executed with `./bf path_to_file.fasta k n nf r` (same order as described in original README).
## Parameters
Parameters are checked within the main function, and program will terminate if any parameter does not match specification.

+ **path_to_file.fasta** must link to an existing file
+ **k** is the kmer size. It should be in the interval [ 1 : 31 ]
+ **n** is the size of the bloom filter. It should be in the interval [ 1 : $2^{34}$ ]
+ **nf** is the number of successive hash functions we will apply on encoded kmers. It should be in the interval [ 1 : 64 ]
+ **r** is the number of randoms kmers we ask if they are in our data structure. Should be positive
## Known issues
+ If the file is not a true fasta file, program will fail to return coherent results