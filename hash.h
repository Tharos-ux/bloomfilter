#ifndef HASH_H
#define HASH_H

#include <cstdint>
#include <iostream>

using namespace std;

/** Hash function that uses xor properties
 * @param x a 64-bits value to hash
 * @return hashed value
 */
uint64_t xorshift64(uint64_t x)
{
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	return x;
}

/** Generate multiple hashes of the same value using sequences of xorshifts.
 * @param x The value to hash
 * @param hashes An array already allocated to fill with hash values
 * @param nb_hashes The number of hash values needed
 * @param max_val The maximum number that a hashed value can be.
 */
void multihash(uint64_t x, uint64_t *hashes, uint64_t nb_hashes, uint64_t max_val)
{
	// Init 64 bits hashes
	hashes[0] = xorshift64(x);
	for (uint64_t i = 1; i < nb_hashes; i++)
		hashes[i] = xorshift64(hashes[i - 1]);

	for (uint64_t i = 0; i < nb_hashes; i++)
		hashes[i] %= max_val + 1;
}

/* Compile:
 *	g++ hash.cpp -o hashtest
 * Exec:
 *	./hashtest 147
 */
int test_hash(int argc, char *argv[])
{
	// Allocate hash values space
	uint64_t hashes[3];
	uint64_t bloom_size = 21354543;

	// Encoded kmer to hash (from command line)
	uint64_t kmer_value = atoi(argv[1]);
	// Hashing
	multihash(kmer_value, hashes, 3, bloom_size - 1);

	// Output the hash values
	cout << "hash list for " << kmer_value << endl;
	for (uint64_t i = 0; i < 3; i++)
		cout << hashes[i] << endl;

	return 0;
}

#endif