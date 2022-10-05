#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <tuple>
#include <vector>
#include <bitset>
#include "hash.h"

using namespace std;

class BloomFilter
{
public:
    // n is not constant so no bitset
    uint64_t hash_func_number;
    uint64_t bloom_size;
    vector<bool> filter;
    /** Constructor for a BloomFilter object.
     * @param n size of filter in bits
     * @param nf number of successive hash functions to apply
     */
    BloomFilter(int n, int nf)
    {
        // number of successive hash functions
        hash_func_number = nf;
        // will store the size of the bloomfilter
        bloom_size = n;
        // filter itself
        filter = vector<bool>(n, false);
    }

    /** Adds a value to the filter
     * @param value_to_hash encoded kmer
     */
    void add_value(uint64_t value_to_hash)
    {
        // defining array for storing hash function results
        uint64_t *hashes = new uint64_t[bloom_size / 8 + 1];
        multihash(value_to_hash, hashes, hash_func_number, bloom_size - 1);
        for (uint64_t i_hash = 0; i_hash < hash_func_number; i_hash++)
        {
            filter[hashes[i_hash]] = true;
        }
        return;
    }

    /** Search for a kmer in the filter
     * A kmer is in if all bit which represents the kmer are
     * @param value_to_search encoded kmer
     * @return A boolean about the presence of the kmer
     */
    bool is_present(uint64_t value_to_search)
    {
        // hashing value we search for
        uint64_t *hashes = new uint64_t[bloom_size / 8 + 1];
        multihash(value_to_search, hashes, hash_func_number, bloom_size - 1);

        // iterating through hashed values
        for (uint64_t i_hash = 0; i_hash < hash_func_number; i_hash++)
        {
            if (!filter[hashes[i_hash]])
            {
                // if one value is not present we are sure kmer is not in our structure
                return false;
            }
        }

        return true;
    }
};

/** Reads char by char a fasta file.
 * @param filename Path to the file to read
 * @param position Current cursor position
 * @return A tuple containing next char and new cursor position
 */
std::tuple<char, int> read_fasta(string filename, long position)
{
    // Opening stream
    std::ifstream fasta_file(filename);
    if (fasta_file.is_open())
    {
        // We assume standard fasta here, so we skip full first line
        if (position == 0)
        {
            fasta_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        else
        {
            fasta_file.ignore(position - 2);
        }
        // Iterating while there's chars in file
        while (fasta_file)
        {
            char fasta_char;
            fasta_char = fasta_file.get();
            if (fasta_char == 'A' || fasta_char == 'T' || fasta_char == 'C' || fasta_char == 'G')
            {
                return std::make_tuple(fasta_char, fasta_file.tellg());
            }
        }
    }
    else
    {
        // we could not read the file
        std::cout << "Error while reading file.\n";
    }
    // if we can't get any acceptable next char we return empty char
    return std::make_tuple('\0', fasta_file.tellg());
}

/** Returns the next clean kmer from the file
 * @param filename Path to the file to read
 * @param k k-mer length
 * @param previous_position Starting point of the previous kmer
 * @return A tuple containing next kmer and cursor position of the start of the kmer
 */
std::tuple<string, long> next_kmer(string filename, int k, long previous_position)
{
    long cursor_previous_kmer = previous_position;
    long cursor;
    char current_char;
    string kmer;

    // getting starting point for next kmer
    tie(current_char, cursor) = read_fasta(filename, cursor_previous_kmer);
    cursor_previous_kmer = cursor;
    kmer += current_char;
    for (int i = 1; i < k; i++)
    {
        tie(current_char, cursor) = read_fasta(filename, cursor);
        if (current_char != '\0')
        {
            kmer += current_char;
        }
        else
        {
            return std::make_tuple("", cursor_previous_kmer);
        }
    }
    // we return the next kmer
    return std::make_tuple(kmer, cursor_previous_kmer);
}

/** Returns the reverse complement
 * @param kmer String reprensenting a kmer
 * @return The DNA reverse complement of the string
 */
string reverse_complement(string kmer)
{
    string reverse_kmer;
    int n = kmer.length();
    // We iterate backwards and substitute letters with conditions
    while (n--)
    {
        if (kmer[n] == 'A')
        {
            reverse_kmer += 'T';
        }
        else if (kmer[n] == 'T')
        {
            reverse_kmer += 'A';
        }
        else if (kmer[n] == 'C')
        {
            reverse_kmer += 'G';
        }
        else
        {
            reverse_kmer += 'C';
        }
    }
    // As strings are cleaned when loaded, we don't have to care for other letters
    return reverse_kmer;
}

/** Returns the lexicographically first chain between a string and its reverse complement
 * @param kmer the string to analyze
 * @return The selected string
 */
string choose_kmer(string kmer)
{
    string reverse_comp = reverse_complement(kmer);
    if (kmer < reverse_comp)
    {
        return kmer;
    }
    else
    {
        return reverse_comp;
    }
}

/** Returns the encoding of the letter
 * @param letter the letter to encode
 * @return The int representing the letter
 */
uint64_t encode_letter(char letter)
{
    if (letter == 'A')
    {
        return 0;
    }
    else if (letter == 'T')
    {
        return 3;
    }
    else if (letter == 'C')
    {
        return 1;
    }
    else
    {
        return 2;
    }
}

/** Return the encoded value of the kmer (between its own and its reverse complement)
 * @param kmer the string to analyze
 * @return the code of the kmer
 */
uint64_t encode_kmer(string kmer)
{
    uint64_t k_code = 0;
    kmer = choose_kmer(kmer);
    int i = kmer.length();
    // We iterate backwards and substitute letters with conditions
    while (i--)
    {
        k_code += encode_letter(kmer[i]) * (4 ^ i);
    }
    // In order to avoid specific case where we happen to have a value of 0
    // (meaning we have AAAAAAA...), we add 1 to all results
    return k_code + 1;
}

/** Returns the length of a text file
 * @param filename Path to the file to read
 * @return Number of chars in the file
 */
long get_file_length(string filename)
{
    long char_count;
    ifstream fasta_file(filename);
    fasta_file.seekg(0, std::ios_base::end);
    char_count = fasta_file.tellg();
    fasta_file.close();

    return char_count;
}

/** Generates a random kmer
 * @param ksize length of generated kmer
 * @return A string of size ksize composed of A, T, C, G
 */
string generate_kmer(int ksize)
{
    string kmer;
    char alpha[4] = {'A', 'T', 'C', 'G'};
    for (int i = 0; i < ksize; i++)
    {
        kmer += alpha[rand() % 4];
    }
    return kmer;
}

/** Generates a random kmer
 * @param filename potential path to a file
 * @param k candidate for kmer size
 * @param n candidate for bloomfilter size
 * @param nf candidate for number of successive hash functions
 * @param r candidate for number of random requests
 * @return A boolean if parameters are acceptable
 */
bool is_correct_parameters(string filename, uint64_t k, uint64_t n, uint64_t nf, uint64_t r)
{
    ifstream ifile;
    ifile.open(filename);
    if (!ifile)
    {
        std::cout << "Specified file doesn't exist ; exiting...";
        return 0;
    }
    else
    {
        std::cout << "filename : " << filename << "\n";
    }

    if (k > 31 || k < 0)
    {
        std::cout << "Invalid value for parameter k ; exiting...";
        return 0;
    }
    else
    {
        std::cout << "k : " << k << "\n";
    }

    if (n > pow(2, 34) || n < 0)
    {
        std::cout << "Invalid value for parameter n ; exiting...";
        return 0;
    }
    else
    {
        std::cout << "n : " << n << "\n";
    }

    if (nf > 64 || nf < 0)
    {
        std::cout << "Invalid value for parameter nf ; exiting...";
        return 0;
    }
    else
    {
        std::cout << "nf : " << nf << "\n";
    }

    if (r < 0)
    {
        std::cout << "Invalid value for parameter r ; exiting...";
        return 0;
    }
    else
    {
        std::cout << "r : " << r << "\n";
    }
    return 1;
}

/** Constructs the bloomfilter for a sequence, and then do a set of random requests
 * @param argc Number of args in vector
 * @param argv A vector containing command line arguments
 */
int main(int argc, char *argv[])
{
    // casting console parameters
    string filename = argv[1];
    uint64_t k = stoi(argv[2]);
    uint64_t n = stoi(argv[3]);
    uint64_t nf = stoi(argv[4]);
    uint64_t r = stoi(argv[5]);

    // verifiying parameters
    if (!is_correct_parameters(filename, k, n, nf, r))
    {
        return 0;
    }

    // creating the bloomfilter object
    BloomFilter my_bf(n, nf);

    long last_char = get_file_length(filename);

    string kmer;
    uint64_t previous_position = 0;
    bool is_correct_kmer = true;
    while (is_correct_kmer)
    {
        tie(kmer, previous_position) = next_kmer(filename, k, previous_position);
        if (kmer != "")
        {
            my_bf.add_value(encode_kmer(kmer));
        }
        else
        {
            // we reached end of file
            is_correct_kmer = false;
        }
    }
    uint64_t good_tries = 0;
    // Executing research requests on the bloom filter
    for (int i = 0; i < r; i++)
    {
        if (my_bf.is_present(encode_kmer(generate_kmer(k))))
        {
            good_tries += 1;
        }
    }
    cout << ">>> Across " << r << " randomly generated k-mers, " << good_tries << " were in filter!";
    return 0;
}
