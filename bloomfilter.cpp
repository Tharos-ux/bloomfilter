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
        // will store the size of the bloomfilter
        bloom_size = n;
        // filter itself
        filter = vector<bool>(n, false);
        // number of successive hash functions
        hash_func_number = nf;
    }

    /** Adds a value to the filter
     * @param value_to_hash encoded kmer
     */
    void add_value(uint64_t value_to_hash)
    {
        // std::cout << "Define?";
        uint64_t *hashes = new uint64_t[bloom_size / 8 + 1];
        multihash(value_to_hash, hashes, hash_func_number, bloom_size - 1);
        // std::cout << "Crash?";
        for (uint64_t i_hash = 0; i_hash < bloom_size / 8 + 1; i_hash++)
        {
            cout << "i_hash = " << i_hash << "\n";
            // getting binary hashed value at index i
            bitset<8> binary_code = bitset<8>(hashes[i_hash]);
            // changing filter positions, iterating through the vector
            for (uint64_t j_boolvector = i_hash * 8 + 1; j_boolvector < i_hash * 8 + 9; j_boolvector++)
            {
                cout << "   j_boolvector = " << j_boolvector << " ; bitset pos = " << j_boolvector - (i_hash * 8) - 1 << " ; bitset content = " << binary_code[j_boolvector - (i_hash * 8) - 1] << "\n";
                if (binary_code[j_boolvector - (i_hash * 8) - 1] == 1)
                {
                    filter[j_boolvector] = true;
                }
                cout << "   >> Returning component " << filter[j_boolvector] << "\n";
            }
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
        for (uint64_t i_hash = 0; i_hash < bloom_size / 8 + 1; i_hash++)
        {
            // getting binary hashed value at index i
            bitset<8> binary_code = bitset<8>(hashes[i_hash]);

            for (uint64_t j_boolvector = i_hash * 8 + 1; j_boolvector < i_hash * 8 + 9; j_boolvector++)
            {
                if (binary_code[j_boolvector - i_hash * 8])
                {
                    if (!filter[j_boolvector])
                    {
                        // getting one value not equal to 1 here means kmer is absent
                        return false;
                    }
                }
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
        std::cout << "Error while reading file.\n";
    }
    return std::make_tuple('\0', fasta_file.tellg());
}

// TODO rework func in order to re-use precalculated kmer if exists
/** Returns the next clean kmer from the file
 * @param filename Path to the file to read
 * @param k k-mer length
 * @param previous_position Starting point of the previous kmer
 * @return A tuple containing next kmer and cursor position of the start of the kmer
 */
std::tuple<string, long> next_kmer(string filename, int k, long previous_position)
{
    // std::cout << "Entering?";
    long cursor_previous_kmer = previous_position;
    long cursor;
    char current_char;
    string kmer;

    // std::cout << "Reading?";
    // getting starting point for next kmer
    tie(current_char, cursor) = read_fasta(filename, cursor_previous_kmer);
    cursor_previous_kmer = cursor;
    kmer += current_char;
    for (int i = 1; i < k; i++)
    {
        tie(current_char, cursor) = read_fasta(filename, cursor);
        kmer += current_char;
    }

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

/** Gets arguments out of command line.
 * @param argc Number of args in vector
 * @param argv A vector containing command line arguments
 */
int main(int argc, char *argv[])
{
    string filename = argv[1];
    int k = stoi(argv[2]);
    int n = stoi(argv[3]);
    int nf = stoi(argv[4]);
    int r = stoi(argv[5]);

    // verifications on parameters
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

    BloomFilter my_bf(n, nf);

    long last_char = get_file_length(filename);

    string kmer;
    long previous_position = 0;
    while (previous_position < last_char - k + 1)
    {
        // std::cout << "Iterate?";
        tie(kmer, previous_position) = next_kmer(filename, k, previous_position);
        // std::cout << "Add?";
        my_bf.add_value(encode_kmer(kmer));
        std::cout << "\nKmer is " << kmer << " at position " << previous_position << " and reverse complement is " << reverse_complement(kmer) << " we select " << choose_kmer(kmer) << " encoded as " << encode_kmer(kmer) << "\n";
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
    std::cout << "Across " << r << " randomly generated k-mers, " << good_tries << " were in filter!";
    return 0;
}
