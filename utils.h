#ifndef __KAST_UTILS_H__
#define __KAST_UTILS_H__

#include "common.h"

#include <array>
#include <numeric>
#include <vector>
#include <unordered_map>

#include <sys/sysinfo.h> // eventually for memory checking

#include <seqan3/alphabet/all.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/utility/math.hpp>
#include <seqan3/utility/views/all.hpp>

template <seqan3::writable_alphabet T>
struct input_traits {};

template<>
struct input_traits<seqan3::dna5> {
    using traits_type = seqan3::sequence_file_input_default_traits_dna;
};

template<>
struct input_traits<seqan3::aa27> {
    using traits_type = seqan3::sequence_file_input_default_traits_aa;
};

struct sequence_file_input_default_traits_aa10 : seqan3::sequence_file_input_default_traits_aa {
    using sequence_alphabet = seqan3::aa10murphy;
    using sequence_legal_alphabet = seqan3::aa27;
};

template<>
struct input_traits<seqan3::aa10murphy> {
    using traits_type = sequence_file_input_default_traits_aa10;
};

template<seqan3::alphabet T>
using sequence_file_type = seqan3::sequence_file_input<typename input_traits<T>::traits_type>;

/*
 * Returns the index of a char, assuming it is from the specified alphabet.
 */
template <seqan3::writable_alphabet T>
unsigned rank_of_char_as(char ch) {
    return seqan3::assign_char_to(ch, T{}).to_rank();
}

/*
 * Prepares a sequence for further processing given the k-mer length.
 * Does nothing by default.
 */
template <seqan3::writable_alphabet T>
void prepare_sequence_inplace(std::vector<T> &sequence, unsigned k) {}

/*
 * Prepares a sequence for further processing given the k-mer length.
 * For seqan3::dna5, turns it into a reverse complement and inserts k 'N' symbols between.
 */
template <>
void prepare_sequence_inplace(std::vector<seqan3::dna5> &sequence, unsigned k);

/*
I need to check that if we are using skip-mers, then we need to check that these are sensible.
  * skip-mer mask must all be 0/1's
  * skip-mer mask should be the same number of characters as the kmer size
  * all skipmers shoud have the same number of 1's across masks
*/
bool parse_mask(modify_string_options const &options, int &effective_klen);

/*
 * Splits the given string into non-whitespace tokens. Any character less than or equal to ' ' is considered whitespace.
 */
std::vector<std::string> split_by_space(std::string const &string);

template <class T>
void safe_increment_impl(std::vector<T> &destination, size_t index, char const *file, unsigned line) {
    if (index >= destination.size()) {
        seqan3::debug_stream << "File " << file << " line " << line << ": Error: Index too large (" << index << " out of " << destination.size() << ")" << std::endl;
        exit(1);
    } else if (destination[index] < std::numeric_limits<T>::max()) {
        ++destination[index];
    } else {
        seqan3::debug_stream << "File " << file << " line " << line << ": Error: Integer overflow detected. Exiting" << std::endl;
        exit(1);
    }
}

#define safe_increment(dest, index) safe_increment_impl(dest, index, __FILE__, __LINE__)

template <seqan3::nucleotide_alphabet T>
double gc_ratio(std::vector<T> const &sequence) {
   // this way it should be faster, if speed is ever an issue here
   std::array<size_t, T::alphabet_size> count;
   for (auto symbol : sequence) {
      ++count[symbol.to_rank()];
   }
   double gc = count[rank_of_char_as<T>('G')]
             + count[rank_of_char_as<T>('C')];
   double at = count[rank_of_char_as<T>('A')]
             + count[rank_of_char_as<T>('T')];
   return gc / (gc + at);
}

template <seqan3::alphabet T>
double gc_ratio(std::vector<T> const &sequence) {
    return 0; // makes no sense for non-nucleotide alphabets, but is required by search printers
}

/*
Perform regular counting
*/
template <seqan3::alphabet T>
void count_kmers(std::vector<unsigned> &kmer_counts, std::vector<T> const &sequence,
                 unsigned const k) {
   seqan3::ungapped ungapped_k{ uint8_t(k) };
   seqan3::shape my_shape{ungapped_k};

   uint64_t kmer_number = seqan3::pow(T::alphabet_size, k);

   kmer_counts.clear();
   kmer_counts.resize(kmer_number, 0);

   auto hash_view = sequence | seqan3::views::kmer_hash(my_shape);
   for (auto hash : hash_view) {
      safe_increment(kmer_counts, hash);
   }
};

// This is a specialization that treats N of dna5 as a non-value
// and computes hashes assuming base 4.
// This is not possible with seqan3 out of the box without writing own code,
// and the computation is reasonably trivial.
template <>
void count_kmers(std::vector<unsigned> & kmer_counts, std::vector<seqan3::dna5> const &sequence, unsigned const k);

/*
Perform counting with a mask array
For AminoAcid and ReducedAminoAcid
*/
template <seqan3::alphabet T>
void count_kmers(std::vector<unsigned> &kmer_counts,
                 std::vector<T> const &sequence,
                 unsigned const k,
                 unsigned const effective_k,
                 std::vector<std::string> const &masks) {
   uint64_t kmer_number = seqan3::pow(T::alphabet_size, effective_k);
   kmer_counts.clear();
   kmer_counts.resize(kmer_number, 0);

   for (std::string const &mask : masks) {
      seqan3::bin_literal bin;
      bin.value = 0;
      for (size_t i = 0; i < mask.length(); ++i) {
         if (mask[i] == '1') bin.value |= size_t(1) << i;
      }
      // Here the thing is.
      //
      // KAST assumes that masks are masking out k-mers, that is, if the mask ends in zeros,
      // we do not have a permission to move the mask so that these zeros fall off the sequence.
      //
      // seqan3, on the other hand, does not allow expressing this in its shapes.
      // So we have to cut the sequence from the end if our mask has zeros at the end.

      // TODO: the current implementation will likely fail if a mask starts with a zero,
      // because shape creation will fail.

      size_t new_size = sequence.size() - k + std::bit_width(bin.value);
      auto hash_view = sequence | std::views::take(new_size) | seqan3::views::kmer_hash(seqan3::shape(bin));
      for (auto hash : hash_view) {
         safe_increment(kmer_counts, hash);
      }
   }
};

// For Dna
template <>
void count_kmers(std::vector<unsigned> &kmer_counts,
                 std::vector<seqan3::dna5> const& sequence,
                 unsigned const k,
                 unsigned const effective_k,
                 std::vector<std::string> const &masks);

// for all others
template <seqan3::writable_alphabet T>
void markov(std::vector<double> &markov_counts,
            std::vector<T> const &sequence,
            unsigned const k,
            unsigned const markov_order) {
   // init the result vector
   uint64_t kmer_number = seqan3::pow(T::alphabet_size, k);
   markov_counts.clear();
   markov_counts.resize(kmer_number, 0);

   // create the background model
   std::vector<unsigned> markov_bg;
   count_kmers(markov_bg, sequence, markov_order + 1);
   unsigned total_count = std::reduce(markov_bg.begin(), markov_bg.end());

   // seqan3 does not have anything for unhash, so we simulate it
   std::vector<unsigned> markov_query;
   std::vector<T> kmer(k);

   for (uint64_t kmer_idx = 0; kmer_idx < kmer_number; ++kmer_idx) {
      count_kmers(markov_query, kmer, markov_order + 1);
      double prob = 1;
      for (size_t i = 0; i < markov_query.size(); ++i) {
         prob *= std::pow((double) (markov_bg[i]) / total_count, markov_query[i]); //seqan3::pow?
      }
      markov_counts[kmer_idx] = prob;

      // this is the "+1" operation on the kmer
      for (size_t i = 0; i < k; ++i) {
         unsigned rank = seqan3::to_rank(kmer[i]);
         if (rank + 1 < T::alphabet_size) {
            kmer[i].assign_rank(rank + 1);
            break;
         } else {
            kmer[i].assign_rank(0);
         }
      }
   }
}

// for DNA sequences
template <>
void markov<>(std::vector<double> &markov_counts,
              std::vector<seqan3::dna5> const &sequence,
              unsigned const k,
              unsigned const markov_order);

template <seqan3::writable_alphabet T>
void populate_counts(modify_string_options const &options, std::vector<T> const &sequence,
                     std::vector<unsigned> &counts, std::vector<double> &markov_counts) {
    if (options.mask.size()) {
        count_kmers(counts, sequence, options.klen, options.effective_length, options.mask);
    } else {
        count_kmers(counts, sequence, options.klen);
    }

    if (options.type == "d2s" || options.type == "D2S" ||
        options.type == "d2star" || options.type == "D2Star" ||
        options.type == "hao" || options.type == "dai")
    {
        if (options.mask.size()) {
            markov(markov_counts, sequence, options.effective_length, options.markov_order);
        } else {
            markov(markov_counts, sequence, options.klen, options.markov_order);
        }
    }
}

double distance_dispatch(std::string const &type,
                         std::vector<unsigned> const &counts_l, std::vector<unsigned> const &counts_r,
                         std::vector<double> const &markov_l, std::vector<double> const &markov_r);


/// Experimental idea

/*
Ignore the order of AA in the kmer for proteins only, so that e.g. TY and YT are the same, 
and for 3mers CMY would be the equivalent to all 3x2x1 combinations 
e.g. CMY, CYM, MCY, MYC, YCM, YMC. This could be useful for 4 mers and 5 mers because it 
starts to reduce the number of kmers by quite a bit e.g. for 4-mers you have 4x3x2 = 24 
redundant combinations, for 5mers its 120. This might make it tractable to use 5-mers. 
E.g with this redundancy we would count, for 3mers 8000/6 = 1333 kmers, for 4mers that 
becomes 6,666 and for 5 mers, 26,666 
(as compared to 3,200,000 kmers without the redundancy for 5 mers).  

I think the best way to do this is to put the kmers in alphabetical order. So in the 
example above, CMY could be the correct order, and CYM and the others will all be reordered to CMY.

*/

/*
template <typename TAlphabet>
int countReducedAlphabet(String<unsigned> & kmerCounts, String<TAlphabet> const & sequence, unsigned const k)
{
   Shape<TAlphabet> myShape;
   resize(myShape, k);
   unsigned long long int kmerNumber = _intPow((unsigned)ValueSize<TAlphabet>::VALUE, weight(myShape));

   seqan2::clear(kmerCounts);
   seqan2::resize(kmerCounts, kmerNumber, 0);

   auto itSequence = begin(sequence);

   for (; itSequence <= (end(sequence) - k); ++itSequence)
   {
      std::cout << itSequence << "\t" << value(itSequence) << "\t";
      DnaString meh;
      seqan2::unhash(meh, seqan2::hash(myShape, itSequence), weight(myShape));
      std::cout << meh << std::endl;
      //cout << seqan::unhash(meh, seqan::hash(myShape, itSequence), itSequence) << endl;
      long long unsigned int hashValue = seqan2::hash(myShape, itSequence);

      safe_increment(kmerCounts, hashValue);
   }
   return 0;
};

*/

#endif // __KAST_UTILS_H__
