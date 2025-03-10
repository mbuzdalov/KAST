#include "utils.h"
#include "distance.h"

/*
I need to check that if we are using skip-mers, then we need to check that these are sensible.
  * skip-mer mask must all be 0/1's
  * skip-mer mask should be the same number of characters as the kmer size
  * all skipmers shoud have the same number of 1's across masks
*/
bool parse_mask(modify_string_options const &options, int &effective_klen) {
   bool first = true;

   for (auto m : options.mask) {
      // all should be the same number of characters as the klen
      if (m.length() != options.klen) {
         seqan3::debug_stream << "ERROR: Mask sizes should be the same size "
                                 "as the K-Mer length." << std::endl;
         return false;
      }

      int counter = 0;

      // checks to see that the mask is only made of 0/1's
      for (unsigned i = 0; i < m.length(); ++i) {
         if (m[i] != '0' && m[i] != '1') {
            seqan3::debug_stream << "ERROR: Masks should only contain 0's or 1's." << std::endl;
            return false;
         }

         if (m[i] == '1') {
            counter++;
         }
      }

      if (first) {
         effective_klen = counter;
         first = false;
      } else {
         if (counter != effective_klen) {
            seqan3::debug_stream << "ERROR: The number of 0's and 1's in each mask "
                                    "should be the same e.g. 10001, 11000, 00011" << std::endl;
            return false;
         }
      }
   }
   return true;
};

/*
 * Splits the given string into non-whitespace tokens. Any character less than or equal to ' ' is considered whitespace.
 */
std::vector<std::string> split_by_space(std::string const &string) {
    std::vector<std::string> result;
    size_t start = 0;
    while (true) {
        while (start < string.length() && string[start] <= ' ') {
            ++start;
        }
        if (start == string.length()) {
            return result;
        }
        size_t end = start;
        while (end < string.length() && string[end] > ' ') {
            ++end;
        }
        result.push_back(string.substr(start, end - start));
        start = end;
    }
}

/*
 * Prepares a sequence for further processing given the k-mer length.
 * For seqan3::dna5, turns it into a reverse complement and inserts k 'N' symbols between.
 */
template <>
void prepare_sequence_inplace(std::vector<seqan3::dna5> &sequence, unsigned k) {
    size_t initial_length = sequence.size();
    sequence.reserve(2 * initial_length + k);
    auto N = seqan3::dna5{}.assign_char('N');
    for (unsigned i = 0; i < k; ++i) {
        sequence.push_back(N);
    }
    while (initial_length) {
        sequence.push_back(sequence[--initial_length].complement());
    }
}

// This is a specialization that treats N of dna5 as a non-value
// and computes hashes assuming base 4.
// This is not possible with seqan3 out of the box without writing own code,
// and the computation is reasonably trivial.
template <>
void count_kmers(std::vector<unsigned> & kmer_counts, std::vector<seqan3::dna5> const &sequence,
                 unsigned const k) {
   uint64_t kmer_number = uint64_t(1) << (2 * k);
   kmer_counts.clear();
   kmer_counts.resize(kmer_number, 0);

   uint64_t curr_hash = 0;
   uint64_t hash_mask = kmer_number - 1;
   uint32_t banned_mask = (uint32_t(1) << k) - 1;
   uint32_t banned_bits = 0;

   for (size_t i = 0; i + 1 < k; ++i) {
      curr_hash <<= 2;
      banned_bits = (banned_bits << 1) & banned_mask;
      char base = sequence[i].to_char();
      if (base == 'N') {
         banned_bits |= 1;
      } else {
         // this makes all indexing compatible to dna4
         curr_hash += rank_of_char_as<seqan3::dna4>(base);
      }
   }

   for (size_t i = k - 1; i < sequence.size(); ++i) {
      curr_hash = (curr_hash << 2) & hash_mask;
      banned_bits = (banned_bits << 1) & banned_mask;
      char base = sequence[i].to_char();
      if (base == 'N') {
         banned_bits |= 1;
      } else {
         // this makes all indexing compatible to dna4
         curr_hash += rank_of_char_as<seqan3::dna4>(base);
      }
      if (!banned_bits) safe_increment(kmer_counts, curr_hash);
   }
};

// For Dna
template <>
void count_kmers(std::vector<unsigned> &kmer_counts,
                 std::vector<seqan3::dna5> const& sequence,
                 unsigned const k,
                 unsigned const effective_k,
                 std::vector<std::string> const &masks) {
   uint64_t kmer_number = uint64_t(1) << (2 * effective_k);
   kmer_counts.clear();
   kmer_counts.resize(kmer_number, 0);

   uint64_t curr_hash = 0;
   uint64_t hash_mask = (uint64_t(1) << (2 * k)) - 1;
   uint32_t banned_mask = (uint32_t(1) << k) - 1;
   uint32_t banned_bits = 0;

   std::vector<uint32_t> parsed_masks;
   for (std::string const &mask : masks) {
      uint32_t bin = 0;
      for (size_t i = 0; i < mask.length(); ++i) {
         // we have these masks reversed here,
         // because they are not used through seqan3 mask API,
         // but in a custom way which wants the opposite order.
         if (mask[i] == '1') bin |= 1 << (mask.length() - 1 - i);
      }
      parsed_masks.push_back(bin);
   }

   for (size_t i = 0; i + 1 < k; ++i) {
      curr_hash <<= 2;
      banned_bits = (banned_bits << 1) & banned_mask;
      char base = sequence[i].to_char();
      if (base == 'N') {
         banned_bits |= 1;
      } else {
         // this makes all indexing compatible to dna4
         curr_hash += rank_of_char_as<seqan3::dna4>(base);
      }
   }
   // Note that here we only count those masked k-mers where:
   //  - the entire k bits, masked or not, are contained in the sequence (e.g. if the mask is 11110100,
   //       we do not count in the k-mer resulting in last 111101 masked bits of the sequence)
   //  - the entire k bits should be N-free, even if all Ns are masked away.
   //
   // Whether this is sound from the bioinf perspective, I don't know,
   // so backward compatibility is preserved just in case.
   for (size_t i = k - 1; i < sequence.size(); ++i) {
      curr_hash = (curr_hash << 2) & hash_mask;
      banned_bits = (banned_bits << 1) & banned_mask;
      char base = sequence[i].to_char();
      if (base == 'N') {
         banned_bits |= 1;
      } else {
         // this makes all indexing compatible to dna4
         curr_hash += rank_of_char_as<seqan3::dna4>(base);
      }
      if (!banned_bits) {
         for (uint32_t mask : parsed_masks) {
            uint64_t copy_hash = curr_hash;
            uint64_t result = 0;
            while (mask) {
               if (mask & 1) {
                  result <<= 2;
                  result |= copy_hash & 3;
               }
               mask >>= 1;
               copy_hash >>= 2;
            }
            safe_increment(kmer_counts, result);
         }
      }
   }
}

// for DNA sequences
template <>
void markov<>(std::vector<double> &markov_counts,
              std::vector<seqan3::dna5> const &sequence,
              unsigned const k,
              unsigned const markov_order) {
   // init the result vector
   uint64_t kmer_number = uint64_t(1) << (2 * k);
   markov_counts.clear();
   markov_counts.resize(kmer_number, 0);

   // create the background model
   std::vector<unsigned> markov_bg;
   count_kmers(markov_bg, sequence, markov_order + 1);
   unsigned total_count = std::reduce(markov_bg.begin(), markov_bg.end());

   // we simulate the unhash and do it for dna4
   std::vector<unsigned> markov_query;
   std::vector<seqan3::dna4> kmer(k);

   for (uint64_t kmer_idx = 0; kmer_idx < kmer_number; ++kmer_idx) {
      // this calls the dna4 version (the generic one), which will do it just right
      count_kmers(markov_query, kmer, markov_order + 1);
      double prob = 1;
      for (size_t i = 0; i < markov_query.size(); ++i) {
         prob *= std::pow((double) (markov_bg[i]) / total_count, markov_query[i]); //seqan3::pow?
      }
      markov_counts[kmer_idx] = prob;

      // this is the "+1" operation on the kmer
      for (size_t i = 0; i < k; ++i) {
         unsigned rank = seqan3::to_rank(kmer[i]);
         if (rank < 3) {
            kmer[i].assign_rank(rank + 1);
            break;
         } else {
            kmer[i].assign_rank(0);
         }
      }
   }
}


double distance_dispatch(std::string const &type,
                         std::vector<unsigned> const &counts_l, std::vector<unsigned> const &counts_r,
                         std::vector<double> const &markov_l, std::vector<double> const &markov_r) {
    if (type == "euclid")           return euler(counts_l, counts_r);
    else if (type == "d2")          return d2(counts_l, counts_r);
    else if (type == "cosine")      return cosine(counts_l, counts_r);
    else if (type == "manhattan")   return manhattan(counts_l, counts_r);
    else if (type == "chebyshev")   return chebyshev(counts_l, counts_r);
    else if (type == "canberra")    return canberra(counts_l, counts_r);
    else if (type == "normalised_canberra")     return normalised_canberra(counts_l, counts_r);
    else if (type == "bc")                      return bray_curtis_distance(counts_l, counts_r);
    else if (type == "ngd")                     return normalised_google_distance(counts_l, counts_r);
    else if (type == "d2s" || type == "D2S")    return d2s(counts_l, counts_r, markov_l, markov_r);
    else if (type == "hao")                         return hao(counts_l, counts_r, markov_l, markov_r);
    else if (type == "d2star" || type == "D2Star")  return d2star(counts_l, counts_r, markov_l, markov_r);
    else if (type == "dai")                         return dai(counts_l, counts_r, markov_l, markov_r);

    seqan3::debug_stream << "Error: distance_dispatch does not know distance type '" << type << "'" << std::endl;
    exit(1);
}
