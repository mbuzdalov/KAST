#ifndef __KAST_UTILS_H__
#define __KAST_UTILS_H__

#include "common.h"

#include <array>
#include <numeric>
#include <vector>
#include <unordered_map>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/utility/math.hpp>
#include <seqan3/utility/views/all.hpp>

// Parse our commandline options
bool parse_command_line(modify_string_options &options, int argc, char const ** argv) {
   seqan3::argument_parser parser{"kast", argc, argv};
   parser.info.version = "1.0.3beta";
   parser.info.date = "Mar 2025";
   //parser.short_description = "Kmer Alignment-free Search Tool";
   parser.add_line("Perform Alignment-free k-tuple frequency comparisons from sequences.");
   parser.add_line("This can be in the form of two input files (e.g. a reference and a query):");
   parser.add_line("   -q query.fasta -r reference.fasta -o results.txt [\\fIOPTIONS\\fP] ");
   parser.add_line("or a single file for pairwise comparisons to be made:");
   parser.add_line("   -p mydata.fasta -o results.txt [\\fIOPTIONS\\fP] ");

   //addOption(parser, ArgParseOption("k", "klen", "Kmer Length.",
   //                                 ArgParseArgument::INTEGER, "INT"));
   //setDefaultValue(parser, "klen", "3");
   //getOptionValue(options.klen, parser, "klen");
   options.klen = 3; //default
   parser.add_option(options.klen, 'k', "klen", "Kmer Length");

   //addOption(parser, ArgParseOption("d", "debug", "Debug Messages."));
   //options.debug = isSet(parser, "debug");
   options.debug = false;
   parser.add_flag(options.debug, 'd', "debug", "Debug Messages");

   //addOption(parser, ArgParseOption("q", "query-file",
   //                                 "Path to the file containing your query \
   //                                 sequence data.\n",
   //                                 ArgParseArgument::INPUT_FILE, "IN"));
   //getOptionValue(options.queryFileName, parser, "query-file");
   options.query_filename = "";
   parser.add_option(options.query_filename, 'q', "query-file",
        "Path to the file containing your query sequence data");

   //addOption(parser, ArgParseOption("r", "reference-file",
   //                                 "Path to the file containing your reference\
   //                                  sequence data.",
   //                                  ArgParseArgument::INPUT_FILE, "IN"));
   //getOptionValue(options.referenceFileName, parser, "reference-file");
   options.reference_filename = "";
   parser.add_option(options.reference_filename, 'r', "reference-file",
        "Path to the file containing your reference sequence data");

   //addOption(parser, ArgParseOption("p", "pairwise-file",
   //                                 "Path to the file containing your sequence \
   //                                 data which you will perform pairwise \
   //                                 comparison on.",
   //                                 ArgParseArgument::INPUT_FILE, "IN"));
   //getOptionValue(options.pairwiseFileName, parser, "pairwise-file");
   options.pairwise_filename = "";
   parser.add_option(options.pairwise_filename, 'p', "pairwise-file",
        "Path to the file containing your sequence data which you will perform pairwise comparison on");

   //addOption(parser, ArgParseOption("i", "interleaved-file",
   //                                 "Path to the file containing your sequence \
   //                                 data which is interleaved.",
   //                                 ArgParseArgument::INPUT_FILE, "IN"));
   //getOptionValue(options.interleavedFileName, parser, "interleaved-file");
   options.interleaved_filename = "";
   parser.add_option(options.interleaved_filename, 'i', "interleaved-file",
        "Path to the file containing your sequence data which is interleaved");

   //setDefaultValue(parser, "markov-order", "0");
   //addOption(parser, ArgParseOption("m", "markov-order", "Markov Order",
   //          ArgParseArgument::INTEGER, "INT"));
   //getOptionValue(options.markovOrder, parser, "markov-order");
   options.markov_order = 0;
   parser.add_option(options.markov_order, 'm', "markov-order", "Markov Order");

   //addOption(parser, ArgParseOption("o", "output-file", "Output file.",
   //          ArgParseArgument::OUTPUT_FILE, "OUT"));
   //getOptionValue(options.outputFileName, parser, "output-file");
   options.output_filename = "";
   parser.add_option(options.output_filename, 'o', "output-file", "Output file");

   //addOption(parser, ArgParseOption("n", "num-hits",
   //          "Number of top hits to return when running a Ref/Query search.\
   //           If you want all the result, enter 0.", ArgParseArgument::INTEGER, "INT"));
   //setDefaultValue(parser, "num-hits", "10");
   //getOptionValue(options.nohits, parser, "num-hits");
   options.nohits = 10;
   parser.add_option(options.nohits, 'n', "num-hits",
        "Number of top hits to return when running a Ref/Query search. If you want all the results, enter 0");

   //addOption(parser, ArgParseOption("t", "distance-type",
   //                                 "The method of calculating the distance \
   //                                 between two sequences. For descriptions of \
   //                                 distance please refer to the wiki. ",
   //                                 ArgParseArgument::STRING, "STR"));
   //setValidValues(parser, "distance-type",
   //setDefaultValue(parser, "distance-type", "d2");
   //getOptionValue(options.type, parser, "distance-type");
   seqan3::value_list_validator<std::string> type_values {
        "d2", "euclid", "d2s", "d2star", "manhattan", "chebyshev",
        "dai", "bc", "ngd", "all", "canberra", "normalised_canberra", "cosine"
   };
   options.type = "d2";
   parser.add_option(options.type, 't', "distance-type",
        "The method of calculating the distance between two sequences. "
        "For descriptions of distance please refer to the wiki", seqan3::option_spec::standard, type_values);

   //addOption(parser, ArgParseOption("sc", "score-cutoff", "Score Cutoff for search mode.",
   //          ArgParseArgument::DOUBLE, "DOUBLE"));
   //if(isSet(parser, "score-cutoff") == false)
   //{
   //   options.score_cutoff = std::numeric_limits<double>::quiet_NaN();
   //}
   //getOptionValue(options.score_cutoff, parser, "score-cutoff");
   options.score_cutoff = std::numeric_limits<double>::quiet_NaN();
   parser.add_option(options.score_cutoff, 'S', "score-cutoff", "Score Cutoff for search mode");

   //addOption(parser, ArgParseOption("s", "sequence-type",
   //          "Define the type of sequence data to work with.",
   //          ArgParseArgument::STRING, "STR"));
   //setValidValues(parser, "sequence-type", "dna aa raa");
   //setDefaultValue(parser, "sequence-type", "dna");
   //getOptionValue(options.sequenceType, parser, "sequence-type");
   seqan3::value_list_validator<std::string> sequence_type_values { "dna", "aa", "raa" };
   options.sequence_type = "dna";
   parser.add_option(options.sequence_type, 's', "sequence-type",
        "The type of sequence data to work with", seqan3::option_spec::standard, sequence_type_values);


   //addOption(parser, ArgParseOption("f", "output-format",
   //          "For Reference/query based usage you can select your output type.",
   //          ArgParseArgument::STRING, "STR"));
   //getOptionValue(options.output_format, parser, "output-format");
   //setValidValues(parser, "output-format", "default tabular blastlike");
   //setDefaultValue(parser, "output-format", "default");
   seqan3::value_list_validator<std::string> output_format_values { "default", "tabular", "blastlike" };
   options.output_format = "default";
   parser.add_option(options.output_format, 'f', "output-format",
        "The output type for Reference/query based usage", seqan3::option_spec::standard, output_format_values);

   //addOption(parser, ArgParseOption("nr", "no-reverse",
   //          "Do not use reverse compliment."));
   //options.noreverse = isSet(parser, "no-reverse");
   options.noreverse = false;
   parser.add_flag(options.noreverse, 'N', "no-reverse",
        "Do not use reverse compliment");

   //addOption(parser, ArgParseOption("gc", "calc-gc",
   //          "Calculate GC content of query/ref in search mode."));
   //options.calcgc = isSet(parser, "calc-gc");
   options.calcgc = false;
   parser.add_flag(options.calcgc, 'G', "calc-gc",
        "Calculate GC content of query/ref in search mode");

   //addOption(parser, ArgParseOption("nh", "no-header",
   //          "Do not print header when performing search mode."));
   //options.noheader = isSet(parser, "no-header");
   options.noheader = false;
   parser.add_flag(options.noheader, 'H', "no-header",
        "Do not print header when performing search mode");

   //addOption(parser, ArgParseOption("mask", "skip-mer",
   //          "Specify binary masks where a zero indicates \
   //          skipping that base and one keeps it. e.g. 01110.",
   //          ArgParseArgument::STRING, "TEXT", true));
   //for(int i = 0; i < getOptionValueCount(parser, "skip-mer"); i++)
   //{
   //   CharString tmpVal;
   //   getOptionValue(tmpVal, parser, "skip-mer", i);
   //   options.mask.push_back(tmpVal);
   //}
   parser.add_option(options.mask, 'M', "skip-mer",
        "Specify binary masks where a zero indicates skipping that base and one keeps it, e.g. 01110");

   //addOption(parser, ArgParseOption("c", "num-cores", "Number of Cores.",
   //          ArgParseArgument::INTEGER, "INT"));
   //setDefaultValue(parser, "num-cores", "1");
   //getOptionValue(options.num_threads, parser, "num-cores");
   options.num_threads = 1;
   parser.add_option(options.num_threads, 'c', "num-cores", "Number of cores");

   //addOption(parser, ArgParseOption("fp", "filter-percent", "In search mode, only match\
   //                                 those results where the query and ref sequence \
   //                                 lengths are within +/- percentage of oneanother.",
   //          ArgParseArgument::DOUBLE, "DOUBLE"));
   //getOptionValue(options.filter_percent, parser, "filter-percent");
   parser.add_option(options.filter_percent, '%', "filter-percent",
        "In search mode, only match those results where the query and ref sequence lengths are within +/- percentage of one another");

   //addOption(parser, ArgParseOption("fb", "filter-bp", "In search mode, only match\
   //                                 those results where the query and ref sequence \
   //                                 lengths are within +/- bp of oneanother.",
   //          ArgParseArgument::INT64, "INT64"));
   //getOptionValue(options.filter_bp, parser, "filter-bp");
   parser.add_option(options.filter_bp, 'B', "filter-bp",
        "In search mode, only match those results where the query and ref sequence lengths are within +/- bp of one another");

   parser.parse();

   /*
   Check to see if markov order is correct.
   */
   if (options.type == "d2s" || options.type == "d2star" ||
       options.type == "dai" || options.type == "hao" ||
       options.type == "all")
   {
      if (options.markov_order < 0) {
         seqan3::debug_stream << "Markov order must be >= 0 " << std::endl;
         return false;
      }

      if (options.mask.size() > 0) {
         if (options.markov_order >= options.effective_length - 1) {
            seqan3::debug_stream << "Markov order must be < effectiveLength-1 " << std::endl;
            return false;
         }
      } else {
         if (options.markov_order >= options.klen - 1) {
            seqan3::debug_stream << "Markov order must be < klen-1 " << std::endl;
            return false;
         }
      }
   }

   if (parser.is_option_set("pairwise-file")) {
      if (parser.is_option_set("reference-file") ||
          parser.is_option_set("query-file"))
      {
         seqan3::debug_stream << "If you are performing a pairwise comparison, "
                                 "you do not need to specify a query (-q) and a "
                                 "reference (-r) file. If you are performing a "
                                 "reference/query based search you do not need to "
                                 "specify a pairwise-file (-p). See kast -h for "
                                 "details." << std::endl;
         return false;
      }
   }

   if (parser.is_option_set("reference-file") == true && !parser.is_option_set("query-file")) {
      seqan3::debug_stream << "You have specified a reference (-r) file but "
                              "not a query (-q) file. See kast -h for details." << std::endl;
      //printHelp(parser);
      return false;
   }

   if (!parser.is_option_set("reference-file") && parser.is_option_set("query-file")) {
      seqan3::debug_stream << "You have specified a query (-q) file but "
                              "not a reference (-r) file. See kast -h for details." << std::endl;
      //printHelp(parser);
      return false;
   }

   if (!parser.is_option_set("reference-file") &&
       !parser.is_option_set("query-file") &&
       !parser.is_option_set("pairwise-file") &&
       !parser.is_option_set("interleaved-file"))
   {
      seqan3::debug_stream << "You have not specifed any input file. "
                              "See kast -h for details." << std::endl;
      //printHelp(parser);
      return false;
   }

   return true;
}

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
