/*
 * KAST - Kmer Alignment-free Search Tool
 * Version 1.0.3beta
 * Original version by Dr. Martin Vickers (martin.vickers@jic.ac.uk)
 * Current version rewritten for seqan3 by Dr. Maxim Buzdalov (mab168@aber.ac.uk)
 */

#include <vector>
#include <string>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "common.h"
#include "utils.h"
#include "pairwise.h"
#include "search.h"

template <seqan3::writable_alphabet T>
int templated_main(modify_string_options const &options) {
    if (options.pairwise_filename != "" && options.type != "all" && options.type != "new") {
        // Running in pairwise mode
        return pairwise_matrix<T>(options);
    } else if (options.pairwise_filename != "" && options.type == "all")  {
        // Running in pairwise all mode
        return pairwise_all_matrix<T>(options);
    } else if (options.interleaved_filename != "") {
        // Running in interleaved mode
        return interleaved<T>(options);
    } else if (options.reference_filename != "" && options.query_filename != "") {
        return query_ref_search<T>(options);
    } else {
        seqan3::debug_stream << "Error: unknown configuration!" << std::endl;
        return 1;
    }
}

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

int main(int argc, char const ** argv) {
    // parse our options
    modify_string_options options;
    if (!parse_command_line(options, argc, argv)) {
        return 1;
    }

    // parse the mask so we know the kmer size
    options.effective_length = options.klen;
    if (!parse_mask(options, options.effective_length)) {
        return 1;
    }

    // dispatch the further work based on the sequence type
    if (options.sequence_type == "aa") {
        return templated_main<seqan3::aa27>(options);
    } else if (options.sequence_type == "raa") {
        return templated_main<seqan3::aa10murphy>(options);
    } else if (options.sequence_type == "dna") {
        return templated_main<seqan3::dna5>(options);
    } else {
        seqan3::debug_stream << "Error: sequence type not found (" << options.sequence_type << ")" << std::endl;
        return 1;
    }
}
