/*
 * KAST - Kmer Alignment-free Search Tool
 * Version 1.0.3beta
 * Original version by Dr. Martin Vickers (martin.vickers@jic.ac.uk)
 * Current version rewritten for seqan3 by Dr. Maxim Buzdalov (mab168@aber.ac.uk)
 */

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
