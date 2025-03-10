#include <fstream>
#include <map>
#include <mutex>
#include <ostream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include "utils.h"
#include "distance.h"

template <seqan3::writable_alphabet T>
int calc_distance(unsigned &rI, unsigned &cI,
                  std::vector<std::vector<double>> &results,
                  std::mutex &location,
                  std::vector<std::vector<unsigned>> const &counts,
                  std::vector<std::vector<double>> const &markov_counts,
                  modify_string_options const &options) {
    while (true) {
        location.lock();
        unsigned row = rI;
        unsigned column = cI;

        // 0 <= row < column <= counts.size()

        if (cI == counts.size()) {
            location.unlock();
            return 0;
        }

        if (++cI == counts.size()) {
            ++rI;
            cI = rI + 1;
        }

        location.unlock();

        double dist = distance_dispatch(options.type, counts[row], counts[column], markov_counts[row], markov_counts[column]);

        results[row][column] = dist;
        results[column][row] = dist;
    }
};

template <seqan3::writable_alphabet T>
int pairwise_matrix(modify_string_options const &options) {
    sequence_file_type<T> pw_file { options.pairwise_filename };
    std::vector<std::string> pwids;
    std::vector<std::vector<T>> pwseqs;

    for (auto record : pw_file) {
        pwids.push_back(record.id());
        pwseqs.push_back(record.sequence());
        if (!options.noreverse) {
            prepare_sequence_inplace(pwseqs.back(), options.klen);
        }
    }

    // set up elements for thread
    std::vector<std::thread> thread_pool;
    std::mutex location;
    unsigned rI = 0;
    unsigned cI = 0;
    unsigned int cores = options.num_threads;

    // store the distance calculation results
    std::vector<std::vector<double>> results(pwseqs.size(), std::vector<double>(pwseqs.size(), 0.0));
    std::vector<std::vector<unsigned>> counts(pwseqs.size());
    std::vector<std::vector<double>> markov_counts(pwseqs.size());

    for (size_t i = 0; i < pwseqs.size(); ++i) {
        populate_counts(options, pwseqs[i], counts[i], markov_counts[i]);
    }

    for (size_t i = 0; i < cores; ++i) {
        thread_pool.emplace_back(calc_distance<T>,
            std::ref(rI), std::ref(cI), std::ref(results), std::ref(location),
            std::ref(counts), std::ref(markov_counts), std::ref(options)
        );
    }

    for (auto &thread : thread_pool) {
        thread.join();
    }

    // Print results
    std::ofstream out_file;
    if (options.output_filename != "") {
        try {
            out_file.open(options.output_filename, std::ios_base::out);
        } catch (std::ofstream::failure const &e) {
            seqan3::debug_stream << "Error: could not open output file " << options.output_filename << std::endl;
            return 1;
        }
    }
    std::ostream &out = options.output_filename == "" ? std::cout : out_file;

    out << pwseqs.size() << std::endl;
    for (size_t i = 0; i < pwseqs.size(); ++i) {
        std::vector<std::string> split = split_by_space(pwids[i]);
        unsigned cut_size = 10;
        split[0].resize(cut_size, ' ');
        out << split[0] << "\t";
        for (size_t j = 0; j < pwseqs.size(); ++j) {
            out << results[i][j] << "\t";
        }
        out << std::endl;
    }

    return 0;
};

template <seqan3::writable_alphabet T>
int interleaved(modify_string_options const &options) {
    sequence_file_type<T> il_file { options.interleaved_filename };
    std::vector<std::string> ids;
    std::vector<std::vector<T>> seqs;
    std::vector<std::vector<unsigned>> counts(2);
    std::vector<std::vector<double>> markov_counts(2);

    // TODO: could be multithreaded
    for (auto record : il_file) {
        ids.push_back(record.id());
        seqs.push_back(record.sequence());
        if (!options.noreverse) {
            prepare_sequence_inplace(seqs.back(), options.klen);
        }
        size_t i = seqs.size() - 1;
        populate_counts(options, seqs[i], counts[i], markov_counts[i]);
        if (ids.size() == 2) {
            double dist = distance_dispatch(options.type, counts[0], counts[1], markov_counts[0], markov_counts[1]);
            std::cout << ids[0] << "\t" << ids[1] << "\t" << options.type << "\t" << dist << std::endl;
            ids.clear();
            seqs.clear();
        }
    }

    if (ids.size() != 0) {
        seqan3::debug_stream << "WARNING: " << options.interleaved_filename
                             << " contains an unequal number of sequences." << std::endl;
    }

    return 0;
}

template <seqan3::writable_alphabet T>
int pairwise_all_matrix(modify_string_options const &options) {
    sequence_file_type<T> pw_file { options.pairwise_filename };
    std::vector<std::string> pwids;
    std::vector<std::vector<T>> pwseqs;

    for (auto record : pw_file) {
        pwids.push_back(record.id());
        pwseqs.push_back(record.sequence());
        if (!options.noreverse) {
            prepare_sequence_inplace(pwseqs.back(), options.klen);
        }
    }

    std::cout << "Q1\tQ2\tEuclid\td2\tcosine\tManhattan\tBC\tNGD\tHao\tdai\tD2S\tD2Star\t"
                 "Chebyshev\tCanberra\tNormalised Canberra\n";

    std::vector<unsigned> counts_i, counts_j;
    std::vector<double> markov_i, markov_j;
    for (size_t i = 0; i < pwids.size(); ++i) {
        populate_counts(options, pwseqs[i], counts_i, markov_i);
        for (size_t j = 0; j < pwids.size(); ++j) {
            populate_counts(options, pwseqs[j], counts_j, markov_j);

            std::cout << pwids[i] << "\t" << pwids[j]
                      << "\t" << euler(counts_i, counts_j)
                      << "\t" << d2(counts_i, counts_j)
                      << "\t" << cosine(counts_i, counts_j)
                      << "\t" << manhattan(counts_i, counts_j)
                      << "\t" << bray_curtis_distance(counts_i, counts_j)
                      << "\t" << normalised_google_distance(counts_i, counts_j)
                      << "\t" << hao(counts_i, counts_j, markov_i, markov_j)
                      << "\t" << dai(counts_i, counts_j, markov_i, markov_j)
                      << "\t" << d2s(counts_i, counts_j, markov_i, markov_j)
                      << "\t" << d2star(counts_i, counts_j, markov_i, markov_j)
                      << "\t" << chebyshev(counts_i, counts_j)
                      << "\t" << canberra(counts_i, counts_j)
                      << "\t" << normalised_canberra(counts_i, counts_j)
                      << std::endl;
        }
    }
    return 0;
};
