#ifndef __KAST_SEARCH_H__
#define __KAST_SEARCH_H__

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

#include "common.h"
#include "distance.h"
#include "utils.h"

template <seqan3::writable_alphabet T>
void print_result(modify_string_options const &options,
                  std::string const &queryid,
                  std::ostream &outfile,
                  std::vector<T> const &queryseq,
                  std::multimap<double, int> &results,
                  std::vector<std::string> const &referenceids,
                  std::vector<std::vector<T>> const &referenceseqs) {
    if (options.output_format == "tabular") {
        std::vector<std::string> split = split_by_space(queryid);
        std::string q_name = split[0];
        outfile << "############################ ";
        outfile << "\t" << q_name << std::endl;

        for (auto p : results) {
            std::vector<std::string> split2 = split_by_space(referenceids[p.second]);
            outfile << p.first << "\t" << referenceseqs[p.second].size();
            outfile << "\t" << split2[0] << std::endl;
        }
    } else if (options.output_format == "blastlike") {
        if (!options.noheader) {
            outfile << "RefID\tQryID\tRefLen\tQryLen\t";
            if (options.calcgc) {
                outfile << "RefGC\tQryGC\t";
            }
            outfile << "HitRank\tScore" << std::endl;
        }

        size_t count = 1;
        for (auto p : results) {
            outfile << referenceids[p.second] << "\t" << queryid << "\t";
            outfile << referenceseqs[p.second].size() << "\t";
            outfile << queryseq.size() << "\t";
            if (options.calcgc) {
               outfile << gc_ratio(referenceseqs[p.second]) << "\t";
               outfile << gc_ratio(queryseq) << "\t" << count << "\t";
            }
            outfile << count << "\t";
            outfile << p.first << std::endl;
            ++count;
        }
    } else {
        outfile << "############################ " << queryid << std::endl;

        for (auto p : results) {
            outfile << referenceids[p.second] << " " << p.first << std::endl;
        }
    }
}

template <seqan3::alphabet T>
struct mt_record_reader {
    std::mutex lock;
    sequence_file_type<T> &file;
    sequence_file_type<T>::iterator iterator;

    mt_record_reader(sequence_file_type<T> &file): file(file), iterator(file.begin()) {}

    bool try_read(std::string &id, std::vector<T> &sequence) {
        lock.lock();
        if (iterator == file.end()) {
            lock.unlock();
            return false;
        } else {
            auto const &record = *iterator;
            id = record.id();
            sequence = record.sequence();
            ++iterator;
            lock.unlock();
            return true;
        }
    }
};

template <seqan3::writable_alphabet T>
int search_thread(modify_string_options const &options,
                  mt_record_reader<T> &qryseq_reader,
                  std::vector<std::string> const &refids,
                  std::vector<std::vector<T>> const &refseqs,
                  std::vector<std::vector<unsigned>> const &refcounts,
                  std::vector<std::vector<double>> const &markovcounts,
                  std::mutex &write, std::ostream &outfile)
{
    std::string queryid;
    std::vector<T> queryseq;
    std::vector<unsigned> querycounts;
    std::vector<double> querymarkov;

    while (qryseq_reader.try_read(queryid, queryseq)) {
        if (!options.noreverse) {
            prepare_sequence_inplace(queryseq, options.klen);
        }

        populate_counts(options, queryseq, querycounts, querymarkov);

        // store the results
        std::multimap<double, int> results;

        for (size_t i = 0; i < refids.size(); ++i) {
            if (options.filter_bp != 0) {
                if (std::abs((int64_t) refseqs[i].size() - (int64_t) queryseq.size()) > options.filter_bp) {
                    continue;
                }
            } else if (options.filter_percent != 0) {
                double ratio = (double) (queryseq.size()) / (double) (refseqs[i].size());
                if (ratio > 1) {
                    ratio = 1 / ratio;
                }
                if (ratio <= options.filter_percent / 100.0) {
                    continue;
                }
            }

            double dist = distance_dispatch(options.type, querycounts, refcounts[i], querymarkov, markovcounts[i]);

            // stores the smallest distance results and corresponding location in refSeq
            /*
                Important:! It's at this point that any decision is made regarding how many
                results are stored.

                The new logic, now that we have a cutoff, is that if we have a cutoff set, it
                overides the num-hits flag.
            */

            if (std::isnan(options.score_cutoff) || dist <= options.score_cutoff) {
                results.insert(std::make_pair(dist, i));
            }

            if (std::isnan(options.score_cutoff) && options.nohits != 0 && results.size() > options.nohits) {
                auto it = results.end();
                results.erase(--it);
            }
        }

        write.lock();
        print_result<T>(options, queryid, outfile, queryseq, results, refids, refseqs);
        write.unlock();
    }

    return 0;
}

template <seqan3::writable_alphabet T>
int query_ref_search(modify_string_options const &options) {
    sequence_file_type<T> refseq_file { options.reference_filename };
    std::vector<std::string> refids;
    std::vector<std::vector<T>> refseqs;

    for (auto record : refseq_file) {
        refids.push_back(record.id());
        refseqs.push_back(record.sequence());
        if (!options.noreverse) {
            prepare_sequence_inplace(refseqs.back(), options.klen);
        }
    }

    // check how much RAM is required to store the reference
    // this check is actually quite useless since there are multiple copies running
    // so we need to replace it with something better
    // if(mem_check(options, refseqs.size(), alphabetType) == 1)
    //    return 1;

    std::vector<std::vector<unsigned>> counts;
    std::vector<std::vector<double>> markov_counts;

    counts.resize(refseqs.size());
    markov_counts.resize(refseqs.size()); // maybe do it lazily

    for (size_t i = 0; i < refseqs.size(); ++i) {
        std::vector<T> const &refseq = refseqs[i];
        populate_counts(options, refseqs[i], counts[i], markov_counts[i]);
    }

    // search thread
    // open up query file
    sequence_file_type<T> qryseq_file { options.query_filename };
    mt_record_reader<T> qryseq_reader(qryseq_file);

    // open up output file
    std::ofstream outfile;

    if (options.output_filename != "") {
        try {
            outfile.open(options.output_filename, std::ios_base::out);
        } catch (std::ofstream::failure const &e) {
            seqan3::debug_stream << "Error: could not open output file '" << options.output_filename << "'" << std::endl;;
            return 1;
        }
    }

    std::mutex write_lock;
    std::vector<std::thread> worker_threads;

    auto outfile_ref = options.output_filename == "" ? std::ref(std::cout) : std::ref(outfile);
    for (unsigned i = 0; i < options.num_threads; ++i) {
        worker_threads.emplace_back(search_thread<T>,
            std::ref(options), std::ref(qryseq_reader),
            std::ref(refids), std::ref(refseqs), std::ref(counts), std::ref(markov_counts),
            std::ref(write_lock), outfile_ref
        );
    }

    for (auto &th : worker_threads) {
        th.join();
    }

   return 0;
}

#endif // __KAST_SEARCH_H__
