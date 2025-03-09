#ifndef __KAST_COMMON_H__
#define __KAST_COMMON_H__

#include <string>
#include <vector>

/*
User defined options struct
*/
struct modify_string_options
{
   int klen;
   int nohits;
   int markov_order;
   std::string type;
   std::string sequence_type;
   std::string output_format;
   bool noreverse;
   bool calcgc;
   bool noheader;
   std::string query_filename;
   std::string reference_filename;
   std::string pairwise_filename;
   std::string interleaved_filename;
   int num_threads;
   bool debug;
   bool lowram;
   bool phylyp = true;
   bool tabout;
   bool blastlike;
   std::string output_filename;
   std::vector<std::string> mask;
   int effective_length;
   double score_cutoff;
   double filter_percent = 0;
   int64_t filter_bp = 0;
};

#endif // __KAST_COMMON_H__
