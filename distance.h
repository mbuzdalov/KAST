#ifndef __KAST_DISTANCE_H__
#define __KAST_DISTANCE_H__

#include <vector>

double d2(std::vector<unsigned> const &kmerCounts1,
          std::vector<unsigned> const &kmerCounts2);

double cosine(std::vector<unsigned> const &kmerCounts1,
              std::vector<unsigned> const &kmerCounts2);

double euler(std::vector<unsigned> const &kmerCounts1,
             std::vector<unsigned> const &kmerCounts2);

double bray_curtis_distance(std::vector<unsigned> const &kmerCounts1,
                            std::vector<unsigned> const &kmerCounts2);

double normalised_google_distance(std::vector<unsigned> const &kmerCounts1,
                                  std::vector<unsigned> const &kmerCounts2);

double chebyshev(std::vector<unsigned> const &kmerCounts1,
                 std::vector<unsigned> const &kmerCounts2);

double canberra(std::vector<unsigned> const &kmerCounts1,
                std::vector<unsigned> const &kmerCounts2);

double normalised_canberra(std::vector<unsigned> const &kmerCounts1,
                           std::vector<unsigned> const &kmerCounts2);

double manhattan(std::vector<unsigned> const &kmerCounts1,
                 std::vector<unsigned> const &kmerCounts2);

double d2s(std::vector<unsigned> const &kmerCounts1,
           std::vector<unsigned> const &kmerCounts2,
           std::vector<double> const &markovCounts1,
           std::vector<double> const &markovCounts2);

double d2star(std::vector<unsigned> const &kmerCounts1,
              std::vector<unsigned> const &kmerCounts2,
              std::vector<double> const &markovCounts1,
              std::vector<double> const &markovCounts2);

double dai(std::vector<unsigned> const &kmerCounts1,
           std::vector<unsigned> const &kmerCounts2,
           std::vector<double> const &markovCounts1,
           std::vector<double> const &markovCounts2);

double hao(std::vector<unsigned> const &kmerCounts1,
           std::vector<unsigned> const &kmerCounts2,
           std::vector<double> const &markovCounts1,
           std::vector<double> const &markovCounts2);

#endif
