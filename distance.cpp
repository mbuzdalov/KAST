#include "distance.h"
#include <cstdint>
#include <cmath>

double d2s(std::vector<unsigned> const &kmerCounts1,
           std::vector<unsigned> const &kmerCounts2,
           std::vector<double> const &markovCounts1,
           std::vector<double> const &markovCounts2)
{

   uint64_t rN = 0;
   uint64_t qN = 0;
   double sum1 = 0.0;
   double sum2 = 0.0;
   double D2S = 0.0;

   for (unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      rN = rN + kmerCounts1[i];
      qN = qN + kmerCounts2[i];
   }

   for (unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      double qCt = kmerCounts2[i] - (qN * markovCounts2[i]);
      double rCt = kmerCounts1[i] - (rN * markovCounts1[i]);

      double dist = sqrt(qCt * qCt + rCt * rCt);

      if (dist == 0) dist = 1;
      D2S = D2S + (qCt * rCt / dist);
      sum1 = sum1 + (qCt * qCt / dist);
      sum2 = sum2 + (rCt * rCt / dist);
   }

   return 0.5 * (1 - D2S / (sqrt(sum1) * sqrt(sum2)));
}

double d2star(std::vector<unsigned> const &kmerCounts1,
              std::vector<unsigned> const &kmerCounts2,
              std::vector<double> const &markovCounts1,
              std::vector<double> const &markovCounts2)
{
   uint64_t rN = 0;
   uint64_t qN = 0;
   double tempQ = 0.0;
   double tempR = 0.0;
   double D_2Star = 0.0;

   for (unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      rN = rN + kmerCounts1[i];
      qN = qN + kmerCounts2[i];
   }

   for (unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      double temp1 = qN * markovCounts2[i];
      double temp2 = rN * markovCounts1[i];
      double cQi_bar = kmerCounts2[i] - temp1;
      double cRi_bar = kmerCounts1[i] - temp2;
      double temp3 = sqrt(temp1 * temp2);

      if (temp1 == 0)
      {
         temp1 = 1.0;
         temp3 = 1.0;
      }
      if (temp2 == 0)
      {
         temp2 = 1.0;
         temp3 = 1.0;
      }

      D_2Star += cQi_bar * cRi_bar / temp3;
      tempQ += cQi_bar * cQi_bar / temp1;
      tempR += cRi_bar * cRi_bar / temp2;
   }
   double temp = D_2Star / (sqrt(tempQ) * sqrt(tempR));
   return 0.5 * (1 - temp);
}

double dai(std::vector<unsigned> const &kmerCounts1,
           std::vector<unsigned> const &kmerCounts2,
           std::vector<double> const &markovCounts1,
           std::vector<double> const &markovCounts2)
{
   double S2kr = 0.0;
   double tempX = 0.0;
   double tempY = 0.0;

   uint64_t qN = 0;
   uint64_t rN = 0;
   uint64_t n = kmerCounts1.size();

   for (unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      rN = rN + kmerCounts1[i];
      qN = qN + kmerCounts2[i];
   }

   for (unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      double rC = kmerCounts1[i];
      double qC = kmerCounts2[i];
      double rP = markovCounts1[i];
      double qP = markovCounts2[i];

      // In the d2-tools implementation of this, for some reason the python code reads
      // a zero from the wordcount files as 1e-10. This has the effect of ensuring that
      // there are never any actual zeros in the calculation. I don't know if this is
      // intentional or accidental, but for now, to ensure this is consistent with
      // d2-tools I've reproduced this situation here.
      if (rC == 0)
         rC = rC + 0.0000000001;
      if (qC == 0)
         qC = qC + 0.0000000001;

      double r_sigma = (float) (rC) * (float) (rP);
      double q_sigma = (float) (qC) * (float) (qP);
      double temp1 = 0.0;
      double temp2 = 0.0;

      if (r_sigma == q_sigma)
      {
         temp1 = 0;
         temp2 = 0;
      }
      else
      {
         double temp = (float) (r_sigma + q_sigma);
         double temp3 = (2 * r_sigma) / (float) (temp);
         double temp4 = (2 * q_sigma) / (float) (temp);

         temp1 = r_sigma * log(temp3);
         temp2 = q_sigma * log(temp4);
      }
      tempX += temp1;
      tempY += temp2;
   }

   if (tempX == 0.0 && tempY == 0.0)
      S2kr = 0.0;
   else
      S2kr = (tempX + tempY) / n + 2 * log(2.0);

   return S2kr;
}

double hao(std::vector<unsigned> const &kmerCounts1,
           std::vector<unsigned> const &kmerCounts2,
           std::vector<double> const &markovCounts1,
           std::vector<double> const &markovCounts2)
{


}

double cosine(std::vector<unsigned> const &kmerCounts1,
              std::vector<unsigned> const &kmerCounts2)
{
   uint64_t sumqCrC = 0;
   uint64_t sumqC2 = 0;
   uint64_t sumrC2 = 0;

   for (unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      sumqCrC += (uint64_t) (kmerCounts1[i]) * (uint64_t) (kmerCounts2[i]);
      sumqC2  += (uint64_t) (kmerCounts1[i]) * (uint64_t) (kmerCounts1[i]);
      sumrC2  += (uint64_t) (kmerCounts2[i]) * (uint64_t) (kmerCounts2[i]);
   }

   double score = sumqCrC / (sqrt(sumqC2) * sqrt(sumrC2));
   return 1 - score;
}

double d2(std::vector<unsigned> const &kmerCounts1,
          std::vector<unsigned> const &kmerCounts2)
{
   return 0.5 * cosine(kmerCounts1, kmerCounts2); // TODO: correct?
}

double euler(std::vector<unsigned> const &kmerCounts1,
             std::vector<unsigned> const &kmerCounts2)
{
   double score = 0.0;
   uint64_t rN = 0;
   uint64_t qN = 0;

   for(unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      rN = rN + kmerCounts1[i];
      qN = qN + kmerCounts2[i];
   }

   for (unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      double rF = kmerCounts1[i] / (double) (rN);
      double qF = kmerCounts2[i] / (double) (qN);
      score = score + pow((rF - qF), 2); // TODO: (rF - qF) * (rF - qF)?
   }

   return pow(score, 0.5); // TODO: sqrt?
}

double bray_curtis_distance(std::vector<unsigned> const &kmerCounts1,
                            std::vector<unsigned> const &kmerCounts2)
{
   double sumMinus = 0.0;
   double sumPlus = 0.0;

   for (unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      sumMinus = sumMinus + std::abs((int64_t) (kmerCounts2[i]) - (int64_t) (kmerCounts1[i]));
      sumPlus = sumPlus + kmerCounts2[i] + kmerCounts1[i];
   }

   return sumMinus / sumPlus;
}

double normalised_google_distance(std::vector<unsigned> const &kmerCounts1,
                                  std::vector<unsigned> const &kmerCounts2)
{
   double score = 0.0;
   double sumqC = 0.0;
   double sumrC = 0.0;
   double sum_min_qr = 0.0;

   for (unsigned i = 0; i < kmerCounts1.size(); i++)
   {
      sumrC += kmerCounts1[i];
      sumqC += kmerCounts2[i];

      if (kmerCounts1[i] < kmerCounts2[i])
         sum_min_qr += kmerCounts1[i];
      else
         sum_min_qr += kmerCounts2[i];
   }

   double sum_max, sum_min;

   if (sumqC > sumrC)
   {
      sum_max = sumqC;
      sum_min = sumrC;
   }
   else
   {
      sum_max = sumrC;
      sum_min = sumqC;
   }

   double sum_all = sumqC + sumrC;

   return (sum_max - sum_min_qr) / (sum_all - sum_min);
}

double chebyshev(std::vector<unsigned> const &kmerCounts1,
                 std::vector<unsigned> const &kmerCounts2)
{
   double score = 0.0;
   uint64_t rN = 0;
   uint64_t qN = 0;
   double temp = 0.0;

   for (unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      rN = rN + kmerCounts1[i];
      qN = qN + kmerCounts2[i];
   }

   for(unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      double rF = kmerCounts1[i] / (double) (rN);
      double qF = kmerCounts2[i] / (double) (qN);
      temp = std::abs(rF - qF);
      if (temp > score)
         score = temp;
   }
   return score;
}

double canberra(std::vector<unsigned> const &kmerCounts1,
                std::vector<unsigned> const &kmerCounts2)
{
   double score = 0.0;
   for (unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      double p1 = std::abs((int64_t) (kmerCounts1[i]) - (int64_t) (kmerCounts2[i]));
      double p2 = (uint64_t) (kmerCounts1[i]) + kmerCounts2[i];
      if (p1 != 0 && p2 != 0)
         score += p1 / p2;
   }
   return score;
}

double normalised_canberra(std::vector<unsigned> const & kmerCounts1,
                           std::vector<unsigned> const & kmerCounts2)
{
   return canberra(kmerCounts1, kmerCounts2) / kmerCounts1.size();
}

double manhattan(std::vector<unsigned> const & kmerCounts1,
                 std::vector<unsigned> const & kmerCounts2)
{
   double score = 0.0;
   uint64_t rN = 0;
   uint64_t qN = 0;

   for (unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      rN = rN + kmerCounts1[i];
      qN = qN + kmerCounts2[i];
   }

   for (unsigned i = 0; i < kmerCounts1.size(); ++i)
   {
      double rF = kmerCounts1[i] / (double) (rN);
      double qF = kmerCounts2[i] / (double) (qN);
      score = score + std::abs(rF - qF);
   }

   return score;
}
