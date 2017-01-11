////////////////////////////////////////////////////////////////////// 
// regress/FancyRegression.cpp 
// (c) 2000-2007 Goncalo Abecasis
// 
// This file is distributed as part of the MERLIN source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile MERLIN.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Tuesday December 18, 2007
// 
 
#include "FancyRegression.h"
#include "MerlinCore.h"
#include "MathStats.h"

#include <math.h>

FancyRegression::FancyRegression() :
   sum("sum_sq"), diff("sum_diff"), y("Y"),
   weights("pi_weights"), pi("pi"), E("E"),
   sigma_ss("sigma_ss"), sigma_dd("sigma_dd"), sigma_sd("sigma_sd"),
   sigma_ds("sigma_ds"), sigma_yy("sigma_yy"), B("B"), H("H"),
   sigma_pi("sigma_pi")
   {
   trait = -1;
   mean = 0.0;
   variance = 1.0;
   heritability = 0.5;
   testRetestCorrel = 1.0;

   // no heterogeneity by default
   // covariate = false;
   // covarType = covarId = 0;
   // covarThreshold = 0.0;

   // Bound estimates of regression parameters
   bounded = true;
   }

void FancyRegression::SetPositionCount(int count)
   {
   scores.Dimension(count);
   scores.Zero();

   information.Dimension(count);
   information.Zero();
   }

void FancyRegression::SetupFamily(RegressKinship & kin)
   {
   Pedigree & ped = *kin.GetPedigree();
   Family * f = kin.GetFamily();

   // Clear data structures
   sum.Clear(); diff.Clear(); y.Clear(); pi.Clear();
   pairs1.Clear(); pairs2.Clear(); weights.Clear();
   repeatCounts.Clear();

   int phenotypes = 0;

   // Check if we have repeated measures for this trait
   int repeats = ped.LookupCovariate(ped.traitNames[trait] + "_repeats");

   // For each non-founder in the family
   for (int j = 0; j < f->count; j++)
      {
      Person & p1 = ped[f->path[j]];

      // Keep track of the number of repeated measurements for each individual
      if (repeats >= 0)
         repeatCounts.Push((int) (p1.covariates[repeats]));
      else
         repeatCounts.Push(1);

      // If not phenotyped, ignore
      if (!p1.isPhenotyped(trait)) continue;

      phenotypes++;

      double trait1 = (p1.traits[trait] - mean) / sqrt(variance);

      // Pair person i with i+1 .. n
      for (int k = j + 1; k < f->count; k++)
         {
         Person & p2 = ped[f->path[k]];

         if (!p2.isPhenotyped(trait)) continue;

         double trait2 = (p2.traits[trait] - mean) / sqrt(variance);
         double kinship = kin.Retrieve(j, k);
         double r = 2 * kinship * heritability;

         // Keep track of squared trait sums and differences
         sum.Push(square(trait1 + trait2) - 2 * (1 + r));
         diff.Push(square(trait1 - trait2) - 2 * (1 - r));

         // Keep track of who is who
         pairs1.Push(j);
         pairs2.Push(k);

         // Keep track of subgroups for testing heterogeneity
         // if (CovariateObserved(p1) && CovariateObserved(p2))
         //   subgroup[Classify(p1) + Classify(p2)].Push(pairs1.Length());
         }
      }

   // Nothing to do for families with no informative pairs
   if (pairs1.Length() == 0)
      {
      maxInfo = 0;
      return;
      }

   // To make sure the model is identified, we only retain one
   // squared difference per phenotyped individual. We want to
   // to have each individual represented at least once, and the
   // following keeps pairs of the first individual with all
   // others plus a single pairing of phenotyped individuals
   // two and three
   //
   diff.Dimension(phenotypes == 2 ? 1 : phenotypes);

   // Stack sums and differences into a single vector
   y = sum;
   y.Stack(diff);

   // Calculate variance covariance matrix for squared sums
   // Note: the function correl(x,y) returns the correl if
   // x <> y and 1.0 if x == y.
   //
   sigma_ss.Dimension(sum.Length(), sum.Length());

   for (int j = 0; j < sum.Length(); j++)
      for (int k = 0; k < sum.Length(); k++)
         sigma_ss[j][k] = 2 * square
            (calc_correl(kin, pairs1[j], pairs1[k]) +
             calc_correl(kin, pairs1[j], pairs2[k]) +
             calc_correl(kin, pairs2[j], pairs1[k]) +
             calc_correl(kin, pairs2[j], pairs2[k]));

   // Calculate variance covariance matrix for squared differences
   //
   sigma_dd.Dimension(diff.Length(), diff.Length());

   for (int j = 0; j < diff.Length(); j++)
      for (int k = 0; k < diff.Length(); k++)
         sigma_dd[j][k] = 2 * square
            (calc_correl(kin, pairs1[j], pairs1[k]) +
             calc_correl(kin, pairs2[j], pairs2[k]) -
             calc_correl(kin, pairs2[j], pairs1[k]) -
             calc_correl(kin, pairs1[j], pairs2[k]));

   // Calculate the variance-covariance matrix for sums vs. differences
   //
   sigma_sd.Dimension(sum.Length(), diff.Length());

   for (int j = 0; j < sum.Length(); j++)
      for (int k = 0; k < diff.Length(); k++)
         sigma_sd[j][k] = 2 * square
            (calc_correl(kin, pairs1[j], pairs1[k]) +
             calc_correl(kin, pairs2[j], pairs1[k]) -
             calc_correl(kin, pairs1[j], pairs2[k]) -
             calc_correl(kin, pairs2[j], pairs2[k]));

   // Transpose it to get differences vs. sum matrix
   //
   sigma_ds.Transpose(sigma_sd);

   // Now we need to stack matrices to get sigma_yy with the
   // following pattern:
   //
   //   +---------------+---------------+
   //   |    sigma_ss   |  sigma_sd     |
   //   +---------------+---------------+
   //   |    sigma_ds   |  sigma_dd     |
   //   +---------------+---------------+
   //

   sigma_yy = sigma_ss;
   sigma_yy.StackLeft(sigma_sd);

   sigma_ds.StackLeft(sigma_dd);
   sigma_yy.StackBottom(sigma_ds);

   // Matrix H has one row per pairs and one column for
   // each squared sum and difference.
   H.Dimension(pairs1.Length(), sum.Length() + diff.Length());
   H.Zero();
   for (int j = 0; j < pairs1.Length(); j++)
      {
      if (j < sum.Length()) H[j][j] = 4.0;
      if (j < diff.Length()) H[j][j + sum.Length()] = -4.0;
      }

   // Matrix sigma_pi is initialized with variance-covariance matrix of
   // pi under the null
   //
   sigma_pi.Dimension(pairs1.Length(), pairs1.Length());
   pi.Dimension(pairs1.Length());
   for (int i = 0; i < pairs1.Length(); i++)
      {
      int person1 = pairs1[i];
      int person2 = pairs2[i];

      sigma_pi[i][i] =
         kin.Retrieve(person1, person2, person1, person2) -
         kin.Retrieve(person1, person2) * kin.Retrieve(person1, person2);

      pi[i] = kin.Retrieve(person1, person2);

      for (int j = i + 1; j < pairs1.Length(); j++)
         {
         int person3 = pairs1[j];
         int person4 = pairs2[j];

         sigma_pi[i][j] = sigma_pi[j][i] =
            kin.Retrieve(person1, person2, person3, person4) -
            kin.Retrieve(person1, person2) * kin.Retrieve(person3, person4);
         }
      }

   // Calculate the matrix inverse ...
   // It probably would be faster to back-substitute rows of B directly
   // using the cholesky decomposition (which should be LU decomposition
   // for non-positive definite matrices anyway!)
   //
   chol.Decompose(sigma_yy);
   chol.Invert();

   temp.Product(H, chol.inv);

   E.Product(temp, y);

   // Dimension unitialized arrays
   local_sigma_pi.Dimension(pairs1.Length(), pairs1.Length());
   local_pi.Dimension(pairs1.Length());

   // Update the maximum possible informativeness
   intermediate.Product(sigma_pi, E);

   maxInfo = E.InnerProduct(intermediate);
   sumMaxInfo += maxInfo;
   }

double FancyRegression::calc_correl(RegressKinship & kin, int person1, int person2)
   {
   // Calculate correlation between two measurements standardized
   // according to population mean and variance, but ignoring
   // measurement error
   return (person1 == person2) ?
      // Variance for each individual depends on the number of
      // repeated measurements and the error variance for each
      // measurement
      (testRetestCorrel + (1.0 - testRetestCorrel)  / repeatCounts[person1]) :
      // Covariance between two individuals just depends on the
      // heritability and kinship
      2.0 * kin.Retrieve(person1, person2) * heritability;
   }

void FancyRegression::Analyse(RegressKinship & kin, int position)
   {
   // Only analyse informative families
   if (pairs1.Length() == 0)
      return;

   // Calculate pi-hat and its variance at marker location
   for (int i = 0; i < pairs1.Length(); i++)
      {
      int person1 = pairs1[i];
      int person2 = pairs2[i];

      local_sigma_pi[i][i] =
         sigma_pi[i][i] -
        (kin.Retrieve(person1, person2, person1, person2) -
         kin.Retrieve(person1, person2) * kin.Retrieve(person1, person2));

      local_pi[i] = kin.Retrieve(person1, person2) - pi[i];

      for (int j = i + 1; j < pairs1.Length(); j++)
         {
         int person3 = pairs1[j];
         int person4 = pairs2[j];

         local_sigma_pi[i][j] = local_sigma_pi[j][i] =
            sigma_pi[i][j] -
           (kin.Retrieve(person1, person2, person3, person4) -
            kin.Retrieve(person1, person2) * kin.Retrieve(person3, person4));
         }
      }

   // Wrap up calculations!
   intermediate.Product(local_sigma_pi, E);

   double info = E.InnerProduct(intermediate);

   if (fabs(info) < 1E-10) return;

   double var_a2 = 1.0 / info;
   double a2 = E.InnerProduct(local_pi) * var_a2;

   scores[position] += a2 * info;
   information[position] += info;
   }

void  FancyRegression::PrintScores(FILE * tablefile, int chr, MerlinPDF & pdf, StringArray & labels)
   {
   printf("Pedigree-Wide Regression Analysis (%s)\n"
          "======================================================\n"
          "%15s %7s %7s %7s %7s %7s\n",
          (const char *) label, "Position", "H2", "Stdev", "Info", "LOD", "pvalue");

   if (pdf.isOpen())
      {
      pdf.PrepareChart();
      pdf.chart.yAxis.SetMin(0.0);
      pdf.chart.yAxis.minMax = 1.0;
      pdf.chart.yAxis.label = "LOD score";
      pdf.chart.title.printf("Regression Analysis for %s", (const char *) label);
      }

   for (int i = 0; i < scores.Length(); i++)
      if (information[i] > 1e-4)
         {
         double h2     = scores[i] / information[i];
         double var    = 1.0 / information[i];
         double chisq  = scores[i] * h2;

         if (bounded && h2 < 0.0)
            h2 = chisq = 0.0;

         if (bounded && h2 > 1.0)
            h2 = 1.0, chisq = information[i];

         double pvalue =(h2 > 0.0) ? chidist(chisq, 1) * 0.5 :
                                     1.0 - chidist(chisq, 1) * 0.5;
         double LOD    = chisq * 0.2171472409516 * (h2 >= 0.0 ? 1.0 : -1.0);

         int digits  = pvalue < 1.5e-4 ? 5 : int(log(pvalue*0.06666666)*-0.4343);

         printf("%15s %7.3f %7.3f %6.1f%% %7.3f %7.*f\n",
                (const char *) labels[i], h2, sqrt(var),
                information[i] / sumMaxInfo * 100., LOD, digits, pvalue);

         if (pdf.isOpen())
            pdf.y[0][i] = LOD;

         if (tablefile != NULL)
            fprintf(tablefile, "%d\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.4g\n",
                    chr, (const char *) labels[i], (const char *) label,
                    h2, sqrt(var), information[i] / sumMaxInfo * 100.,
                    LOD, pvalue);
         }
      else
         {
         printf("%15s %7s %7s %7s %7s %7s\n",
                (const char *) labels[i], "na", "na", "na", "na", "na");

         if (tablefile != NULL)
            fprintf(tablefile, "%d\t%s\t%s\tna\tna\tna\tna\tna\n",
                    chr, (const char *) labels[i], (const char *) label);
         }

   if (pdf.isOpen())
      pdf.DrawChart();

   printf("\n");
   }

bool FancyRegression::CovariateObserved(Person & p)
   {
   switch (covarType)
      {
      case COVAR_SEX :
         return p.sex == 1 || p.sex == 2;
      case COVAR_USER :
         return p.covariates[covarId] != _NAN_;
      case COVAR_TRAIT :
         return p.traits[covarId] != _NAN_;
      case COVAR_AFFECTION :
         return p.affections[covarId] == 1 || p.affections[covarId] == 2;
      }

   // Keeps the compiler happy!
   return -1;
   }

int FancyRegression::Classify(Person & p)
   {
   switch (covarType)
      {
      case COVAR_SEX :
         return p.sex - 1;
      case COVAR_USER :
         return p.covariates[covarId] > covarThreshold;
      case COVAR_TRAIT :
         return p.traits[covarId] > covarThreshold;
      case COVAR_AFFECTION :
         return p.affections[covarId] - 1;
      }

   // Keeps the compiler happy!
   return -1;
   }


 
