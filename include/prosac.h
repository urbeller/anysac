#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <limits>
#include <cassert>
#include <numeric>
#include "utils.h"

// Taken from : https://github.com/RomanJuranek/librectify/blob/master/src/prosac.h


namespace AnySac
{

  static const float chi2_table[20] = {
    std::numeric_limits<float>::infinity(), 6.6348966 , 5.41189443, 4.70929225, 4.21788459,
    3.84145882, 3.5373846 , 3.28302029, 3.06490172, 2.8743734,
    2.70554345, 2.55422131, 2.41732093, 2.29250453, 2.17795916,
    2.07225086, 1.97422609, 1.88294329, 1.79762406, 1.71761761};

  static inline
    int niter_RANSAC(double p, // probability that at least one of the random samples picked up by RANSAC is free of outliers
        double epsilon, // proportion of outliers
        int s, // sample size
        int Nmax) // upper bound on the number of iterations (-1 means INT_MAX)
    {
      // compute safely N = ceil(log(1. - p) / log(1. - exp(log(1.-epsilon) * s)))
      double logarg, logval, N;
      if (Nmax == -1) 
      {
        Nmax = INT_MAX;
      }
      
      assert(Nmax >= 1);
      
      if (epsilon <= 0.) 
      {
        return 1;
      }
      
      // logarg = -(1-epsilon)^s
      logarg = -exp(s*log(1.-epsilon)); 
      logval = log(1.+logarg); 
      N = log(1.-p) / logval;
      if (logval  < 0. && N < Nmax) 
      {
        return (int)ceil(N);
      }

      return Nmax;
    }

  template< typename DriverType, typename ModelType>
    class Prosac
    {
      public:

        float eta {0.05};
        float beta {0.01};
        float psi {0.02};
        float p_good_sample {0.9};
        float max_outlier_proportion {0.5};

        Prosac(DriverType &driver) : m_driver(driver)
        {
          chi2_value = chi2(2*psi);          
        }

               
       

      private:
        std::mt19937 rng;
        float chi2_value;

        DriverType &m_driver;



        int Imin(int m, int n) 
        {
          double mu = n*beta;
          double sigma = sqrt(n*beta*(1-beta));
          return (int)ceil(m + mu + sigma*sqrt(chi2_value));
        }

        float chi2(float p)
        {
          int i = floor(clip(p, 0.01f, 0.2f) * 100);
          return chi2_table[i];
        }


      public:
        // Solve the problem.
        // Input: a model object (driver).
        // Returns: number of inliers.
        int solve(ModelType &bestModel)
        {
            const int N = m_driver.data_size();
            const int m = m_driver.model_size();
            std::vector<int> sample(m);

          if(N == m)
          {
              std::iota(sample.begin(), sample.end(), 0);
              bool status = m_driver.fit_model(sample, bestModel);
              return status;
          }
            
         
          const int T_N = niter_RANSAC(p_good_sample, max_outlier_proportion, m, -1);

          std::vector<unsigned char> inliersFlags(N);

          int n_star; // termination length (see sec. 2.2 Stopping criterion)
          int I_n_star; // number of inliers found within the first n_star data points
          int I_N_best; // best number of inliers found so far (store the model that goes with it)
          const int I_N_min = (1.-max_outlier_proportion)*N; // the minimum number of total inliers
          int t; // iteration number
          int n; // we draw samples from the set U_n of the top n data points
          double T_n; // average number of samples {M_i}_{i=1}^{T_N} that contain samples from U_n only
          int T_n_prime; // integer version of T_n, see eq. (4)
          int k_n_star; // number of samples to draw to reach the maximality constraint
          int i;
          const double logeta0 = log(eta);

          n_star = N;
          I_n_star = 0;
          I_N_best = 0;
          t = 0;
          n = m;
          T_n = T_N;
          
          for(i = 0; i < m; i++) 
          {
            T_n *= (double)(n-i)/(N-i);
          }

          T_n_prime = 1;
          k_n_star = T_N;


          while(((I_N_best < I_N_min) || t <= k_n_star) && t < T_N)
          {
            // Choice of the hypothesis generation set
            t = t + 1;

            if ((t > T_n_prime) && (n < n_star)) 
            {
              double T_nplus1 = (T_n * (n+1)) / (n+1-m);
              n = n+1;
              T_n_prime = T_n_prime + ceil(T_nplus1 - T_n);
              T_n = T_nplus1;
            }
           

            // Draw semi-random sample (note that the test condition from Algorithm1 in the paper is reversed):
            if (t > T_n_prime) 
            {
              AnySac::choiceKnuth(n, m, rng, sample);
            }
            else 
            {
              AnySac::choiceKnuth(n-1, m-1, rng, sample);
              sample[m-1] = n;
            }



            if (!m_driver.check_subset(sample))
            {
              continue;
            }

            ModelType model;
            bool ok = m_driver.fit_model(sample, model);
            if(!ok)
              continue;

              float inliersErrors = 0;
            int I_N = m_driver.count_inliers(model, inliersFlags, inliersErrors);


            if (I_N > I_N_best) 
            {
              int n_best; // best value found so far in terms of inliers ratio
              int I_n_best; // number of inliers for n_best


              I_N_best = I_N;
              bestModel = model;
               // printf("--> total = %f avg=%f  n=%d\n", inliersErrors, inliersErrors / float(I_N_best), I_N_best);

              // INSERT: Store the best model
              //p_best = p_t;
              //best_score = err;

              n_best = N;
              I_n_best = I_N;

              int n_test; // test value for the termination length
              int I_n_test; // number of inliers for that test value
              double epsilon_n_best = (double)I_n_best/n_best;

              for(n_test = N, I_n_test = I_N; n_test > m; n_test--) { 
                assert(n_test >= I_n_test);

                if (( I_n_test * n_best > I_n_best * n_test ) &&
                    ( I_n_test > epsilon_n_best*n_test + sqrt(n_test*epsilon_n_best*(1.-epsilon_n_best)*2.706) )) 
                {
                  if (I_n_test < Imin(m,n_test)) 
                  {
                    break; // jump out of the for(n_test) loop
                  }

                  n_best = n_test;
                  I_n_best = I_n_test;
                  epsilon_n_best = (double)I_n_best/n_best;
                }

                I_n_test -= (int)inliersFlags[n_test-1];
              } // for(n_test ...

              // is the best one we found even better than n_star?
              if ( I_n_best * n_star > I_n_star * n_best ) 
              {
                double logarg;
                assert(n_best >= I_n_best);
                // update all values
                n_star = n_best;
                I_n_star = I_n_best;
                k_n_star = niter_RANSAC(1.-eta, 1.-I_n_star/(double)n_star, m, T_N);
              }
              } // if (I_N > I_N_best)
            } // while(t <= k_n_star ...

            return I_N_best;
        }

    };
}

