#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <limits>
#include <cassert>
#include <numeric>
#include "anysac_utils.h"

// Taken from : https://github.com/RomanJuranek/librectify/blob/master/src/prosac.h

#define PROSAC_DISABLE_N_STAR_OPTIMIZATION

namespace AnySac
{

  static const float chi2_table[20] = {
    std::numeric_limits<float>::infinity(), 6.6348966 , 5.41189443, 4.70929225, 4.21788459,
    3.84145882, 3.5373846 , 3.28302029, 3.06490172, 2.8743734,
    2.70554345, 2.55422131, 2.41732093, 2.29250453, 2.17795916,
    2.07225086, 1.97422609, 1.88294329, 1.79762406, 1.71761761};

  template< typename DriverType, typename ModelType>
    class Prosac
    {
      public:

        int maxIter {200000};
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
        float chi2_value;

        DriverType &m_driver;
        std::mt19937 m_rgn;


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

        inline int uniformInt(int n)
        {
            return std::uniform_int_distribution<>(0, n) (m_rgn);
        }
        
        inline void deal_k_near_n(int n, int k, std::vector<int> &a)
        {
            assert(k <= (int)a.size());

            /* Warning: modifies k and n */
            /* Algorithm: go though all candidates from n-1 to 0, and pick each one with probability k/n */
            while((n > k) && (k > n/2)) {
                /* each number has probability k/n of being picked up */
                if (k > uniformInt(n-1)) {
                    /* pick this one up */
                    k--;
                    n--;
                    a[k]= n;
                } else {
                    /* don't pick this one */
                    n--;
                }
            }
            if (n == k) {
                /* we've got k numbers to sample from a set of k, easy... */
                for(n=n-1; n>=0; n--) {
                    a[n] = n;
                }
                k = 0;
            }
            if (k > 0) {
                assert(k <= n/2);
                ranksb(n, k, a);                /* reduced to ranksb */
            }
        }

        void ranksb(int n, int k, std::vector<int> &a)
        {
            assert(k <= (int)a.size());
            int i, l, m = 0, p, r = 0, s, x, m0 = 0, ds;
            

            if (k == 0)
                return;
            if (k == 1) 
            {
                a[0] = uniformInt(n-1);
                return;
            }
            
            
            /* Partition [0 : n-1] into k intervals:
             Store the least element of the I'th interval in a[i] */
            for (i = 0; i < k; ++i) 
            {
                a[i] = i * n / k;
            }
            
            /* Using a uniformly distributed random variable in the
             range 0 <= x < n, make k selections of x such that
             x lands in the remaining portion of some interval.
             At each successful selection, reduce the remaining
             portion of the interval by 1. */
            /* If k is close to n (say, bigger than n/2), the
             while loop may take many iterations. For this reason,
             it is better to use another algorithm in these
             situations (see deal_k_near_n() below). */
            for(i = 0; i < k;  ++i) {
                do {
                    x = uniformInt(n-1) + 1;
                    l = (x * k - 1) / n;
                    // l is the choosen index in a.
                    // a vector has size k.
                    // x <= n ensures that l < k:
                    // x == n:
                    //   l = (n*k-1)/n < n*k/n (=k)
                    // x == n+1:
                    //   l = ((n+1) * k - 1 ) /n
                    //     = (n*k + (k-1))/n >= n*k/n (=k)
                    assert(l<k);
                } while ((std::size_t)x <= a[l]);
                ++a[l];
            }
            
            /* Collect the least elements of any interval which
             absorbed a selection in the previous step into the
             low-order indices of a. */
            p = -1;
            for (i = 0; i < k; ++i) 
            {
                m = a[i];
                a[i] = 0;
                if (m != i * n / k) 
                {
                    /* A non-empty partition */
                    ++p;
                    a[p] = m;
                }
            }
            /* Allocate space for each non-empty partition starting
             from the high-order indices.  At the last position
             in each partition's segment, store the interval index
             of the partitions's successor. */
            s = k-1;
            for(; p >=0; p--) 
            {
                l = (a[p] * k - 1) / n;
                ds = a[p] - l * n / k;
                a[p] = 0;
                a[s] = l + 1;
                s -= ds;
            }
            
            for(l = k-1; l >= 0; l--) 
            {
                /* ranksb each of the sub-problems */
                x = a[l];
                if (x != 0) 
                {
                    /* Start a new bin */
                    r = l;
                    m0 = (x - 1) * n / k;
                    m = x * n / k - m0;
                    /* The order of arithmetic operations is important!
                     The same rounding errors must be produced in each
                     computation of a boundary. */
                }
                
                /* m0 is the least element of the current (l'th)
                 interval.  m is the count of the number of
                 unselected members of the interval. */
                x = m0 + uniformInt(m-1);
                assert(x >= 0 && x < n);
                
                /* Bubble Merge the (x-base_l)'th unselected member
                 of the current interval into the current interval's
                 segment (a [l..r]). */
                
                i = l;
                while (i < r && (std::size_t)x >= a[i+1]) 
                {
                    a[i] = a[i+1];
                    ++x;
                    ++i;
                }
                assert(x >= 0 && x < n);
                a[i] = x;
                --m;
            }
        }
        
        inline void deal(int n, int k, std::vector<int> &a)
        {
            assert(k <= n);
            assert(k <= (int)a.size());
            if (k <= n/2) {
                ranksb(n, k, a);
            } else {
                deal_k_near_n(n, k, a);
            }
        }
        
        int Prosac_Imin(int m, int n, double beta)
        {
            assert(beta > 0. && beta < 1.);
            const double mu = n * beta;
            const double sigma = std::sqrt(n*beta*(1-beta));
            // Imin(n) (equation (8) can then be obtained with the Chi-squared test with P=2*psi=0.10 (Chi2=2.706)
            return (int)std::ceil(m + mu + sigma * std::sqrt(2.706));
        }
        
      public:
        // Solve the problem.
        // Input: a model object (driver).
        // Returns: number of inliers.
        int solve(ModelType &bestModel, std::vector<unsigned char> &mask)
        {
            const int N = m_driver.data_size();
            const int N_draw = N;
            const int m = m_driver.model_size();
            std::vector<int> sample(m);

            std::random_device rd;
            m_rgn = std::mt19937(rd());
            
          if(N < m)
          {
              printf("%s::%d : No enough data points.\n", __FUNCTION__, __LINE__);
              return 0;
          }
            
          if(N == m)
          {
              std::iota(sample.begin(), sample.end(), 0);
              bool status = m_driver.fit_model(sample, bestModel);
              return status;
          }
            
          assert(beta > 0. && beta < 1.);
         
          const int T_N = (p_good_sample >= 1.0) ?  std::numeric_limits<int>::max() : niter_RANSAC(p_good_sample, max_outlier_proportion, m, -1);
          const int t_max = maxIter > 0 ? maxIter : T_N;
          
            int n_star = N; // termination length (see sec. 2.2 Stopping criterion)
            int I_n_star = 0; // number of inliers found within the first n_star data points
            int I_N_best = 0; // best number of inliers found so far (store the model that goes with it)
            const int I_N_min = (1. - max_outlier_proportion) * N; // the minimum number of total inliers
            int t = 0; // iteration number
            int n = m; // we draw samples from the set U_n of the top n data points
            double T_n = T_N; // average number of samples {M_i}_{i=1}^{T_N} that contain samples from U_n only
            int T_n_prime = 1; // integer version of T_n, see eq. (4)

            for(int i = 0; i < m; ++i) 
            {
                T_n *= (double)(n - i) / (N - i);
            }
            
            int k_n_star = T_N; // number of samples to draw to reach the maximality constraint
            
            std::vector<unsigned char> inliersFlags(N);

            while (((I_N_best < I_N_min) || t <= k_n_star) && t < T_N && t <= t_max) 
            {
                int I_N; // total number of inliers for that sample
                
                // Choice of the hypothesis generation set
                t = t + 1;
                
                // from the paper, eq. (5) (not Algorithm1):
                // "The growth function is then deﬁned as
                //  g(t) = min {n : T′n ≥ t}"
                // Thus n should be incremented if t > T'n, not if t = T'n as written in the algorithm 1
                if ((t > T_n_prime) && (n < n_star)) 
                {
                    double T_nplus1 = (T_n * (n+1)) / (n+1-m);
                    n = n+1;
                    T_n_prime = T_n_prime + std::ceil(T_nplus1 - T_n);
                    T_n = T_nplus1;
                }
                
                // Draw semi-random sample (note that the test condition from Algorithm1 in the paper is reversed):
                if (t > T_n_prime) 
                {
                    // during the finishing stage (n== n_star && t > T_n_prime), draw a standard RANSAC sample
                    // The sample contains m points selected from U_n at random
                    deal(n, m, sample);
                }
                else
                {
                    // The sample contains m-1 points selected from U_{n−1} at random and u_n
                    deal(n - 1, m - 1, sample);
                    sample[m - 1] = n - 1;
                }
                
           
                ModelType model;
                bool ok = m_driver.fit_model(sample, model);
                if(!ok)
                    continue;
                                
                   
                float inliersErrors = 0;
                 I_N = m_driver.count_inliers(model, inliersFlags, inliersErrors);
                
                    
                    if (I_N > I_N_best) 
                    {
                        int n_best; // best value found so far in terms of inliers ratio
                        int I_n_best; // number of inliers for n_best
                        int I_N_draw; // number of inliers withing the N_draw first data
                        
                        
                        bestModel = model;
                        I_N_best = I_N;
                                            
                        
                        // Select new termination length n_star if possible, according to Sec. 2.2.
                        // Note: the original paper seems to do it each time a new sample is drawn,
                        // but this really makes sense only if the new sample is better than the previous ones.
                        n_best = N;
                        I_n_best = I_N_best;
#ifndef PROSAC_DISABLE_N_STAR_OPTIMIZATION
                        I_N_draw = std::accumulate(isInlier.begin(), isInlier.begin() + N_draw, 0);

                        int n_test; // test value for the termination length
                        int I_n_test; // number of inliers for that test value
                        double epsilon_n_best = (double)I_n_best/n_best;
                        
                        for (n_test = N, I_n_test = I_N_draw; n_test > m; n_test--) {
                            // Loop invariants:
                            // - I_n_test is the number of inliers for the n_test first correspondences
                            // - n_best is the value between n_test+1 and N that maximizes the ratio I_n_best/n_best
                            assert(n_test >= I_n_test);
                            
                            // * Non-randomness : In >= Imin(n*) (eq. (9))
                            // * Maximality: the number of samples that were drawn so far must be enough
                            // so that the probability of having missed a set of inliers is below eta=0.01.
                            // This is the classical RANSAC termination criterion (HZ 4.7.1.2, eq. (4.18)),
                            // except that it takes into account only the n first samples (not the total number of samples).
                            // kn_star = log(eta0)/log(1-(In_star/n_star)^m) (eq. (12))
                            // We have to minimize kn_star, e.g. maximize I_n_star/n_star
                            //printf("n_best=%d, I_n_best=%d, n_test=%d, I_n_test=%d\n",
                            //        n_best,    I_n_best,    n_test,    I_n_test);
                            // a straightforward implementation would use the following test:
                            //if (I_n_test > epsilon_n_best*n_test) {
                            // However, since In is binomial, and in the case of evenly distributed inliers,
                            // a better test would be to reduce n_star only if there's a significant improvement in
                            // epsilon. Thus we use a Chi-squared test (P=0.10), together with the normal approximation
                            // to the binomial (mu = epsilon_n_star*n_test, sigma=sqrt(n_test*epsilon_n_star*(1-epsilon_n_star)).
                            // There is a significant difference between the two tests (e.g. with the findSupport
                            // functions provided above).
                            // We do the cheap test first, and the expensive test only if the cheap one passes.
                            if (( I_n_test * n_best > I_n_best * n_test ) &&
                                ( I_n_test > epsilon_n_best * n_test + std::sqrt(n_test * epsilon_n_best * (1. - epsilon_n_best) * 2.706) )) {
                                if (I_n_test < Prosac_Imin(m,n_test,beta)) {
                                    // equation 9 not satisfied: no need to test for smaller n_test values anyway
                                    break; // jump out of the for(n_test) loop
                                }
                                n_best = n_test;
                                I_n_best = I_n_test;
                                epsilon_n_best = (double)I_n_best / n_best;
                            }
                            
                            // prepare for next loop iteration
                            I_n_test -= isInlier[n_test - 1];
                        } // for(n_test ...
#endif // #ifndef PROSAC_DISABLE_N_STAR_OPTIMIZATION
                        
                        // is the best one we found even better than n_star?
                        if ( I_n_best * n_star > I_n_star * n_best ) 
                        {
                            assert(n_best >= I_n_best);
                            // update all values
                            n_star = n_best;
                            I_n_star = I_n_best;
                            k_n_star = niter_RANSAC(1. - eta, 1. - I_n_star / (double)n_star, m, T_N);
                        }
                    } // if (I_N > I_N_best)
                
            } // while(t <= k_n_star ...
          
            float inliersError;
            int maxCount = m_driver.count_inliers( bestModel, mask, inliersError);
            return maxCount;
        }

    };
}

