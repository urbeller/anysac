#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <cassert>

namespace AnySac
{

  template <typename T>
    T clip(const T& n, const T& lower, const T& upper) 
    {
      return std::max(lower, std::min(n, upper));
    }

    
    inline int niter_RANSAC(double p, // probability that at least one of the random samples picked up by RANSAC is free of outliers
                     double epsilon, // proportion of outliers
                     int s, // sample size
                     int Nmax = -1) // upper bound on the number of iterations (-1 means INT_MAX)
    {
        // compute safely N = ceil(log(1. - p) / log(1. - exp(log(1.-epsilon) * s)))
        if (Nmax == -1) {
            Nmax = std::numeric_limits<int>::max();
        }
        assert(Nmax >= 1);
        if (epsilon <= 0.) {
            return 1;
        }
        // logarg = -(1-epsilon)^s
        double logarg = -std::pow(1.-epsilon, s); /* use -exp(s*log(1.-epsilon) if pow is not avail. */
        // logval = log1p(logarg)) = log(1-(1-epsilon)^s)
        double logval = std::log(1. + logarg); // C++/boost version: logval = boost::math::log1p(logarg)
        double N = std::log(1. - p) / logval;
        if (logval  < 0. && N < Nmax) {
            // for very big N, log1p(x) is more precise than log(1+x), so N values may differ
            assert(N > 1e8 || std::ceil(N) == std::ceil(std::log(1. - p) / std::log(1. - std::exp(std::log(1. - epsilon) * s))));
            return (int)std::ceil(N);
        }
        return Nmax;
    }
    
   inline void choiceKnuth
    (
     int N,    // size of set sampling from
     int n,        // size of each sample
     std::mt19937 &rng,
     std::vector<int> &dst  // output, zero-offset indicies to selected items
    )
    {
      int t = 0; // total input records dealt with
      int m = 0; // number of items selected so far
      double u;
      auto uniform = std::uniform_real_distribution<float>(0, 1);
        
      while (m < n)
      {
        u = uniform(rng); // call a uniform(0,1) random number generator
        if ( (N - t)*u >= n - m )
        {
          t++;
        }
        else
        {
          dst[m] = t;
          t++; m++;
        }
      }
    }

   template<class T>
   struct GaussianNoiser
   {
     GaussianNoiser(T mean, T stddev) : generator(rd()), distribution(mean, stddev) {}

     T sample()
     {
        return distribution(generator);
     }

      private:
         std::random_device rd;
         std::mt19937 generator;
         std::normal_distribution<T> distribution;
   };

}
