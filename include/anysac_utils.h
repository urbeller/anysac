#pragma once

#include <iostream>
#include <vector>
#include <random>

namespace AnySac
{

  template <typename T>
    T clip(const T& n, const T& lower, const T& upper) 
    {
      return std::max(lower, std::min(n, upper));
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
