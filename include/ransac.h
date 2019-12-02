#ifndef RANSAC_H_
#define RANSAC_H_

#include <iostream>
#include <vector>
#include <random>
#include <limits>

using namespace std;

namespace AnySac
{
	// A driver must define the following methods:
	//
	// int model_size() const;
	// int data_size() const;
	// bool fit_model(const vector<int> &indices, vector<ModelType> &models);
	// int count_inliers(const ModelType &model);
	// bool check_subset(const vector<int> &indices);
	// bool has_generator();
	// bool get_subset(vector<int> &indices);

	template< typename DriverType, typename ModelType>
		class Ransac
		{
			public:

				Ransac( DriverType &k , double conf_=0.99) : 
					driver(k), 
					maxAttempts(300), 
					maxIters(300), 
					conf(conf_),
					rnd_gen((unsigned int)time(0)),
					rnd_dist(0, k.data_size() - 1)

			{
			}

				// Solve the problem.
				// Input: a model object (driver).
				// Returns: number of inliers.
				int solve(ModelType &model)
				{

					int         modelSize = driver.model_size();
					int         dataSize = driver.data_size();
					int         niter    = maxIters;
					int         maxCount = 0;
					vector<int> ndices(modelSize);
					vector<ModelType> models;
					bool has_generator = driver.has_generator();

					for (int iter=0; iter<niter; iter++)
					{
						if ( dataSize >= modelSize )
						{


							bool found = (has_generator) ? driver.get_subset(ndices) : get_subset(ndices);

							if( !found )
							{
								if( iter == 0 )
									return 0;

								break;
							}

							// We got a valid subset. Estimate a model.
							if( !driver.fit_model ( ndices, models ) )
							{
								//return 0;
							}

							// Find the best model.
							for(unsigned int m = 0, len = models.size(); m<len; m++)
							{
								int count = driver.count_inliers( models[m] );

								if( count > std::max(maxCount, modelSize-1) )
								{
									maxCount = count;

									model = models[m];
									niter = update_niters(conf, double(dataSize-count)/dataSize, modelSize, niter);
								}
							}

						} 
						else 
						{
							cerr<<__FUNCTION__<<": dataSize is smaller than model size.\n";
							break;
						}

					}

					// Update the inliers count and the inliers mask.
					maxCount = driver.count_inliers( model );
					return maxCount;
				}

			private:

				// Pick a random combination of indices.
				bool get_subset( vector<int> &pickedIndex)
				{

					int j, i=0, iter=0;
					int dataSize  = driver.data_size();
					int modelSize = driver.model_size();

					for (; iter< maxAttempts; iter++)
					{

						for (i=0; i<modelSize && iter<maxAttempts; )
						{
							pickedIndex[i] = rnd_dist(rnd_gen);

							for (j = 0; j < i; j++ )
								if (pickedIndex[j] == pickedIndex[i])
									break;

							if (j < i)
								continue;

							i++;
						}

						if (i==modelSize && !driver.check_subset(pickedIndex))
							continue;

						break;
					}

					return i == modelSize && iter < maxAttempts;
				}


				// Update the minimal number of ransac iterations.
				int update_niters( double p, double ep, int modelPoints, int maxIter) const
				{

					p   = std::max(p, 0.0);
					p   = std::min(p, 1.0);
					ep  = std::max(ep, 0.0);
					ep  = std::min(ep, 1.0);

					// avoid inf's & nan's
					double num = std::max(1.0 - p, std::numeric_limits<double>::min() );
					double denom = 1. - pow(1.0 - ep,modelPoints);
					if( denom < std::numeric_limits<double>::min() )
						return 0;

					num = log(num);
					denom = log(denom);

					return denom >= 0 || -num >= maxIter*(-denom) ?
						maxIter : round(num/denom);
				}


				DriverType &driver;
				int         maxAttempts;
				int         maxIters;
				double      conf;

				// For combinations generation.
				std::mt19937 rnd_gen;
				std::uniform_int_distribution<> rnd_dist;
		};
}

#endif
