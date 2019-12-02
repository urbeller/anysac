#include "ransac.h"

#include <iostream>
#include <vector>
#include <random>

using namespace std;


// Struct for 2d points.
struct Pt2d
{
	Pt2d(double x_, double y_) : x(x_), y(y_){}
	Pt2d() :x(0.), y(0.){}

	double x, y;

	double dist(const Pt2d &other) const { return sqrt( pow(x - other.x,2) + pow(y - other.y,2));}
};

class LineDriver
{

public:
	LineDriver(const vector<Pt2d> &pts_, double th_) : pts(pts_), threshold(th_),
			rnd_gen((unsigned int)time(0)), rnd_dist(0, pts.size() - 1)
	{}

	// Required methods.
	int model_size() const;
	int data_size() const;
	bool fit_model(const vector<int> &indices, vector<Pt2d> &models);
	int count_inliers(const Pt2d &model);
	bool check_subset(const vector<int> &indices);
	bool get_subset( vector<int> &pickedIndex);
	bool has_generator() const
	{
		return true; // set this to false to use Ransac<> generator.
	}

private:
	const vector<Pt2d> &pts;
	double threshold;
	const int modelSize = 2;

	// To show an example of a user-defined indices-generator.
	std::mt19937 rnd_gen;
	std::uniform_int_distribution<> rnd_dist;
};

int main(int argc, char **argv)
{

	std::random_device rd;  
	std::mt19937 gen(rd()); 
	std::uniform_real_distribution<> dis(-100.0, 100.0);

	// Generate noisy points on a line.
	vector<Pt2d> pts(100);
	double slope = 10, b = 5;

	for (int i = 0; i < pts.size(); ++i)
	{
		pts[i].x = dis(gen);

		double noise = (dis(gen) / 200.0) ;
		pts[i].y = pts[i].x * slope + b +  noise;
	}


	// Init the driver.
	LineDriver ld(pts, 0.01);
	
	// Solve the problem.
	Pt2d sol;
	AnySac::Ransac<LineDriver, Pt2d> sac(ld);
	sac.solve(sol);

	cout << sol.x << " " << sol.y << endl;
}



	int LineDriver::model_size() const 
	{ 
		return modelSize;
	}

	int LineDriver::data_size() const 
	{ 
		return pts.size();
	}
	
	bool LineDriver::fit_model(const vector<int> &indices, vector<Pt2d> &models)
	{
		const Pt2d &p0 = pts[indices[0]];
		const Pt2d &p1 = pts[indices[1]];

		if (p0.dist(p1) < threshold)
			return false;

		models.resize(1);

		double slope = (p1.y - p0.y) / (p1.x - p0.x);
		double b = p0.y - slope * p0.x;

		models[0] = Pt2d(slope, b);

		return true;
	}


	int LineDriver::count_inliers(const Pt2d &model)
	{
		int inliers = 0;

		for (int i = 0; i < pts.size(); ++i)
		{
			const Pt2d &p = pts[i];

			double dist = fabs(p.x * model.x + model.y - p.y);
			inliers += (dist < threshold);
		}

		return inliers;
	}

	bool LineDriver::check_subset(const vector<int> &indices)
	{
		// Consider only first-n indices.
		if(indices.size() < modelSize || indices[0] == indices[1])
			return false;

		const Pt2d &p0 = pts[indices[0]];
		const Pt2d &p1 = pts[indices[1]];

		if (p0.dist(p1) < threshold)
			return false;
		
		return true;
	}



bool LineDriver::get_subset( vector<int> &indices)
{
	indices.resize(2, 0);
	int max_iters = 10;

	int niter = 0;
	while (indices[0] == indices[1] && niter < max_iters)
	{
		indices[0] = rnd_dist(rnd_gen);
		indices[1] = rnd_dist(rnd_gen);
		++niter;
	}

	return niter < max_iters;
}
