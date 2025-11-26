#include "prosac.h"
#include "ransac.h"

#include <iostream>
#include <vector>
#include <random>

using namespace std;

#define TOLERANCE 1e-10
#define POW2(X) ((X) * (X))



template<class T>
inline bool isZero(T x)
{
  return (x < TOLERANCE) && (x > -TOLERANCE);
}

struct Vec2f
  {
    Vec2f() : x(0), y(0) {}
    Vec2f(float _x, float _y) : x(_x), y(_y){}

    inline float normSq() const
    {
      return x * x + y * y;
    }

    inline float distSq(const Vec2f &p) const
    {
      float dx = x - p.x;
      float dy = y - p.y;

      return dx * dx + dy * dy;
    }

    Vec2f& operator=(const Vec2f& v)
    {
      x = v.x;
      y = v.y;
      return *this;
    }

    Vec2f operator+(const Vec2f& v) const
    {
      return Vec2f(x + v.x, y + v.y);
    }

    Vec2f operator-(const Vec2f& v) const
    {
      return Vec2f(x - v.x, y - v.y);
    }

	  Vec2f operator*(const Vec2f& v) const
	  {
		  return Vec2f(x * v.x, y * v.y);
	  }
	  
    template<class T>
      Vec2f operator*(T s) const
      {
        return Vec2f(s * x, s * y);
      }

    template<class T>
      friend Vec2f operator*(T s, const Vec2f& v);

    friend std::ostream& operator<<(std::ostream& os, const Vec2f& v);

    float dot(const Vec2f& v) const
    {
      return x * v.x + y * v.y;
    }

    float x;
    float y;
  };

  struct Mat33f
  {
    Mat33f(){}
    Mat33f(float mm[3][3])
    {
       for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
          m[i][j] = mm[i][j];
    }
    
    Mat33f(const Mat33f& mm)
    {
       for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
          m[i][j] = mm.m[i][j];
    }
    
    static Mat33f eye(float v = 1.0f)
    {
      Mat33f m;
      
      for(int i = 0; i < 3; ++i)
          m.m[i][i] = v;
          
      return m;
    }
      

     Mat33f operator*(const Mat33f& other) const
    {
      Mat33f tmp;
      for(int i = 0; i < 3; ++i)
      {
          for(int j = 0; j < 3; ++j)
          {
              tmp.m[i][j] = 0;
              for(int k = 0; k < 3; ++k)
                  tmp.m[i][j] += this->m[i][k] * other.m[k][j]; 
          }
      }
      
      return tmp;
    }
      
    friend std::ostream& operator<<(std::ostream& os, const Mat33f& m);

    float m[3][3] = {{0}};
  };


void normalizingValues(const vector<Vec2f> &pts, Vec2f &center, float &scale)
{
  center.x = center.y = 0.;

  for(auto p:pts)
  {
    center = center + p;  
  }

  center.x /= pts.size();
  center.y /= pts.size();

  scale = 0;
  for(auto p:pts)
  {
    scale += sqrt(POW2(p.x - center.x) + POW2(p.y - center.y));
  }

  scale /= pts.size();
  
  scale = sqrt(2) / scale;
}


void normalizePts(const vector<Vec2f> &src, vector<Vec2f> &srcn, Vec2f &shift, float &scale)
{
  int size = src.size();
  srcn.resize(size);
  normalizingValues(src, shift, scale);
  for(int i = 0; i < size; ++i)
  {
    srcn[i].x = scale * (src[i].x - shift.x);
    srcn[i].y = scale * (src[i].y - shift.y);
  }

}


bool findExactHomography(const std::vector<Vec2f> &_src, const std::vector<Vec2f> &_dst, Mat33f &H, bool normalize = true)
{

  constexpr int size = 4;
  std::vector<Vec2f> src(size), dst(size);
  Vec2f shift1, shift2;
  float scale1, scale2;

  if(normalize)
  {
    normalizePts(_src, src ,shift1, scale1);
    normalizePts(_dst, dst ,shift2, scale2);
  }
  else
  {
    src = _src;
    dst = _dst;
  }

    double Kc[3][3];
    Kc[0][0] = dst[0].x; Kc[1][0] = dst[0].y; Kc[2][0] = 1.0;
    Kc[0][1] = dst[1].x; Kc[1][1] = dst[1].y; Kc[2][1] = 1.0;
    Kc[0][2] = dst[2].x; Kc[1][2] = dst[2].y; Kc[2][2] = 1.0;
    
    
    double adjKc[3][3];
    adjKc[0][0] = dst[1].y - dst[2].y;  adjKc[0][1] = dst[2].x - dst[1].x;  adjKc[0][2] = dst[1].x * dst[2].y - dst[2].x * dst[1].y ;
    adjKc[1][0] = dst[2].y - dst[0].y;  adjKc[1][1] = dst[0].x - dst[2].x;  adjKc[1][2] = dst[2].x * dst[0].y - dst[0].x * dst[2].y ;
    adjKc[2][0] = dst[0].y - dst[1].y;  adjKc[2][1] = dst[1].x - dst[0].x;  adjKc[2][2] = dst[0].x * dst[1].y - dst[1].x * dst[0].y ;
    
    double adjKp[3][3];
    adjKp[0][0] = src[1].y - src[2].y;  adjKp[0][1] = src[2].x - src[1].x;  adjKp[0][2] = src[1].x * src[2].y - src[2].x * src[1].y ;
    adjKp[1][0] = src[2].y - src[0].y;  adjKp[1][1] = src[0].x - src[2].x;  adjKp[1][2] = src[2].x * src[0].y - src[0].x * src[2].y ;
    adjKp[2][0] = src[0].y - src[1].y;  adjKp[2][1] = src[1].x - src[0].x;  adjKp[2][2] = src[0].x * src[1].y - src[1].x * src[0].y ;
    
    
    double diagP[3];
    diagP[0] = adjKp[0][0] * src[3].x + adjKp[0][1] * src[3].y + adjKp[0][2] ;
    diagP[1] = adjKp[1][0] * src[3].x + adjKp[1][1] * src[3].y + adjKp[1][2] ;
    diagP[2] = adjKp[2][0] * src[3].x + adjKp[2][1] * src[3].y + adjKp[2][2] ;

    if( isZero(diagP[0]) || isZero(diagP[1]) || isZero(diagP[2]))
        return false;
    
    double diagC[3];
    diagC[0] = adjKc[0][0] * dst[3].x + adjKc[0][1] * dst[3].y + adjKc[0][2] ;
    diagC[1] = adjKc[1][0] * dst[3].x + adjKc[1][1] * dst[3].y + adjKc[1][2] ;
    diagC[2] = adjKc[2][0] * dst[3].x + adjKc[2][1] * dst[3].y + adjKc[2][2] ;
    
    
    diagC[0] /= diagP[0];
    diagC[1] /= diagP[1];
    diagC[2] /= diagP[2];
    
    for(int i = 0; i < 3; ++i)
    {
        Kc[i][0] *= diagC[0];
        Kc[i][1] *= diagC[1];
        Kc[i][2] *= diagC[2];
    }


    // Compute H[2][2] for normalization.
    double v_22 = Kc[2][0] * adjKp[0][2] + Kc[2][1] * adjKp[1][2] + Kc[2][2] * adjKp[2][2] ;
    if(isZero(v_22))
      return false;



    v_22 = 1.0 / v_22; 

    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
        {
            double v = 0;
            for(int k = 0; k < 3; ++k)
                v += Kc[i][k] * adjKp[k][j];

            H.m[i][j] = float(v * v_22);
        }
    

    if(normalize)
    {
      Mat33f tm, tsol;
      tm = Mat33f::eye();
      tm.m[0][0] = scale1; tm.m[0][2] = -shift1.x * scale1;
      tm.m[1][1] = scale1; tm.m[1][2] = -shift1.y * scale1;
      tsol = H * tm;

      tm.m[0][0] = 1 / scale2; tm.m[0][2] = shift2.x;
      tm.m[1][1] = 1 / scale2; tm.m[1][2] = shift2.y;
      H = tm * tsol;
    }

    return true;
}

class HomogDriver
{

public:
  HomogDriver(const vector<Vec2f> &src, const vector<Vec2f> &dst, float tol, bool normalize) : m_src(src), m_dst(dst), m_tol(tol), m_normalize(normalize)
  {}

  void setNormalization(bool nrm)
  {
    m_normalize = nrm;
  }

  // Required methods.
  int model_size() const
  {
    return m_modelSize;
  }

  int data_size() const
  {
    return m_src.size();
  }

  bool fit_model(const vector<int> &indices, Mat33f &model)
  {

    std::vector<Vec2f> sp(4);
    std::vector<Vec2f> dp(4);

    for(int i = 0; i < 4; ++i)
    {
      sp[i] = m_src[indices[i]];
      dp[i] = m_dst[indices[i]];
    }

    bool status = findExactHomography(sp, dp, model, m_normalize);
    return status;
  }

  int count_inliers(const Mat33f &model, std::vector<unsigned char> &flags)
  {
    const float tol2 =  m_tol * m_tol;
    int ngood = 0;
  
    int sz = m_src.size();
    flags.resize(sz);

    for(int i = 0; i < sz; ++i)
    {
      float x = m_src[i].x;
      float y = m_src[i].y;
      float xp = model.m[0][0] * x + model.m[0][1] * y + model.m[0][2];
      float yp = model.m[1][0] * x + model.m[1][1] * y + model.m[1][2];
      float zp = model.m[2][0] * x + model.m[2][1] * y + model.m[2][2];
      
      if(!isZero(zp))
      {
        xp /= zp;
        yp /= zp;

        float err = POW2(xp - m_dst[i].x) + POW2(yp - m_dst[i].y);
        flags[i] = (err < tol2);
        ngood += flags[i];
      }
    }

    return ngood;
  }

  bool check_subset(const vector<int> &indices)
  {
    return true;
  }


  bool get_subset( vector<int> &indices )
  {
    return false;
  }

  bool has_generator() const
  {
    return false; // set this to false to use Ransac<> generator.
  }

private:
  static const int m_modelSize = 4;
  const vector<Vec2f> &m_src, &m_dst;
  float m_tol = 2.0f;
  bool m_normalize = true;
};

std::ostream& operator<<(std::ostream& os, const Vec2f& v)
{
  os << "(" << v.x << "," << v.y << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const Mat33f& m)
  {
    for(int i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 3; ++j)
        os << m.m[i][j] << " ";

      os << std::endl;
    }

    return os;
  }


int main(int argc, char **argv)
{

  std::random_device rd;  
  std::mt19937 gen(rd()); 
  std::uniform_real_distribution<> dis(-100.0, 100.0);

  int size = 1000;
  std::vector<Vec2f> src(size), dst(size);
  float mm[3][3] = {{10,0,100},{0,10,200},{0,0,1}};
  Mat33f realH(mm);
  AnySac::GaussianNoiser<float> gn(0, 0.5);

  for(int i = 0; i < size; ++i)
  {
    float x = drand48() * 1000;
    float y = drand48() * 500;
    src[i].x = x;
    src[i].y = y;

    float xp = realH.m[0][0] * x + realH.m[0][1] * y + realH.m[0][2];
    float yp = realH.m[1][0] * x + realH.m[1][1] * y + realH.m[1][2];
    float zp = realH.m[2][0] * x + realH.m[2][1] * y + realH.m[2][2];

    dst[i].x = xp / zp + gn.sample();
    dst[i].y = yp / zp + gn.sample();

    if(i > int(size * 0.8))
    {
      dst[i].x = drand48() * 1000;
      dst[i].y = drand48() * 500;
    }


  }


  std::vector<Vec2f> srcn(size), dstn(size);
  Vec2f shift1, shift2;
  float scale1, scale2;
  normalizePts(src, srcn ,shift1, scale1);
  normalizePts(dst, dstn ,shift2, scale2);


  // Init the driver.
  float tol = 10.0f;
  
  // Solve the problem.
  Mat33f sol, tsol, tm;
  HomogDriver hd(src, dst, tol, true); 
  AnySac::Prosac<HomogDriver, Mat33f> psac(hd);
  AnySac::Ransac<HomogDriver, Mat33f> rsac(hd);
  psac.solve(sol);
  cout << sol << std::endl;

  rsac.solve(sol);
  cout << sol << std::endl;


  return 0;
}



