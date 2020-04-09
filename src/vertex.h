#ifndef vertex_H
#define vertex_H

// #include "utility/vector.h"
#include "global.h"
#include <array>
#include <vector>

double sum2(const momentum &);
double norm2(const momentum &);

namespace ver {

  enum channel { I = 0, T, U, S };

  const int MAXSIGMABIN = 100000;

  class weightMatrix {
    // 2x2 matrix of weight; Direct/Exchange and chain/Lver/Rver/other
  public:
    weightMatrix() { SetZero(); }
    void SetZero() {
      for (auto &i : _Weight)
        i = 0.0;
    }
    double Sum() {
      double sum = 0;
      for (auto &i : _Weight)
        sum += i;
      return sum;
    }

    double Abs() {
      double sum = 0;
      // for (auto &i : _Weight)
      //   sum += fabs(i);
      // sum = fabs(_Weight[DIR] + _Weight[EX] / 2.0) + fabs(_Weight[EX] / 2.0);
      sum = fabs(_Weight[DIR] + _Weight[EX] / 2.0);
      return sum;
    }
    double &operator()(int dir) { return _Weight[dir]; }

  private:
    array<double, 2> _Weight;
  };

  class fermi {
  public:
    fermi();
    double Green(double Tau, const momentum &Momentum, spin Spin, int GType,
                 double Scale = 0);

  private:
    // beyond which the expense sigma function will be called
    double UpperBound, LowerBound;
    double DeltaK;
    double UpperBound2, LowerBound2; // lower upbound for better sigma
    double DeltaK2;
    double PhyGreen(double Tau, const momentum &Mom, int GType, double Scale = 0);
    double FockSigma(const momentum &Mom);
    double BuildFockSigma();
    double Fock(double k);
    // warning: this function only works for T=0 and 3D!!!!
    double GetSigma(double k);
    double Sigma[MAXSIGMABIN];
    double Sigma2[MAXSIGMABIN];
  };

  class verTensor {
  public:
    verTensor();
    ~verTensor();
    void Initialize();
    double &Interaction(int Angle, int ExtQ, int Dir);
    double &Estimator(int Order, int Angle, int ExtQ, int Dir);

  private:
    double *_Estimator;
    double *_Interaction;
    int QIndex;
    int AngleIndex;
    int OrderIndex;
  };

  class delta{
    // load F from file, measure \Delta, and save \Delta to file
  public:
    delta();
    ~delta();
    void Initialize();

    void LoadF();
    void Save(bool Simple);

    void Measure(const momentum &InL, const momentum &Tran,
                 const int QIndex, int Freq,int Order,
                 ver::weightMatrix &Weight, double Factor);

    double &DeltaVal(int Freq, int Angle, int ExtQ, int Dir);
    double &FVal(int tau, int ExtQ);
    double F(double tau, int ExtQ);
  private:
    double Normalization;
    double PhyWeightT;
    int AngleIndex;
    int FreqIndex;
    vector<double> _data;
    vector<double> _f;
    vector<double> _taulist;
  };

  class verQTheta {
  public:
    verQTheta();
    void Interaction(const array<momentum *, 4> &LegK, double Tau, bool IsRenorm,
                     int CTOrder, double &WeightDir, double &WeightEx);

    void Measure(const momentum &InL, const momentum &InR, const int QIndex,
                 int Order, double Tau, int Channel, ver::weightMatrix &Weight,
                 double Factor);
    void Update(double Ratio, int Order);
    void Save(bool Simple = false);
    void ClearStatis();
    void LoadWeight();
    void ResetIRScale(int IRScaleBin);

    array<verTensor, 4> Chan;
    delta Delta;
    // double TauBasis[TauBinSize][TauBasisNum];

    double Normalization;
    double PhyWeightT;
    double PhyWeightI;
  };

  double Angle3D(const momentum &K1, const momentum &K2);
  double Index2Angle(const int &Index, const int &AngleNum);
  int Angle2Index(const double &Angle, const int &AngleNum);
  void _TestAngleIndex();
  void _TestAngle2D();

  double Index2Mom(const int &Index);
  int Mom2Index(const double &K);

  double Index2Scale(const int &Index);
  int Scale2Index(const double &Scale);

  double Index2Tau(const int &Index);
  int Tau2Index(const double &Tau);

}; // namespace ver
#endif
