#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

//The constants we're refering
//placed in header so that i can access them from main
//as well as from MPC.cpp
const size_t c_N = 10;
const double c_dt = 0.1;
const double c_Lf = 2.67;
const double c_ref_v = 100.;
const double c_delay = 0.1;



class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
