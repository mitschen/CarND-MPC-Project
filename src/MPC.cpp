#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t const N = 25; //following suggestion of 19/9
double const dt = 0.5; //following suggestion of 19/9

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// Both the reference cross track and orientation errors are 0.
// The reference velocity is set to 40 mph.
double ref_v = 40;


const int idx_x_start(0);
const int idx_y_start(N);
const int idx_psi_start(2*N);
const int idx_v_start(3*N);
const int idx_cte_start(4*N);
const int idx_epsi_start(5*N);
const int idx_delta_start(6*N);
const int idx_acc_start(7*N-1);


class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // Cost function
    // TODO: Define the cost related the reference state and
    // any anything you think may be beneficial.

    // The part of the cost based on the reference state.
    AD<double> const ref_v(35.);
    for (int t(0); t < N; t++) {
      fg[0] += CppAD::pow(vars[idx_cte_start + t], 2);
      fg[0] += CppAD::pow(vars[idx_epsi_start + t], 2);
      fg[0] += CppAD::pow(vars[idx_y_start + t] - ref_v, 2);
    }

    // Minimize the use of actuators.
    for (int t = 0; t < N - 1; t++) {
      fg[0] += CppAD::pow(vars[idx_delta_start + t], 2);
      fg[0] += CppAD::pow(vars[idx_acc_start + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int t = 0; t < N - 2; t++) {
      fg[0] += CppAD::pow(vars[idx_delta_start + t + 1] - vars[idx_delta_start + t], 2);
      fg[0] += CppAD::pow(vars[idx_acc_start + t + 1] - vars[idx_acc_start + t], 2);
    }



    //Initialize the init-state
    fg[1 + idx_x_start] = vars[idx_x_start];
    fg[1 + idx_y_start] = vars[idx_y_start];
    fg[1 + idx_psi_start] = vars[idx_psi_start];
    fg[1 + idx_v_start] = vars[idx_v_start];
    fg[1 + idx_cte_start] = vars[idx_cte_start];
    fg[1 + idx_epsi_start] = vars[idx_epsi_start];


    for (int t(1); t < N ; t++) {
      // psi, v, delta at time t
      AD<double> psi0 = vars[idx_psi_start + t - 1];
      AD<double> v0 = vars[idx_v_start + t - 1];
      AD<double> delta0 = vars[idx_delta_start + t - 1];

      // psi at time t+1
      AD<double> psi1 = vars[idx_psi_start + t];

      // how psi changes
      fg[1 + idx_psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
    }

    for (int t(1); t < N; t++) {
      // The state at time t+1 .
      AD<double> x1 = vars[idx_x_start + t];
      AD<double> y1 = vars[idx_y_start + t];
      AD<double> psi1 = vars[idx_psi_start + t];
      AD<double> v1 = vars[idx_v_start + t];
      AD<double> cte1 = vars[idx_cte_start + t];
      AD<double> epsi1 = vars[idx_epsi_start + t];

      // The state at time t.
      AD<double> x0 = vars[idx_x_start + t - 1];
      AD<double> y0 = vars[idx_y_start + t - 1];
      AD<double> psi0 = vars[idx_psi_start + t - 1];
      AD<double> v0 = vars[idx_v_start + t - 1];
      AD<double> cte0 = vars[idx_cte_start + t - 1];
      AD<double> epsi0 = vars[idx_epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[idx_delta_start + t - 1];
      AD<double> a0 = vars[idx_acc_start + t - 1];

      AD<double> f0 = coeffs[0] + coeffs[1] * x0;
      AD<double> psides0 = CppAD::atan(coeffs[1]);

      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // Recall the equations for the model:
      // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
      // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
      // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
      // v_[t] = v[t-1] + a[t-1] * dt
      // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
      // epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
      fg[1 + idx_x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + idx_y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + idx_psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + idx_v_start + t] = v1 - (v0 + a0 * dt);
      fg[1 + idx_cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + idx_epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }

  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
//  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  int const c_no_states(6); //x,y, psi, v, cte, epsi
  int const c_no_actuators(2); //delta, accelaration
  size_t const n_vars(c_no_states * N + c_no_actuators * (N-1));
  // TODO: Set the number of constraints
  size_t n_constraints = 6*N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  //set the initial state vector
  vars[idx_x_start]     = state[0]; //set x
  vars[idx_y_start]     = state[1]; //set y
  vars[idx_psi_start]   = state[2]; //set psi
  vars[idx_v_start]     = state[3]; //set v
  vars[idx_cte_start]   = state[4]; //set cte
  vars[idx_epsi_start]  = state[5]; //set epsi



  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
  //set limits for delta which is 25°
  double const c_pos_noLimit(std::numeric_limits<double>::max());
  double const c_neg_noLimit(std::numeric_limits<double>::min());
  for(int i(0), _maxI(6*N); i<_maxI; i++)
  {
    vars_lowerbound[i] = c_neg_noLimit;
    vars_upperbound[i] = c_pos_noLimit;
  }
  double const rad25Deg(0.436332); //25/(M_PI/180)
  for(int i(idx_delta_start), _maxI(7*N); i<_maxI; i++)
  {
    vars_lowerbound[i] = -rad25Deg;
    vars_upperbound[i] = rad25Deg;
  }
  //set constraints for acceleration which is -1 <-> 1
  for(int i(idx_acc_start), _maxI(8*N-1); i < _maxI; i++)
  {
    vars_lowerbound[i] = -1.;
    vars_upperbound[i] = 1.;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = vars[i];
    constraints_upperbound[i] = vars[i];
  }

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.


  //we'll provide the steering angle and the velocity and furthermore
  //the x, y coordinates of our trajectory
  vector<double> result(N+2);
  result.push_back(solution.x[idx_delta_start]);
  result.push_back(solution.x[idx_acc_start]);
  for(int i(0); i<N; i++)
  {
    result.push_back(solution.x[idx_x_start+i*N]);
    result.push_back(solution.x[idx_y_start+i*N]);
  }

  return result;
}
