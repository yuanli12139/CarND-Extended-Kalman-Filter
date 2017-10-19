#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
  	cout << "Invalid estimation or ground_truth data" << endl;
  	return rmse;
  }

  //accumulate squared residuals
  VectorXd residuals(4);
  residuals << 0, 0, 0, 0;
  for (int i = 0; i < estimations.size(); ++i) {
  	VectorXd residual = (estimations[i] - ground_truth[i]).array() * (estimations[i] - ground_truth[i]).array();
  	residuals += residual;
  }

  //calculate the mean
  VectorXd mean = residuals / estimations.size();

  //calculate the squared root
  rmse = mean.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);

  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float d = sqrt(px * px + py * py);
  float d_2 = d * d;
  float d_3 = d_2 * d;

  //check division by zero
  if (d) {
  	//compute the Jacobian matrix
  	Hj << px / d,                 py / d,                 0,      0, 
          -py / d_2,              px / d_2,               0,      0,
          py*(vx*py-vy*px) / d_3, px*(vy*px-vx*py) / d_3, px / d, py / d;
  } 
  else 
  	cout << "CalculateJacobian() - Error - Division by Zero" << endl;

  return Hj;
}
