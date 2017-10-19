#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &Hj_in, MatrixXd &R_laser_in, MatrixXd &R_radar_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  Hj_ = Hj_in;
  R_laser_ = R_laser_in;
  R_radar_ = R_radar_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  // For this project, however, we do not need to use the f function or Fj.
  // If we had been using a non-linear model in the prediction step, we would need to replace the F matrix with its Jacobian, Fj. 
  // However, we are using a linear model for the prediction step. So, for the prediction step, we can still use the regular Kalman filter equations and the F matrix rather than the extended Kalman filter equations.
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_laser_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  // Note: The equation y = z − Hx' for the Kalman filter does not become y = z − Hjx for the extended Kalman filter. 
  // Instead, for extended Kalman filters, we'll use the h function directly to map predicted locations x′ from Cartesian to polar coordinates.
  VectorXd hx = VectorXd(3);
  float px = x_(0), py = x_(1), vx = x_(2), vy = x_(3);

  float rho = sqrt(px * px + py * py);
  float phi = atan2(py, px);
  float rho_dot = (px * vx + py * vy) / rho;
  hx << rho, phi, rho_dot;

  VectorXd y = z - hx;
  while (y(1) > M_PI)
    y(1) -= 2 * M_PI;
  while (y(1) < -M_PI)
    y(1) += 2 * M_PI;

  Hj_ = tools.CalculateJacobian(x_);
  MatrixXd S = Hj_ * P_ * Hj_.transpose() + R_radar_;
  MatrixXd K = P_ * Hj_.transpose() * S.inverse();

  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * Hj_) * P_;
}

