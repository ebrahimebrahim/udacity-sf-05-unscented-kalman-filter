#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  n_x_=5;
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 6;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  n_aug_ = n_x_+2;
  lambda_ = 3-n_x_;
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(1/(2*(lambda_+n_aug_)));
  weights_(0)=lambda_/(lambda_+n_aug_);
  is_initialized_ = false;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_){
    if (meas_package.sensor_type_==MeasurementPackage::LASER){
      const double px = meas_package.raw_measurements_(0);
      const double py = meas_package.raw_measurements_(1);
      x_ << px, py, 0, 0, 0; // Assume 0 for initial v, psi, psidot, with very high uncertainty
      P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
            0, std_laspy_ * std_laspy_, 0, 0, 0,
            0, 0, 50, 0, 0,
            0, 0, 0, 4*M_PI, 0,
            0, 0, 0, 0, M_PI;
    }
    else if (meas_package.sensor_type_==MeasurementPackage::RADAR) {
      const double rho = meas_package.raw_measurements_(0);
      const double phi = meas_package.raw_measurements_(1);
      const double rhodot = meas_package.raw_measurements_(2);
      const double px = rho*cos(phi);
      const double py = rho*sin(phi);
      x_ << px, py,
            rhodot, // Assume that there's zero tangential velocity initially, with high uncertainty
            0, 0; // Assume 0 for initial psi, psidot, with very high uncertainty

      // Using some taylor expansions we can translate radar position uncertainties into px,py terms
      double std_px = std_radr_ * abs(cos(phi)) + abs(py)*std_radphi_ + abs(sin(phi))*std_radr_*std_radphi_;
      double std_py = std_radr_ * abs(sin(phi)) + abs(px)*std_radphi_ + abs(cos(phi))*std_radr_*std_radphi_;
      
      // There might be correlation between px and py belief, but let's pretend there isn't initially.

      P_ << std_px * std_px, 0, 0, 0, 0,
            0, std_py * std_py, 0, 0, 0,
            0, 0, 50, 0, 0,
            0, 0, 0, 4*M_PI, 0,
            0, 0, 0, 0, M_PI;
    }
    is_initialized_ = true;
    return;
  }

  Prediction(meas_package.timestamp_ - time_us_);

  if (meas_package.sensor_type_==MeasurementPackage::LASER)
    UpdateLidar(meas_package);
  else if (meas_package.sensor_type_==MeasurementPackage::RADAR)
    UpdateRadar(meas_package);
  
  time_us_ = meas_package.timestamp_;

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}