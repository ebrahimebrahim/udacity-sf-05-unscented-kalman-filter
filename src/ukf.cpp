#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

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
  lambda_ = 3-n_aug_;

  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(1/(2*(lambda_+n_aug_)));
  weights_(0)=lambda_/(lambda_+n_aug_);
  
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

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
      x_ << px, py, 5, 0, 0; // Assume initial v, psi, psidot, with very high uncertainty
      P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
            0, std_laspy_ * std_laspy_, 0, 0, 0,
            0, 0, 50, 0, 0,
            0, 0, 0, M_PI*M_PI, 0,
            0, 0, 0, 0, 0.25 * M_PI * M_PI;
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
            0, 0, 0, M_PI*M_PI, 0,
            0, 0, 0, 0, 0.25 * M_PI * M_PI;
    }
    is_initialized_ = true;
    return;
  }

  if (meas_package.sensor_type_==MeasurementPackage::LASER &&  !use_laser_)
    return;
  if (meas_package.sensor_type_==MeasurementPackage::RADAR && !use_radar_)
    return;

  Prediction((meas_package.timestamp_ - time_us_)/1e6);

  if (meas_package.sensor_type_==MeasurementPackage::LASER)
    UpdateLidar(meas_package);
  else if (meas_package.sensor_type_==MeasurementPackage::RADAR)
    UpdateRadar(meas_package);
  
  time_us_ = meas_package.timestamp_;

}

// Force an angle into the (-pi, pi] range.
static double standardize_angle(double angle) {
  const double revs = (angle+M_PI)/(2*M_PI);
  const double revs_residual = revs-floor(revs);
  return (revs_residual * 2 * M_PI)-M_PI;
}

// This function f is the process model
static Matrix<double, 5, 1> f(Matrix<double, 7, 1> xsig_aug, double dt){

  const double px=xsig_aug(0);
  const double py=xsig_aug(1);
  const double v=xsig_aug(2);
  const double psi=xsig_aug(3);
  const double psidot=xsig_aug(4);
  const double nu_a=xsig_aug(5);
  const double nu_psidd=xsig_aug(6);
  Matrix<double, 5, 1> step;
  if (abs(psidot)>0.0001){
      step << 
            (v/psidot)*(sin(psi+psidot*dt)-sin(psi)) + 0.5*dt*dt*cos(psi)*nu_a,// px
            (v/psidot)*(-cos(psi+psidot*dt)+cos(psi)) + 0.5*dt*dt*sin(psi)*nu_a,// py
            0+dt*nu_a,// v
            psidot*dt+0.5*dt*dt*nu_psidd,// psi
            0+dt*nu_psidd// psidot
      ;
  }
  else {
      step << 
            v*cos(psi*dt) + 0.5*dt*dt*cos(psi)*nu_a,// px
            v*sin(psi*dt) + 0.5*dt*dt*sin(psi)*nu_a,// py
            0+dt*nu_a,// v
            psidot*dt+0.5*dt*dt*nu_psidd,// psi
            0+dt*nu_psidd// psidot
      ;
  }
  Matrix<double, 5, 1> x_out = xsig_aug.head(5) + step;
  return x_out;
}

// Here delta_t is in seconds. Don't forget to convert!
void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  int n_extra = n_aug_ - n_x_;

  // Augment x and P to get x_aug and P_aug

  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_)= x_;
  x_aug.tail(n_extra) << 0,0 ;

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.block(0,0,n_x_,n_x_) = P_;
  P_aug.block(n_x_,n_x_,n_extra,n_extra) << std_a_*std_a_ , 0.0 ,
                                            0.0 , std_yawdd_*std_yawdd_;

  // Compute sigma points for x_aug and P_aug

  MatrixXd sqrt_of_covariance = P_aug.llt().matrixL();
  MatrixXd pm_mat = sqrt(lambda_+n_aug_) * sqrt_of_covariance; // we shall add and subtract this matrix
  
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1); // Columns of this matrix are the sigma pts
  Xsig_aug.col(0)=x_aug;
  
  for (int i = 0; i<n_aug_; ++i) {
      Xsig_aug.col(i+1) =     x_aug + pm_mat.col(i);
      Xsig_aug.col(n_aug_+i+1) = x_aug - pm_mat.col(i);
  }

  // Send sigma points through process model
  
  for (int i = 0; i<2*n_aug_+1; ++i)
      Xsig_pred_.col(i) = f(Xsig_aug.col(i), delta_t);
  
  // Predict new x and P by combining processed sigma points

  x_.fill(0.0);
  for (int i=0; i<2*n_aug_+1; ++i){
      x_ += weights_(i) * Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  for (int i=0; i<2*n_aug_+1; ++i){
      VectorXd diff = Xsig_pred_.col(i) - x_;
      diff(3) = standardize_angle(diff(3));
      P_ += weights_(i) * diff * diff.transpose();
  }
}

// This function is the lidar measurement model, mapping from internal states into radar meas. space
static Matrix<double, 2, 1> h_lidar(Matrix<double, 5, 1> x){
    const double px=x(0);
    const double py=x(1);
    const double v=x(2);
    const double psi=x(3);
    const double psidot=x(4);
    
    Matrix<double, 2, 1> z;
    z << px,py;
    return z;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  int n_z = 2;

  // This vector consists of the actual measured p_x, p_y
  VectorXd z = meas_package.raw_measurements_;

  // The columns of this matrix will be the images of our sigma points in the lidar measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // z_pred and S represent the gaussian distribution of our beief, expressed in lidar measurement terms
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z,n_z);

  // Map sigma points over to lidar measurement space
  for (int i = 0; i< 2*n_aug_+1; ++i)
      Zsig.col(i) = h_lidar(Xsig_pred_.col(i));

  // Compute z_pred and S from Zsig

  z_pred.fill(0.0);
  for (int i=0; i<2*n_aug_+1; ++i)
      z_pred += weights_(i) * Zsig.col(i);
  
  S.fill(0.0);
  for (int i=0; i<2*n_aug_+1; ++i) {
      VectorXd diff = Zsig.col(i) - z_pred;
      S += weights_(i) * diff * diff.transpose();
  }
  // Incorporate lidar measurement uncertainty into S
  S(0,0) += std_laspx_ * std_laspx_;
  S(1,1) += std_laspy_ * std_laspy_;

  // Compute cross-correlation matrix and Kalman gain

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i=0; i<2*n_aug_+1; ++i) {
      VectorXd xdiff = Xsig_pred_.col(i) - x_;
      xdiff(3) = standardize_angle(xdiff(3));
      VectorXd zdiff = Zsig.col(i) - z_pred;
      Tc += weights_(i) * xdiff * zdiff.transpose(); 
  }

  MatrixXd K = Tc * S.inverse();

  // Update state mean and covariance matrix

  x_ = x_ + K*(z - z_pred);
  P_ = P_ - K * S * K.transpose();


}

// This function is the radar measurement model, mapping from internal states into radar meas. space
static Matrix<double, 3, 1> h_radar(Matrix<double, 5, 1> x){
    const double px=x(0);
    const double py=x(1);
    const double v=x(2);
    const double psi=x(3);
    const double psidot=x(4);
    
    const double rho=sqrt(px*px+py*py);
    const double phi=atan2(py,px);
    const double rhodot=(px*cos(psi)*v + py*sin(psi)*v)/rho;
    
    Matrix<double, 3, 1> z;
    z << rho,phi,rhodot;
    return z;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  int n_z = 3;

  // This vector consists of the actual measured rho, phi, rhodot
  VectorXd z = meas_package.raw_measurements_;

  // The columns of this matrix will be the images of our sigma points in the radar measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // z_pred and S represent the gaussian distribution of our beief, expressed in radar measurement terms
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z,n_z);

  // Map sigma points over to radar measurement space
  for (int i = 0; i< 2*n_aug_+1; ++i)
      Zsig.col(i) = h_radar(Xsig_pred_.col(i));

  // Compute z_pred and S from Zsig

  z_pred.fill(0.0);
  for (int i=0; i<2*n_aug_+1; ++i) 
      z_pred += weights_(i) * Zsig.col(i);
  
  S.fill(0.0);
  for (int i=0; i<2*n_aug_+1; ++i) {
      VectorXd diff = Zsig.col(i) - z_pred;
      diff(1) = standardize_angle(diff(1));
      S += weights_(i) * diff * diff.transpose();
  }
  // Incorporate radar measurement uncertainty into S
  S(0,0) += std_radr_*std_radr_;
  S(1,1) += std_radphi_*std_radphi_;
  S(2,2) += std_radrd_*std_radrd_;

  // Compute cross-correlation matrix and Kalman gain

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i=0; i<2*n_aug_+1; ++i) {
      VectorXd xdiff = Xsig_pred_.col(i) - x_;
      xdiff(3) = standardize_angle(xdiff(3));
      VectorXd zdiff = Zsig.col(i) - z_pred;
      zdiff(1) = standardize_angle(zdiff(1));
      Tc += weights_(i) * xdiff * zdiff.transpose(); 
  }

  MatrixXd K = Tc * S.inverse();

  // Update state mean and covariance matrix

  x_ = x_ + K*(z - z_pred);
  P_ = P_ - K * S * K.transpose();

}