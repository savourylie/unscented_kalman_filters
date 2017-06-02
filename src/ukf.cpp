#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 10;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 10;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  double rho;
  double phi;
  double rho_dot;

  // Laser components
  double px;
  double py;

  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
	  double previous_timestamp_ = measurement_pack.timestamp_;
    // first measurement
    x_ << 0, 0, 0, 0, 0;
    P_ << 1, 0, 0, 0, 0,
		  		0, 1, 0, 0, 0,
		  		0, 0, 1, 0, 0,
		  		0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      rho = measurement_pack.raw_measurements_[0];
      phi = measurement_pack.raw_measurements_[1];
      rho_dot = measurement_pack.raw_measurements_[2];

      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
      x_(2) = 0;
      x_(3) = 0;
      x_(4) = 0;
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];

      x_(0) = px;
      x_(1) = py;
      x_(2) = 0;
      x_(3) = 0;
      x_(4) = 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  double delta_t = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; 
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Generate sigma points
  VectorXd x_aug = VectorXd(7);
  x_aug.head(5) = x;
  x_aug(5) = 0;
  x_aug(5 + 1) = 0;

  MatrixXd P_aug = MatrixXd(7, 7);
  MatrixXd Q = MatrixXd(2, 2);
  MatrixXd Xsig_aug = MatrixXd(7, 2 * 7 + 1);
  MatrixXd Xsig_pred = MatrixXd(5, 2 * 7 + 1);

  P_aug.setZero();
  P_aug.topLeftCorner(n_x, n_x) = P;

  Q(0, 0) = pow(std_a, 2);
  Q(0, 1) = 0;
  Q(1, 0) = 0;
  Q(1, 1) = pow(std_yawdd, 2);

  MatrixXd A = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug;

  for(int i = 0; i < n_aug; ++i) {
      Xsig_aug.col(i + 1) = x_aug + sqrt(lambda + n_aug) * A.col(i);
      Xsig_aug.col(i + n_aug + 1) = x_aug - sqrt(lambda + n_aug) * A.col(i);
  }

  // Predict sigma points
  double v_k, psi_k, psi_dot_k, nu_a_k, nu_psi_ddot_k;

  for (int i = 0; i < 2*n_aug + 1; ++i) {
        VectorXd temp_vector_x = VectorXd(5);
        VectorXd temp_vector_noise = VectorXd(5);
        
        v_k = Xsig_aug.col(i)(2);
        psi_k = Xsig_aug.col(i)(3);
        psi_dot_k = Xsig_aug.col(i)(4);
        nu_a_k = Xsig_aug.col(i)(5);
        nu_psi_ddot_k = Xsig_aug.col(i)(6);
        
        temp_vector_noise <<
            0.5 * pow(delta_t, 2) * cos(psi_k) * nu_a_k,
            0.5 * pow(delta_t, 2) * sin(psi_k) * nu_a_k,
            delta_t * nu_a_k,
            0.5 * pow(delta_t, 2) * nu_psi_ddot_k,
            delta_t * nu_psi_ddot_k;
            
        if (Xsig_aug.col(i)(4) == 0) {
            temp_vector_x <<
            v_k * cos(psi_k) * delta_t,
            v_k * sin(psi_k) * delta_t,
            0,
            0,
            0;
            
        }
      
        else {
            temp_vector_x <<
            (v_k / psi_dot_k) * (sin(psi_k + psi_dot_k * delta_t) - sin(psi_k)),
            (v_k / psi_dot_k) * (-cos(psi_k + psi_dot_k * delta_t) + cos(psi_k)),
            0,
            psi_dot_k * delta_t,
            0;
        }
        
        Xsig_pred.col(i) = Xsig_aug.col(i).head(n_x) + temp_vector_x + temp_vector_noise;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
