#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using std::vector;

double normalize_atan(double theta);

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
  std_a_ = 7;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.95;

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
  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);;

  // Prediction sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
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
	  time_us_ = meas_package.timestamp_;
    // first measurement
    x_ << 0.2, 0.4, 0, 0, 0;
    P_ << 0.1, 0, 0, 0, 0,
		  		0, 0.1, 0, 0, 0,
		  		0, 0, 0.1, 0, 0,
		  		0, 0, 0, 0.1, 0,
          0, 0, 0, 0, 0.1;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      rho = meas_package.raw_measurements_[0];
      phi = meas_package.raw_measurements_[1];
      rho_dot = meas_package.raw_measurements_[2];

      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
      x_(2) = 0;
      x_(3) = 0;
      x_(4) = 0;
    }

    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      px = meas_package.raw_measurements_[0];
      py = meas_package.raw_measurements_[1];

      x_(0) = px;
      x_(1) = py;
      x_(2) = 0;
      x_(3) = 0;
      x_(4) = 0;
    }

    std::cout << "Origin x:" << "\n" << x_ << "\n";
    std::cout << "Origin P:" << "\n" << P_ << "\n";

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; 
  
  // std::cout << "Time Now: " << meas_package.timestamp_ << std::endl;
  // std::cout << "Last measurement: " << time_us_ << std::endl;
  // std::cout << "delta_t: " << delta_t << std::endl;

  if (delta_t > 0.001) {
      while (delta_t > 0.1)
        {
            const double dt = 0.05;
            Prediction(dt);
            delta_t -= dt;
        }

        Prediction(delta_t);
    // std::cout << "Prediction x: " << "\n" << x_ << "\n";
    // std::cout << "Prediction P: " << "\n" << P_ << "\n";
  }

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
    // std::cout << "Radar x: " << "\n" << x_ << "\n";
    // std::cout << "Radar P: " << "\n" << P_ << "\n";
  } 

  else
  {
    UpdateLidar(meas_package);
    // std::cout << "Lidar x: " << "\n" << x_ << "\n";
    // std::cout << "Lidar P: " << "\n" << P_ << "\n";
  }

  time_us_ = meas_package.timestamp_;

  std::cout << "psi: " << x_(3) << std::endl;
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
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_ + 1) = 0;

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Q = MatrixXd(2, 2);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;

  Q(0, 0) = std_a_ * std_a_;
  Q(0, 1) = 0;
  Q(1, 0) = 0;
  Q(1, 1) = std_yawdd_ * std_yawdd_;
  P_aug.bottomRightCorner(2, 2) = Q;

  MatrixXd A = P_aug.llt().matrixL();
//   MatrixXd A = P_aug.ldlt().matrixL();

  Xsig_aug.col(0) = x_aug;

  for(int i = 0; i < n_aug_; ++i) {
      Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
      Xsig_aug.col(i + n_aug_ + 1) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);

    //   // Normalize Xsig_aug
    //   while (Xsig_aug.col(i)(3) > M_PI) Xsig_aug.col(i)(3) -= 2.*M_PI;
    //   while (Xsig_aug.col(i)(3) < -M_PI) Xsig_aug.col(i)(3) += 2.*M_PI;
  }

  // std::cout << "current x: \n" << x_ << std::endl;
  // std::cout << "current P: \n" << P_ << std::endl;
  // std::cout << "current Q: \n" << Q << std::endl;

  // std::cout << "x_aug: \n" << x_aug << std::endl;
  // std::cout << "P_aug: \n" << P_aug << std::endl;
  // std::cout << "Xsig_aug: \n" << Xsig_aug << std::endl;

  // Predict sigma points
  double v_k, psi_k, psi_dot_k, nu_a_k, nu_psi_ddot_k;

  for (int i = 0; i < 2*n_aug_ + 1; ++i) {
        VectorXd temp_vector_x = VectorXd(5);
        VectorXd temp_vector_noise = VectorXd(5);
        
        v_k = Xsig_aug.col(i)(2);
        psi_k = Xsig_aug.col(i)(3);
        psi_dot_k = Xsig_aug.col(i)(4);
        nu_a_k = Xsig_aug.col(i)(5);
        nu_psi_ddot_k = Xsig_aug.col(i)(6);
        
        temp_vector_noise <<
            0.5 * delta_t * delta_t * cos(psi_k) * nu_a_k,
            0.5 * delta_t * delta_t * sin(psi_k) * nu_a_k,
            delta_t * nu_a_k,
            0.5 * delta_t * delta_t * nu_psi_ddot_k,
            delta_t * nu_psi_ddot_k;
            
        if (fabs(psi_dot_k) < 0.001) {
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
        
        Xsig_pred_.col(i) = Xsig_aug.col(i).head(n_x_) + temp_vector_x + temp_vector_noise;

        // // Normalize psi
        // while (Xsig_pred_.col(i)(3) > M_PI) Xsig_pred_.col(i)(3) -= 2.*M_PI;
        // while (Xsig_pred_.col(i)(3) < -M_PI) Xsig_pred_.col(i)(3) += 2.*M_PI;        

  }

  // Predict mean / covariance matrix using sigma points
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
      weights_(i) = 1 / (2 * (lambda_ + n_aug_));
  }

  //   predict state mean

  VectorXd x_mean = VectorXd(n_x_);
  x_mean.fill(0.0);

  for (int i = 0; i < 2*n_aug_ + 1; i++) {
      x_mean = x_mean + weights_(i) * Xsig_pred_.col(i);
  }
  
  x_ = x_mean;

  //   std::cout << x << "\n";
  
  //predict state covariance matrix
  MatrixXd P_cal = MatrixXd(n_x_, n_x_);
  P_cal.fill(0.0);
  
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
      // state difference
      MatrixXd x_diff = Xsig_pred_.col(i) - Xsig_pred_.col(0);
     
      //angle normalization
      while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
      while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

      P_cal = P_cal + weights_(i) * x_diff * x_diff.transpose();
  }

  P_ = P_cal;

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

  double px_z = meas_package.raw_measurements_[0];
  double py_z = meas_package.raw_measurements_[1];

  VectorXd z = VectorXd(2);

  z <<
   px_z,
   py_z;

//   std::cout << "Lidar z: " << "\n" << z << "\n";

  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  double px, py, v, psi, psi_dot;
  
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
      px = Xsig_pred_.col(i)(0);
      py = Xsig_pred_.col(i)(1);
      v = Xsig_pred_.col(i)(2);
      psi = Xsig_pred_.col(i)(3);
      psi_dot = Xsig_pred_.col(i)(4);
      
      VectorXd temp_vector = VectorXd(n_z);
      
      temp_vector(0) = px;
      temp_vector(1) = py;
      
      Zsig.col(i) = temp_vector;
  }
  
  //calculate mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  
//   std::cout << "Lidar z: " << "\n" << z << "\n";
//   std::cout << "Lidar z_pred: " << "\n" << z_pred << "\n";
//   std::cout << "Lidar z - z_pred: " << "\n" << z - z_pred << "\n";
//   std::cout << "Lidar Zsig: \n" << Zsig << "\n";

  //calculate measurement covariance matrix S
  
  MatrixXd Lidar_R = MatrixXd(n_z, n_z);
  Lidar_R <<
    std_laspx_ * std_laspx_, 0,
    0, std_laspy_ * std_laspy_;
  
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
      S = S + weights_(i) * (Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose();
  }
  
  S = S + Lidar_R;
//   std::cout << "Lidar S: \n" << S << "\n";

  MatrixXd z_diff = z - z_pred;
  MatrixXd temp_vector2 = z_diff.transpose() * S.inverse();
  // std::cout << "Temp vector \n" << temp_vector2 << "\n";
  // std::cout << "z_diff vector \n" << z_diff << "\n";
  double nis_lidar_i = (temp_vector2 * z_diff).value();
//   std::cout << "NIS Lidar: \n" << nis_lidar_i << "\n";
  nis_lidar_.push_back(nis_lidar_i);

  // std::cout << "Lidar R: \n" << Lidar_R << "\n";
  // std::cout << "Lidar S: \n" << S << "\n";

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

   //calculate cross correlation matrix
  for (int i = 0; i < 2*n_aug_ + 1; ++i) {
        // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - Xsig_pred_.col(0);
     
      //angle normalization
      while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
      while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

      Tc = Tc + weights_(i) * x_diff * (Zsig.col(i) - z_pred).transpose();
      // std::cout << "Lidar WTF Tc: \n" << Tc << "\n";
  }
  
  // std::cout << "Lidar Xsig_pred: \n" << Xsig_pred_ << "\n";
  // std::cout << "Lidar Zsig: \n" << Zsig << "\n";
  // std::cout << "Lidar Tc: \n" << Tc << "\n";
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  // std::cout << "Lidar K: \n" << K << "\n";
  
  //update state mean and covariance matrix
  x_ = x_ + K * (z - z_pred);    
  P_ = P_ - K * S * K.transpose();
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
  
  // Get RADAR measurement
  double rho = meas_package.raw_measurements_[0];
  double phi = meas_package.raw_measurements_[1];
  double rho_dot = meas_package.raw_measurements_[2];

  VectorXd z = VectorXd(3);

  z <<
   rho,
   phi,
   rho_dot;

  // std::cout << "Radar z: " << "\n" << z << "\n";

  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  double px, py, v, psi, psi_dot;
  
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
      px = Xsig_pred_.col(i)(0);
      py = Xsig_pred_.col(i)(1);
      v = Xsig_pred_.col(i)(2);
      psi = Xsig_pred_.col(i)(3);
      psi_dot = Xsig_pred_.col(i)(4);
      
      VectorXd temp_vector = VectorXd(n_z);
      
      temp_vector(0) = sqrt(px*px + py*py);
      temp_vector(1) = normalize_atan(atan2(py, px));
      // std::cout << "Radar atan: " << temp_vector(1) << "\n";
      temp_vector(2) = (px * cos(psi) * v + py * sin(psi) * v) /  sqrt(px*px + py*py);
      
      Zsig.col(i) = temp_vector;
  }
  
  //calculate mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  
//   std::cout << "Radar z: " << "\n" << z << "\n";
//   std::cout << "Radar z_pred: " << "\n" << z_pred << "\n";
//   std::cout << "Radar z - z_pred: " << "\n" << z - z_pred << "\n";
//   std::cout << "Radar Zsig: \n" << Zsig << "\n";

  //calculate measurement covariance matrix S
  
  MatrixXd Radar_R = MatrixXd(n_z, n_z);
  
  Radar_R <<
    std_radr_ * std_radr_, 0, 0,
    0, std_radphi_ * std_radphi_, 0,
    0, 0, std_radrd_ * std_radrd_;
  
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
      S = S + weights_(i) * (Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose();
  }
  
  S = S + Radar_R;
  
//   // Normalize S phi
//   for (int i = 0; i < 3; i++) {
//     S.col(i)(1) = normalize_atan(S.col(i)(1));
//   }

//   std::cout << "Radar S: \n" << S << "\n";

  MatrixXd z_diff = z - z_pred;

  MatrixXd temp_vector2 = z_diff.transpose() * S.inverse();
  // std::cout << "Temp vector \n" << temp_vector2 << "\n";
  // std::cout << "z_diff vector \n" << z_diff << "\n";
  double nis_radar_i = (temp_vector2 * z_diff).value();
//   std::cout << "NIS Radar: \n" << nis_radar_i << "\n";
  nis_radar_.push_back(nis_radar_i);

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
   //calculate cross correlation matrix
  for (int i = 0; i < 2*n_aug_ + 1; ++i) {
     // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - Xsig_pred_.col(0);
     
    //   //angle normalization
    //   while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    //   while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

      Tc = Tc + weights_(i) * x_diff * (Zsig.col(i) - z_pred).transpose();
      // std::cout << "Radar WTF Tc: \n" << Tc << "\n";
  }
  
  // std::cout << "Radar Xsig_pred: \n" << Xsig_pred_ << "\n";
  // std::cout << "Radar Zsig: \n" << Zsig << "\n";
  // std::cout << "Radar Tc: \n" << Tc << "\n";
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  // std::cout << "Radar K: \n" << K << "\n";
  
  //update state mean and covariance matrix
  x_ = x_ + K * (z - z_pred);    
  P_ = P_ - K * S * K.transpose();
}

double normalize_atan(double theta) {
	// Normalize tan_phi so it's between -pi and pi
	while (theta > M_PI || theta <= -M_PI) {
		if (theta > M_PI) {
			theta -= 2*M_PI;
		}

		if (theta <= -M_PI) {
			theta += 2*M_PI;
		}
	}

	return theta;
}