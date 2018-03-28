#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  n_x = 5;
  n_aug = 7;
  lambda = 3 - n_aug;
  P_ << 1.0,0.0,0.0,0.0,0.0,
      0.0,1.0,0.0,0.0,0.0,
      0.0,0.0,1.0,0.0,0.0,
      0.0,0.0,0.0,1.0,0.0,
      0.0,0.0,0.0,0.0,1.0;

  weights = VectorXd(2*n_aug+1);
  double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  cout<<"begin"<<endl;
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (is_initialized_ == false)
  {
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    float radius = meas_package.raw_measurements_[0];
    float phi = meas_package.raw_measurements_[1];
    float rho_dot = meas_package.raw_measurements_[2];
    float px = radius*cos(phi);
    float py = radius*sin(phi);
    x_ << px,py,rho_dot,0.0,0.0;
  } 
  else
  {
    float px = meas_package.raw_measurements_[0];
    float py = meas_package.raw_measurements_[1];
    x_ << px,py,0.0,0.0,0.0;
  }
  is_initialized_ = true;
  // std::cout<<"initialized"<<std::endl;
  previous_timestamp_ = meas_package.timestamp_;
  return ;
}
  // std::cout<<"Start"<<std::endl;
  double dt = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;
  Prediction(dt);
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    //Radar State Update

    UpdateRadar(meas_package);
    std::cout<<"X State"<<std::endl;
    std::cout<<x_<<std::endl;
    std::cout<<"P Matrix"<<std::endl;
    std::cout<<P_<<std::endl;
    
  }
  else{
    UpdateLidar(meas_package);
     std::cout<<"X State"<<std::endl;
    std::cout<<x_<<std::endl;
    std::cout<<"P Matrix"<<std::endl;
    std::cout<<P_<<std::endl;
  }
}

void UKF::GenerateAugmentedSigmaPoints(MatrixXd *Xsig_aug_out,double delta_t)
{
  cout<<"AugmentedPoints"<<endl;
  //Augmeneted Sigma Points

      VectorXd x_aug = VectorXd(7);

      //create augmented state covariance
      MatrixXd P_aug = MatrixXd(7, 7);

      //create sigma point matrix
      

      MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
      Xsig_aug.fill(0.0);
    x_aug.fill(0.0);
    x_aug.head(n_x)= x_;
    x_aug.tail(2).setZero();
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x, n_x) = P_;
    P_aug.bottomRightCorner(2,2) << std_a_*std_a_,0,
                                    0,std_yawdd_*std_yawdd_;

    MatrixXd L = P_aug.llt().matrixL();
    // cout<<"Square Root P"<<endl;
    // cout<<L<<endl;

    Xsig_aug.col(0) = x_aug;

    for(int i=0;i<n_aug;i++)
    {
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda+n_aug)*L.col(i);
        Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug)*L.col(i);
    }
    // cout<<Xsig_aug<<endl;
    *Xsig_aug_out = Xsig_aug;

}


void UKF::PredictSigmaPoints(MatrixXd *Xsig_pred_out,double delta_t,MatrixXd &Xsig_aug)
{
  cout<<"PredictSigmaPoints"<<endl;
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred.fill(0.0);
  Xsig_pred = Xsig_aug.topLeftCorner(n_x,2*n_aug+1);
  MatrixXd X_correction = Xsig_aug.bottomLeftCorner(2,2*n_aug+1);
  X_correction.fill(0.0);

  for(int i=0;i<2*n_aug+1;i++)
      {
          // float yaw_rate = Xsig_pred.col(i).row(4).value();
          // float velocity = Xsig_pred.col(i).row(2).value();
          // float px = Xsig_pred.col(i).row(0).value();
          // float py = Xsig_pred.col(i).row(1).value();
          // float yaw = Xsig_pred.col(i).row(3).value();
          // float noise_a = X_correction.col(i).row(0).value();
          // float noise_yaw = X_correction.col(i).row(1).value();
          // // std::cout<<"test123" <<std::endl;
          // MatrixXd noise = MatrixXd(5,1);
          // MatrixXd rate  = MatrixXd(5,1);
          // // if(noise_a)
          // // {
          // //     std::cout<<"Noise a"<<std::endl;
          // //     std::cout<<noise_a<<std::endl;
          // // }
          // noise << (0.5)*(delta_t*delta_t)*noise_a*cos(yaw),
          //                     (0.5)*(delta_t*delta_t)*sin(yaw)*noise_a,
          //                     delta_t*noise_a,
          //                     (0.5)*(delta_t*delta_t)*noise_yaw,
          //                     delta_t*noise_yaw;
          // cout<<"noise"<<endl;
          // cout<<noise<<endl;
          // // if(noise_a)
          // // {
          // //     std::cout<<cos(yaw)<<std::endl;
          // //     std::cout<<noise<<std::endl;
          // // }
          // if(yaw_rate >0.001)
          // {
          //   cout<<"Rate"<<endl;
          //     // std::cout<<"Yaw Rate not zeros"<<std::endl;
          //     rate << (velocity/yaw_rate)*(sin(yaw+yaw_rate*delta_t) - sin(yaw)),
          //                     (velocity/yaw_rate)*(-cos(yaw+yaw_rate*delta_t) + cos(yaw)),
          //                     0,
          //                     yaw_rate*delta_t,
          //                     0;
          //     cout<<rate<<endl;
          //     Xsig_pred.col(i) = Xsig_pred.col(i)+rate+noise;
          // }
          // else
          // {
          //  cout<<"Rate"<<endl;

          //     rate << (velocity)*cos(yaw)*delta_t,
          //                     (velocity)*sin(yaw)*delta_t,
          //                     0,
          //                     yaw_rate*delta_t,
          //                     0;
          //   cout<<rate<<endl;

          //     Xsig_pred.col(i) = Xsig_pred.col(i)+rate+noise;
          // }


    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      // cout<<"delta"<<delta_t<<endl;
      // cout<<"second equation"<<endl;

        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }
    // cout<<"Interation "<<i<<endl;
    // cout<<px_p<<" "<<py_p<<endl;
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
      }
  
  // cout<<Xsig_pred<<endl;
  *Xsig_pred_out= Xsig_pred;
}


void UKF::PredictMeanCovariance(VectorXd *X,MatrixXd *PP, MatrixXd &Xsig_pred)
{
  // cout<<"PredictMeanCovariance"<<endl;
  VectorXd x = VectorXd(n_x);
  // x.fill(0.0);
  // cout<<"weights"<<endl;
  // cout<<weights<<endl;

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x, n_x);

 x.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points

    // cout<<"temp Value"<<weights(i)*Xsig_pred.col(i)<<endl;
    x = x + weights(i)*Xsig_pred.col(i);
    // cout<<"X at iter "<<i<<endl;
    // cout<<x<<endl;
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }

  *X = x;
  *PP = P;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
  
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
cout<<"Start Prediction"<<endl;

MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
Xsig_aug.fill(0.0);
GenerateAugmentedSigmaPoints(&Xsig_aug,delta_t);
MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
Xsig_pred.fill(0.0);
PredictSigmaPoints(&Xsig_pred,delta_t,Xsig_aug); 

// cout<<"X pred\n"<<Xsig_pred<<endl;
PredictMeanCovariance(&x_,&P_,Xsig_pred);

Xsig_pred_ = Xsig_pred;

// cout<<"End Prediction"<<endl;
// cout<<"P after prediction\n"<<P_<<endl;
// cout<<"X Mean \n"<<x_<<endl;

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


  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
  Zsig.fill(0.0);

    VectorXd  z = VectorXd(n_z);
    z = meas_package.raw_measurements_;
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);



  //transform sigma points into measurement space
  for(int i=0;i<n_aug*2+1;i++)
  {
      float x_state = Xsig_pred_.col(i).row(0).value();
      float y_state = Xsig_pred_.col(i).row(1).value();
      // if(x_state > 0.0001 && y_state > 0.0001)
      // {
          Zsig.col(i) << x_state,
                          y_state;

      // }
      // else
      // {
      //     float x_state = 0.0;
      //     float y_state = 0.0;
      //     Zsig.col(i) << x_state,
      //                   y_state;
                        
      // }
      
  }
  
 MatrixXd R = MatrixXd(2,2);
 R.fill(0.0);
 R(0,0) = std_laspx_*std_laspx_;
 R(1,1) = std_laspy_*std_laspy_;
 
for(int i=0;i<2*n_aug+1;i++)
{
    z_pred = z_pred + weights(i)*Zsig.col(i);
}

for(int i = 0;i< 2*n_aug+1;i++)
{
    VectorXd diff = Zsig.col(i) - z_pred;
    S = S + weights(i)*diff*diff.transpose();
}
// cout<<"S before adding R"<<endl;
  // cout<<S<<endl;

S = S + R;

// cout<<"S lidar"<<endl;
// cout<<S<<endl;

MatrixXd Tc = MatrixXd(n_x, n_z);


  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    if (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    if (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

// cout<<"Kalman gain lidar"<<endl;
// cout<<K<<endl;
  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
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

  int n_z = 3;

   MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
   Zsig.fill(0.0);

    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }
  // cout<<"Z Prediction "<<z_pred<<endl;
  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    if (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    if (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }
  // cout<<"S before adding R"<<endl;
  // cout<<S<<endl;

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

// cout<<"S radar"<<endl;
// cout<<S<<endl;

  VectorXd z = VectorXd(n_z);

  z = meas_package.raw_measurements_;

  MatrixXd Tc = MatrixXd(n_x, n_z);
Tc.fill(0.0);
   for(int i=0;i<2*n_aug+1;i++)
  {
      // VectorXd updatediff = Xsig_pred_.col(i) - x_;
      // VectorXd measurementdiff = Zsig.col(i) - z_pred;
      // if(updatediff(3)<-M_PI) updatediff(3)+=2.*M_PI;
      // if(updatediff(3)>M_PI) updatediff(3)-=2.*M_PI;
      // if(measurementdiff(1)<-M_PI) measurementdiff(1)+=2.*M_PI;
      // if(measurementdiff(1)>M_PI) measurementdiff(1)-=2.*M_PI;
      
      VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
// cout<<"Z diff \n"<<z_diff<<endl;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // cout<<"X diff \n"<<x_diff<<endl;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
   Tc = Tc + weights(i) * x_diff * z_diff.transpose();

  //  cout<<"TC \n"<<Tc<<endl;

  }

  // cout<<"Kalman Gain Before Tc*S\n"<<Tc<<endl;
  MatrixXd Kalmangain = Tc * S.inverse();
  // cout<<"kalman gain radar"<<endl;
  // cout<<Kalmangain<<endl;
  VectorXd z_diff = z-z_pred;
  if (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  if (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  x_  = x_ + Kalmangain * (z - z_pred);
  
  // std::cout << "TEst2 "<<std::endl;

  P_ = P_ - Kalmangain * S * Kalmangain.transpose();

}
