#include "kalman_filter.h"
#include <math.h>
#include <iostream>
#include "tools.h"


using std::vector;

using std::cout;
using std::endl;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  Tools tools;

  cout << "EKF: px, py, vx, vy" << endl;

  float rho = sqrt(px*px + py*py);
  float phi = atan2(py, px);
  float rho_dot;

  if (rho>0){
       rho_dot = (px*vx + py*vy)/rho;     

  }

  cout << "EKF: rhodot checked " << endl;
  
  VectorXd h = VectorXd(3);
  h << rho,
       phi,
       rho_dot;


  
  VectorXd y = z - h;


  if(y[1]>M_PI){
	y[1]-=2*M_PI;
  }
  if(y[1]<-M_PI){
        y[1]+=2*M_PI;
  }
  cout << "EKF: y checked " << endl;
  MatrixXd H_j = tools.CalculateJacobian(x_);

  MatrixXd Ht = H_j.transpose();


  cout << "EKF: Ht checked " << endl;
  MatrixXd S = H_j * P_ * Ht + R_;

  cout << "EKF: S checked " << endl;  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  	
  cout << "EKF: gain checked " << endl;
  x_ = x_ + (K * y);
  long x_size = x_.size();
  cout << "EKF: size checked " << endl;
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_j) * P_;

  cout << "EKF: P checked " << endl;



}
