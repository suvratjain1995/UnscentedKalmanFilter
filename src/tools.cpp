#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  
  VectorXd RMSE = VectorXd(4);
  RMSE<< 0.0,0.0,0.0,0.0;
  if(estimations.size() != 0 && estimations.size() == ground_truth.size())
  {
    for(int i=0;i<estimations.size();i++)
    {
      VectorXd temp = estimations[i] - ground_truth[i];
      RMSE = RMSE.array()+temp.array()*temp.array();
    }
    RMSE = RMSE.array()/estimations.size();
    RMSE = RMSE.array().sqrt();
    cout<<"done"<<endl;
    return RMSE;
  }
  else
  {
    cout<<"error in rmse"<<endl;
    return RMSE;
  }           
  /**
  TODO:
    * Calculate the RMSE here.
  */
}