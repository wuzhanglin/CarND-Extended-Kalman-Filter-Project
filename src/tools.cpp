#include <iostream>
#include "tools.h"

using std::cout;
using std::vector;

using Eigen::VectorXd;
using Eigen::MatrixXd;

Tools::Tools()
{
}

Tools::~Tools()
{
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
    // TODO:
    // Calculate the RMSE here
    
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;
    
    if (estimations.size() == 0)
    {
        cout << "ERROR: Tools::CalculateRMSE, The estimations vector is empty" << endl;
        return rmse;
    }
    
    if (ground_truth.size() == 0)
    {
        cout << "ERROR: Tools::CalculateRMSE, The ground-truth vector is empty" << endl;
        return rmse;
    }
    
    if (estimations.size() != ground_truth.size())
    {
        cout << "ERROR: Tools::CalculateRMSE, The ground-truth and estimations vectors must have the same size" << endl;
        return rmse;
    }
    
    for (unsigned int i = 0; i < estimations.size(); ++i)
    {
        VectorXd diff = estimations[i] - ground_truth[i];
        diff = diff.array() * diff.array();
        rmse += diff;
    }
    
    rmse = rmse / estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
    // TODO:
    // Calculate a Jacobian here
    
    MatrixXd Hj(3, 4);
    if (x_state.size() != 4)
    {
        cout << "ERROR: Tools::CalculateJacobian, The state vector must have size 4" << endl;
        return Hj;
    }
    
    // Recover state parameters
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);
    
    // Pre-compute a set of terms to avoid repeated calculation
    double c1 = px*px + py*py;
    double c2 = sqrt(c1);
    double c3 = (c1 * c2);
    
    // Check division by zero
    if (fabs(c1) < 0.0001)
    {
        cout << "ERROR: Tools::CalculateJacobian, Division by zero" << endl;
        return Hj;
    }
    
    // Compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
          -(py/c1), (px/c1), 0, 0,
          py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
    
    return Hj;
}
