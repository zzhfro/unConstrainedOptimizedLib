#ifndef BASE_ALG_H
#define BASE_ALG_H
#include <iostream>
#include<Eigen/Eigen> 

 Eigen::VectorXd getGradient(double (*fToSolved)(const Eigen::VectorXd &input),const Eigen::VectorXd &xPosition);
 Eigen::MatrixXd getHessian(double (*fToSolved)(const Eigen::VectorXd &input),const Eigen::VectorXd &xPosition); 
 bool IfConvergence(double (*fToSolved)(const Eigen::VectorXd &input),const Eigen::VectorXd &xPosition);




#endif