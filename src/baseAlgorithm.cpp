#include "baseAlgorithm.h"
/**
 * zzh 2024.1
 * This file contains the implementation of some basic functions
 * such as get the derivation,hessian matrix of the function ,convergence judgment,etc.
 * 
*/
Eigen::VectorXd getGradient(double (*fToSolved)(const Eigen::VectorXd &input),const Eigen::VectorXd &xPosition)
{
   double delta=1e-10;
   int n=xPosition.size();
   Eigen::VectorXd gradResult(n);
   for(int i=0;i<n;++i)
   {
        Eigen::VectorXd xPlus= xPosition;
        Eigen::VectorXd xMinus = xPosition;
        xPlus(i) =xPlus(i)+ delta;
        xMinus(i)=xMinus(i)-delta;
        gradResult(i) =(fToSolved(xPlus) - fToSolved(xMinus)) / (2 * delta);
   }
   return gradResult;
}
/**
 * 
*/
Eigen::MatrixXd getHessian(double (*fToSolved)(const Eigen::VectorXd &input), const Eigen::VectorXd &xPosition) 
{
    int n = xPosition.size();
    Eigen::MatrixXd hessian(n, n);
    double epsilon = 1e-5;
    double f = fToSolved(xPosition);
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            Eigen::VectorXd x1 = xPosition;
            Eigen::VectorXd x2 = xPosition;
            Eigen::VectorXd x3 = xPosition;
            Eigen::VectorXd x4 = xPosition;
            
            if(i==j)
            {         
             x1(i)=x1(i)+epsilon;
             x2(i)=x2(i)-epsilon; 
             hessian(i, j) = (fToSolved(x1) + fToSolved(x2) - 2 * fToSolved(xPosition)) / (epsilon * epsilon);
            }
            else
            {
             x1(i)=x1(i)+epsilon,x1(j)=x1(j)+epsilon;
             x2(i)=x2(i)-epsilon,x2(j)=x2(j)-epsilon; 
             x3(i)=x3(i)+epsilon,x3(j)=x3(j)-epsilon;
             x4(i)=x4(i)-epsilon,x4(j)=x4(j)+epsilon; 
             hessian(i, j)=(fToSolved(x1) + fToSolved(x2) - fToSolved(x3) - fToSolved(x4)) / (4*epsilon * epsilon);
            }
            hessian(j, i) = hessian(i, j);
        }
    }

    return hessian;
}

//currently it can no be applied to non convex function
 bool IfConvergence(double (*fToSolved)(const Eigen::VectorXd &input),const Eigen::VectorXd &xPosition)
{
  Eigen::VectorXd gradResult=getGradient(fToSolved,xPosition);
  double gnorm=gradResult.cwiseAbs().maxCoeff();
  //double xnorm=xPosition.cwiseAbs().maxCoeff();
  if(gnorm <=1e-7)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}