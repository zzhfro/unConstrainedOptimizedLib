/**
 * this code is written for practicing optimization also a small optiamziation library 
 * 
 * zzh 2023.9
 * i use eigen for review and matplotlib (python) for visualization
*/
#include <iostream>
//#include <Python.h> //i will write a python script for visualization
#include <Eigen/Eigen>
//#include "optmMethod.cpp" //just fot convenience i make fix it in the future
using namespace std;

double f(Eigen::VectorXd input)
{
  int N=input.size();
  double result=0;
  for(int i=0;i<N/2;++i)
  {
    result=result+100*pow(pow(input(2*i),2)-input(2*i+1),2)+pow(input(2*i)-1,2);
  } 
  return result;
}
Eigen::VectorXd getGradient(double (*fToSolved)(Eigen::VectorXd input),const Eigen::VectorXd xPosition,double delta=1e-5)
{
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
bool IfConvergence(double (*fToSolved)(Eigen::VectorXd input),const Eigen::VectorXd xPosition)
{
  Eigen::VectorXd gradResult=getGradient(fToSolved,xPosition);
  double norm=gradResult.norm();
  if(norm<1e-5)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

int main(int argc,char **argv)
{  
  Eigen::VectorXd xInitial(4); 
  xInitial<<0,0,0,0;
  Eigen::VectorXd optimalX=xInitial; 
  
  double alpha=0.2;  // Initialize alpha as a small positive value
  double c=0.1;   // Chooses a small value for c to adjust the step size alpha
  cout<<getGradient(f,xInitial)<<endl;
  
  while(! IfConvergence(f,optimalX))
  {
   Eigen::VectorXd grad=getGradient(f,optimalX); 
   double fPlus=f(optimalX-alpha*grad);
   if(fPlus>f(xInitial)+c*alpha*grad.dot(grad)) // check if the update would increase the function
   {
    alpha=alpha/2;  // if so, reduce step size by half
   }
   else
   {
    optimalX=optimalX-alpha*grad; // update the position by moving in the direction opposite to the gradient
   }
  }
   cout<<getGradient(f,optimalX)<<endl;
  cout<<optimalX<<endl;
return 0;
}



