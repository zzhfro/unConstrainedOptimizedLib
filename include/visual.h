#ifndef VISUAL_H
#define VISUAL_H


#include <cstdlib>
#include <iostream>
#include<Eigen/Eigen> 

typedef double(* fToSolved)(Eigen::VectorXd input);
class visual
{
  public:

  void visualOpt(fToSolved f,Eigen::MatrixXd);

  
};



#endif