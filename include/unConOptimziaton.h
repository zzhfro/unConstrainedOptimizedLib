#ifndef MyOptimzaiton_H
#define MyOptimzaiton_H

#include <iostream>
#include <Eigen/Eigen>
#include <chrono>
using namespace std;

typedef double(* fToSolved)(const Eigen::VectorXd &input);


class OptimizationMethod {
public:
    virtual Eigen::VectorXd  optimize(fToSolved f, const Eigen::VectorXd &xInitial){} ;

    virtual ~OptimizationMethod() {}
};

class  SGDMethod: public OptimizationMethod {
public:
    Eigen::VectorXd optimize(fToSolved f, const Eigen::VectorXd &xInitial) override ;
};

class  NewtonMethod: public OptimizationMethod {
public:
    Eigen::VectorXd optimize(fToSolved f, const Eigen::VectorXd &xInitial) override ;
};
/**
 * quasi-NewtonMethod is much more practical
 * 
*/
class LBFGS:public OptimizationMethod{
public:
    int LOvertonLsearch(fToSolved f, Eigen::VectorXd &optimalX,const Eigen::VectorXd &d,double c1=1e-4,double c2=0.9);
    Eigen::MatrixXd CautiousUpate(const Eigen::VectorXd &gradientK,const Eigen::VectorXd &gradientKPlus,const Eigen::VectorXd &deltaX ,const Eigen::MatrixXd &B);
    Eigen::VectorXd optimize(fToSolved f, const Eigen::VectorXd &xInitial) override ;
};



class Optimizer {
private:
    OptimizationMethod* method;
public:
    Optimizer(OptimizationMethod* optMethod) : method(optMethod) {}

    Eigen::VectorXd performOptimization(fToSolved f, Eigen::VectorXd xInitial) 
    {
        return method->optimize(f, xInitial);
    }
};

#endif