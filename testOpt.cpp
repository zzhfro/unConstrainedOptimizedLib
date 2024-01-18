#include "unConOptimziaton.h"
//#define SGDTEST
using namespace std;
double fTest(const Eigen::VectorXd &input)
{
  int N=input.size();
  double result=0;
  for(int i=0;i<N/2;++i)
  {
    result=result+100*pow(pow(input(2*i),2)-input(2*i+1),2)+pow(input(2*i)-1,2);
  } 
  return result;
}

double fTest2(const Eigen::VectorXd &input)
{
  return 0.01*(input(0))*(input(0))+input(1)*input(1)+10*input(2)*input(2);
}

double fTest3(const Eigen::VectorXd &x)
{
  const int n = x.size();
        double fx = 0.0;
        for (int i = 0; i < n; i += 2)
        {
            const double t1 = 1.0 - x(i);
            const double t2 = 10.0 * (x(i + 1) - x(i) * x(i));
           
            fx += t1 * t1 + t2 * t2;
        }
        return fx;
}
int main(int argc,char **argv)
{
#ifdef NEWTONTEST
  NewtonMethod NW;
  
  Optimizer Opt(&NW);
  int N=200;
  Eigen::VectorXd xInitial(N); 
  for (int i = 0; i < N; i += 2)
        {
            xInitial(i) = -5;
            xInitial(i + 1) = 1.0;
        }
  Eigen::VectorXd optimalX; 
  auto start = std::chrono::high_resolution_clock::now();
  optimalX= Opt.performOptimization(fTest3,xInitial);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << "NewtonMethod优化运行时间: " << duration.count() << " 微秒" << std::endl;
#endif
  #ifdef SGDTEST
   SGDMethod Sgd;
   Optimizer Opt2(&Sgd);
    Eigen::VectorXd xInitial3(N); 
  for (int i = 0; i < N; i += 2)
        {
            xInitial3(i) = -1.2;
            xInitial3(i + 1) = 1.0;
        }

   auto start2 = std::chrono::high_resolution_clock::now();
   optimalX = Opt2.performOptimization(fTest3,xInitial);
   auto end2 = std::chrono::high_resolution_clock::now();
   auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
   std::cout << "SGDMethod优化运行时间: " << duration2.count() << " 微秒" << std::endl;
  #endif
  
  int N=200;
  LBFGS Lbfgs;
  Optimizer Opt3(&Lbfgs);
  Eigen::VectorXd xInitial2(N); 
  Eigen::VectorXd optimalX; 
  for (int i = 0; i < N; i += 2)
        {
            xInitial2(i) = -5.0;
            xInitial2(i + 1) = 1.0;
        }
  

  auto start3 = std::chrono::high_resolution_clock::now();
  optimalX = Opt3.performOptimization(fTest3,xInitial2);
  auto end3 = std::chrono::high_resolution_clock::now();
  auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3);
  std::cout << "LBFGS优化运行时间: " << duration3.count() << " 微秒" << std::endl;
  return 0;
}