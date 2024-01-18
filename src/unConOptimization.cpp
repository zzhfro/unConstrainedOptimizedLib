//namespcae uncon to do 
#define LimitedMemory
//#define LbfgsDebug
#include "unConOptimziaton.h"
#include "baseAlgorithm.h"
using namespace std;
 Eigen::VectorXd SGDMethod::optimize(fToSolved f,const Eigen::VectorXd &xInitial) 
 {
   double c=0.02,alpha=1;
   double iteration=0;
   Eigen::VectorXd optimalX= xInitial;
   Eigen::VectorXd gradient(optimalX.size());
   while(!IfConvergence(f,optimalX))
   {
     gradient=getGradient(f,optimalX);
     if(f(optimalX)-f(optimalX-alpha*gradient)<-c*alpha*gradient.dot(gradient))
     {
       alpha=alpha/2; 
     }
     else
     {
        ++iteration;
        optimalX=optimalX-alpha*gradient;
#ifdef SGD;
         cout<<"------------"<<endl;
         cout<<optimalX<<endl;
         cout<<iteration<<endl;
#endif;
     }
   }
   cout<<"----final result--------"<<endl;
   cout<<optimalX<<endl;
   cout<<"----iteration number"<<endl;
   cout<<iteration<<endl;
   cout<<"-----end-----"<<endl;
   return optimalX;
 }

 Eigen::VectorXd NewtonMethod::optimize(fToSolved f,const Eigen::VectorXd &xInitial) 
 {
   double iteration=0;
   Eigen::VectorXd optimalX= xInitial;
   Eigen::VectorXd gradient(optimalX.size());
   while(!IfConvergence(f,optimalX))
   {
     gradient=getGradient(f,optimalX);
     
     Eigen::MatrixXd Hessian=getHessian(f,optimalX);
     Eigen::MatrixXd HessianInverse=Hessian.inverse();
     optimalX=optimalX-HessianInverse*gradient;
     ++iteration;

#ifdef Newton;
     cout<<"------------"<<endl;
     cout<<optimalX<<endl;
     cout<<iteration<<endl;
#endif;
   }
   cout<<"----final result--------"<<endl;
   cout<<optimalX.transpose()<<endl;
   cout<<"----iteration number"<<endl;
   cout<<iteration<<endl;
   cout<<"-----end-----"<<endl;
   return optimalX;
 }
 inline int LBFGS::LOvertonLsearch(fToSolved f, Eigen::VectorXd &optimalX,const Eigen::VectorXd &d,double c1,double c2)
 {
  double l=0,u=1e20,alpha=1.0;
  double fX=f(optimalX),fXPast=fX;
  Eigen::VectorXd x=optimalX;
  Eigen::VectorXd gradient=getGradient(f,optimalX);
  Eigen::VectorXd gradientPast=gradient;
  while(true)
  {
    //更新x fx,g
     optimalX=x+alpha*d;
    fX=f(optimalX);
    gradient=getGradient(f,optimalX);
    if((fXPast-fX)<-c1*alpha*d.dot(gradientPast))
    {
      u=alpha;
    }
    else if(d.dot(gradient)<c2*d.dot(gradientPast))
    {
      l=alpha;
    }
    else
    {
      break;
    }
    if(u<1e20)
    {
      alpha=(l+u)*0.5;
    }
    else
    {
      alpha=2*l;
    }
  }
  return 1;
 }
 inline Eigen::MatrixXd LBFGS::CautiousUpate(const Eigen::VectorXd &gradientK,const Eigen::VectorXd &gradientKPlus,const Eigen::VectorXd &deltaX ,const Eigen::MatrixXd &B)
 {
   double theta=1e-6;
   Eigen::MatrixXd I = Eigen::MatrixXd::Identity(B.rows(), B.cols());
   Eigen::VectorXd deltaGadient=gradientKPlus-gradientK;
   if(deltaGadient.dot(deltaX)>theta*gradientK.norm()*deltaX.dot(deltaX))
   {
    return (I -deltaX*deltaGadient.transpose()
           /deltaGadient.dot(deltaX))*B*(I-deltaGadient*deltaX.transpose()/(deltaGadient.dot(deltaX)))+
           deltaX*deltaX.transpose()/deltaGadient.dot(deltaX);
   }
   else
   {
     return B;
   }
 }
 /**
  * d<-limitedMemory(gk,xk) //cautious  opdate
  * t<-inexact linesearch
  * x update->xk
  * g update->gk
  * k++
  * 
  * notation:s(k)=x(k+1)-xk
  *          y(k)=g(k+1)-gk
  *          p(k)=1/<s(k),y(k)>
 */
 
 Eigen::VectorXd LBFGS::optimize(fToSolved f,const Eigen::VectorXd &xInitial)
 {
   double iteration=0;
   
   Eigen::VectorXd optimalX= xInitial;
   Eigen::VectorXd gradientKPlus=getGradient(f,optimalX);
   Eigen::MatrixXd B=Eigen::MatrixXd::Identity(optimalX.size(),optimalX.size());
   Eigen::VectorXd gradientK=getGradient(f,optimalX);
   Eigen::VectorXd d=-getGradient(f,optimalX);
   
   //variables,limited memory
   int endPosition=0;int currentEndPosition=0,position=0;
   int maxLength=8,currentLength=0;


   
   Eigen::MatrixXd s=Eigen::MatrixXd::Zero(maxLength,optimalX.size());
   Eigen::MatrixXd y=Eigen::MatrixXd::Zero(maxLength,optimalX.size());
   Eigen::VectorXd p=Eigen::VectorXd::Zero(maxLength);
   Eigen::VectorXd a=Eigen::VectorXd::Zero(maxLength);
   Eigen::VectorXd pastd;
   //std::vector<Eigen::VectorXd> s(maxLength),y(maxLength);
   //std::vector<double> p(maxLength),a(maxLength);
   
   //label the last element's position in the vector
   
   double updatebound=1e-6;
   Eigen::VectorXd lastOptimalX,lastGradient;
   Eigen::VectorXd deltaX,deltaGradient;
   
   //something about debug funciton 
   // it will be a general funciton in the future
   char userInput;
   int status;
   while(true)
   {
#ifdef LimitedMemory

        
       /*
        std::cout << "Press Enter to execute optimized function (Q to quit): ";
        std::cin.get(userInput);
        if(userInput=='\n')
        {
          */
          
       
       
     //lbfgs do not store Bk explicitly
      
      
      lastOptimalX=optimalX,lastGradient=gradientK;
      status=LOvertonLsearch(f,optimalX,d);
      gradientK=getGradient(f,optimalX);
      
      /**
       * the ifConvergence function will test the inf norm of the cost function
       * so if it pass the test break
      */
      if(IfConvergence(f,optimalX))
      {
        break;
      }

      deltaX=optimalX-lastOptimalX,deltaGradient=gradientK-lastGradient;
      
      /**
       * if satisfy cautious update store this       
      */
     endPosition=endPosition%maxLength;
     currentEndPosition=endPosition;
     d=-gradientK;
     
    if(deltaGradient.dot(deltaX)>updatebound*gradientK.norm()*deltaX.squaredNorm())
    { 
       
     s.row(currentEndPosition)=deltaX,y.row(currentEndPosition)=deltaGradient;
     p(currentEndPosition)=1/(s.row(currentEndPosition).dot(y.row(currentEndPosition)));
        
     ++endPosition;
     ++currentLength;
     currentLength=currentLength<maxLength?currentLength:maxLength;
    }
      //get d
      for(int i=0;i<currentLength;++i)
      {
        position=(currentEndPosition+currentLength-i)%currentLength;   
        a(position)=p(position)*(s.row(position).dot(d));
        d+=(-a(position))*y.row(position);      
      }
      
      double gama=p(currentEndPosition)*(y.row(currentEndPosition).squaredNorm());
      d*=1/gama;

      for(int i=0;i<currentLength;++i)
      {

        position=(currentEndPosition+i+1)%currentLength;  
        double beta=p(position)*(y.row(position).dot(d));
        d+=s.row(position)*(a(position)-beta);
      }  

      
     

        
      ++iteration;    

#else
       if(IfConvergence(f,optimalX))
      {
        break;
      }
      gradientK=getGradient(f,optimalX);
      d=-B*gradientK;
      double alpha=LOvertonLsearch(f,optimalX,d);
      optimalX=optimalX+alpha*d;
      gradientKPlus=getGradient(f,optimalX);
      B=CautiousUpate(gradientK,gradientKPlus,alpha*d,B);
      ++iteration;
#endif

#ifdef LbfgsDebug;
     cout<<"------LbfgsDebug------"<<endl;
     cout<<optimalX.transpose()<<endl;
     cout<<iteration<<endl;
     cout<<"gradientINfnorm"<<endl;
     cout<<gradientK.cwiseAbs().maxCoeff()<<endl;
#endif;
      // }
  }
   cout<<"----LBFGS result--------"<<endl;
   cout<<optimalX.transpose()<<endl;
   cout<<"----iteration number"<<endl;
   cout<<iteration<<endl;
   cout<<"-----end-----"<<endl;
   return optimalX;
}
