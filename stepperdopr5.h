#include<C:\Users\windows7\projects\NS_Integrator\nr3.h>

template <class D> 
struct StepperDopr5 : StepperBase {
/*Dormand-Prince fifth-order Runge-Kutta step with monitoring of local truncation error to ensure
accuracy and adjust stepsize.*/
typedef D Dtype; //Make the type of derivs available to odeint.
VecDoub k2,k3,k4,k5,k6;
VecDoub rcont1,rcont2,rcont3,rcont4,rcont5;
VecDoub dydxnew;
StepperDopr5(VecDoub_IO &yy, VecDoub_IO &dydxx, Doub &xx,
const Doub atoll, const Doub rtoll, bool dens);
void step(const Doub htry,D &derivs);
void dy(const Doub h,D &derivs);
void prepare_dense(const Doub h,D &derivs);
Doub dense_out(const Int i, const Doub x, const Doub h);
Doub error();
struct Controller {
Doub hnext,errold;
bool reject;
Controller();
bool success(const Doub err, Doub &h, VecDoub yout, Doub pbound);
};
Controller con;
};

template <class D> 
StepperDopr5<D>::StepperDopr5(VecDoub_IO &yy,VecDoub_IO &dydxx,Doub &xx,
const Doub atoll,const Doub rtoll,bool dens) :
StepperBase(yy,dydxx,xx,atoll,rtoll,dens), k2(n),k3(n),k4(n),k5(n),k6(n),
rcont1(n),rcont2(n),rcont3(n),rcont4(n),rcont5(n),dydxnew(n) {
/*Input to the constructor are the dependent variable y[0..n-1] and its derivative dydx[0..n-1]
at the starting value of the independent variable x. Also input are the absolute and relative
tolerances, atol and rtol, and the boolean dense, which is true if dense output is required.*/
EPS=numeric_limits<Doub>::epsilon();
}

template <class D> 
void StepperDopr5<D>::step(const Doub htry,D &derivs) {
/*Attempts a step with stepsize htry. On output, y and x are replaced by their new values, hdid
is the stepsize that was actually accomplished, and hnext is the estimated next stepsize.*/
//int i =1;
Doub h=htry; //Set stepsize to the initial trial value.
for (;;) {
    
   //if(i == 20) break; 
   dy(h,derivs);           //Take a step.
   Doub err=error(); 
            //Evaluate accuracy.
   if (con.success(err,h,yout,pbound)) break;  
   //cout<< "step rejected "<<endl;                  //Step rejected. Try again with reduced h set by controller.
   if (abs(h) <= abs(x)*EPS) 
    throw("stepsize underflow in StepperDopr5");
   //i++;
}
if (dense)                     //Step succeeded. Compute coefficients for dense output.
prepare_dense(h,derivs); 
dydx=dydxnew;                  //Reuse last derivative evaluation for next step.
y=yout;
xold=x;                         //Used for dense output.
x += (hdid=h);
hnext=con.hnext;
}

template <class D>
void StepperDopr5<D>::dy(const Doub h,D &derivs) {
/*Given values for n variables y[0..n-1] and their derivatives dydx[0..n-1] known at x, use the
fifth-order Dormand-Prince Runge-Kutta method to advance the solution over an interval h and
store the incremented variables in yout[0..n-1]. Also store an estimate of the local truncation
error in yerr using the embedded fourth-order method.*/
static const Doub c2=0.2,c3=0.3,c4=0.8,c5=8.0/9.0,a21=0.2,a31=3.0/40.0,
a32=9.0/40.0,a41=44.0/45.0,a42=-56.0/15.0,a43=32.0/9.0,a51=19372.0/6561.0,
a52=-25360.0/2187.0,a53=64448.0/6561.0,a54=-212.0/729.0,a61=9017.0/3168.0,
a62=-355.0/33.0,a63=46732.0/5247.0,a64=49.0/176.0,a65=-5103.0/18656.0,
a71=35.0/384.0,a73=500.0/1113.0,a74=125.0/192.0,a75=-2187.0/6784.0,
a76=11.0/84.0,e1=71.0/57600.0,e3=-71.0/16695.0,e4=71.0/1920.0,
e5=-17253.0/339200.0,e6=22.0/525.0,e7=-1.0/40.0;
VecDoub ytemp(n);
Int i;
//cout<<yout[0]<<" "<<yout[1]<<" "<<x<<" "<<h<<" "<<endl;
for (i=0;i<n;i++)              //First step.
ytemp[i]=y[i]+h*a21*dydx[i];
derivs(x+c2*h,ytemp,k2);       //Second step.
for (i=0;i<n;i++)
ytemp[i]=y[i]+h*(a31*dydx[i]+a32*k2[i]);
//cout<< ytemp[0] << " "<< ytemp[1] << " "<< h << endl;
derivs(x+c3*h,ytemp,k3);       //Third step.
for (i=0;i<n;i++)
ytemp[i]=y[i]+h*(a41*dydx[i]+a42*k2[i]+a43*k3[i]);
//cout<< ytemp[0] << " "<< ytemp[1] << " "<< h << endl;
derivs(x+c4*h,ytemp,k4);          //Fourth step.
for (i=0;i<n;i++)
ytemp[i]=y[i]+h*(a51*dydx[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]);
//cout<< ytemp[0] << " "<< ytemp[1] << " "<< h << endl;
derivs(x+c5*h,ytemp,k5);         //Fifth step.
for (i=0;i<n;i++)
ytemp[i]=y[i]+h*(a61*dydx[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]);
//cout<< ytemp[0] << " "<< ytemp[1] << " "<< h << endl;
Doub xph=x+h;
derivs(xph,ytemp,k6);            //Sixth step.
for (i=0;i<n;i++)               //Accumulate increments with proper weights.
yout[i]=y[i]+h*(a71*dydx[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
//cout<< ytemp[0] << " "<< ytemp[1] << " "<< h << endl;
derivs(xph,yout,dydxnew);      //Will also be first evaluation for next step.
for (i=0;i<n;i++) {
//Estimate error as difference between fourth- and fifth-order methods.
yerr[i]=h*(e1*dydx[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*dydxnew[i]);
}
}

template <class D>
void StepperDopr5<D>::prepare_dense(const Doub h,D &derivs) {
//Store coefficients of interpolating polynomial for dense output in rcont1...rcont5.
VecDoub ytemp(n);
static const Doub d1=-12715105075.0/11282082432.0,
d3=87487479700.0/32700410799.0, d4=-10690763975.0/1880347072.0,
d5=701980252875.0/199316789632.0, d6=-1453857185.0/822651844.0,
d7=69997945.0/29380423.0;
for (Int i=0;i<n;i++) {
rcont1[i]=y[i];
Doub ydiff=yout[i]-y[i];
rcont2[i]=ydiff;
Doub bspl=h*dydx[i]-ydiff;
rcont3[i]=bspl;
rcont4[i]=ydiff-h*dydxnew[i]-bspl;
rcont5[i]=h*(d1*dydx[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+
d7*dydxnew[i]);
}
}

template <class D> 
Doub StepperDopr5<D>::dense_out(const Int i,const Doub x,const Doub h) {
//Evaluate interpolating polynomial for y[i] at location x, where xold <= x <= xold + h.
Doub s=(x-xold)/h;
Doub s1=1.0-s;
return rcont1[i]+s*(rcont2[i]+s1*(rcont3[i]+s*(rcont4[i]+s1*rcont5[i])));
}

template <class D> 
Doub StepperDopr5<D>::error() {
//Use yerr to compute norm of scaled error estimate. A value less than one means the step was successful.
Doub err=0.0,sk;
for (Int i=0;i<n;i++) {
sk=atol+rtol*MAX(abs(y[i]),abs(yout[i]));
err += SQR(yerr[i]/sk);
}
return sqrt(err/n);
}

template <class D> 
StepperDopr5<D>::Controller::Controller() : reject(false), errold(1.0e-4) {}
//Step size controller for fifth-order Dormand-Prince method.
template <class D>
bool StepperDopr5<D>::Controller::success(const Doub err,Doub &h,VecDoub yout,Doub pbound) {
/*Returns true if err <= 1, false otherwise. If step was successful, sets hnext to the estimated
optimal stepsize for the next step. If the step failed, reduces h appropriately for another try.*/
static const Doub beta=0.04,alpha=0.2-beta*0.75,safe=0.9,minscale=0.2,
maxscale=10.0;
//Set beta to a nonzero value for PI control. beta = 0:04–0.08 is a good default.
Doub scale;
if (err <= 1.0 && yout[0] >= -pbound) {                         //Step succeeded. Compute hnext.
if (err == 0.0)
scale=maxscale;
else {                                     //PI control if beta != 0.
scale=safe*pow(err,-alpha)*pow(errold,beta);
if (scale<minscale) scale=minscale;        //Ensure minscale <= hnext/h <= maxscale.
if (scale>maxscale) scale=maxscale;
}
if (reject)                                //Don’t let step increase if last one was rejected.
hnext=h*MIN(scale,1.0); 
else
hnext=h*scale;
errold=MAX(err,1.0e-4);                    //Bookkeeping for next call.
reject=false;
return true;
} 
else if(yout[0] < -pbound){
    scale=minscale;
    h *= scale;
    reject=true;
    return false;

}

else {                                            //Truncation error too large, reduce stepsize.
    scale=MAX(safe*pow(err,-alpha),minscale);
    h *= scale;
    reject=true;
    return false;
}
}                       

