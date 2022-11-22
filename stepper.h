

struct StepperBase { 
//Base class for all ODE algorithms.
Doub &x;
Doub xold; //Used for dense output.
VecDoub &y,&dydx;
Doub atol,rtol,pbound = 3.17e+28;
bool dense;
Doub hdid; //Actual stepsize accomplished by the step routine.
Doub hnext; //Stepsize predicted by the controller for the next step.
Doub EPS;
Int n,neqn; //neqn D n except for StepperStoerm.
VecDoub yout,yerr; //New value of y and error estimate.
StepperBase(VecDoub_IO &yy, VecDoub_IO &dydxx, Doub &xx, const Doub atoll,
const Doub rtoll, bool dens) : x(xx),y(yy),dydx(dydxx),atol(atoll),
rtol(rtoll),dense(dens),n(y.size()),neqn(n),yout(n),yerr(n) {}
/*Input to the constructor are the dependent variable vector y[0..n-1] and its derivative
dydx[0..n-1] at the starting value of the independent variable x. Also input are the
absolute and relative tolerances, atol and rtol, and the boolean dense, which is true
if dense output is required.*/
};