#include<C:\Users\windows7\projects\NS_Integrator\nr3.h>
#include<C:\Users\windows7\projects\NS_Integrator\OdeInt.h>
#include<C:\Users\windows7\projects\NS_Integrator\stepper.h>
#include<C:\Users\windows7\projects\NS_Integrator\stepperdopr5.h>

struct test_eqn {
//Doub eps;
//rhs_van(Doub epss) : eps(epss) {}
void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
dydx[0]= y[0] - x;
}
};

int main(){
const Int nvar=1;
const Doub atol=1.0e-5, rtol=atol, h1=0.01, hmin=0.0, x1=-1.0;
VecDoub ystart(nvar);
ystart[0]=-0.25;
Output out(0.01); //Dense output at 20 points plus x1.
test_eqn d; //Declare d as a rhs_van object.
Odeint<StepperDopr5<test_eqn> > ode(ystart,x1,atol,rtol,h1,hmin,out,d);
ode.integrate();
cout<<endl<<ystart[0]<<endl<<out.count<<endl<<ode.nok<<endl<<ode.nbad<<endl;
for (Int i=0;i<out.count;i++)
cout << out.xsave[i] << " " << out.ysave[0][i] <<endl;

}