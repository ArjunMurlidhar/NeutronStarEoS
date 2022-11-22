#include<C:\Users\windows7\projects\NS_Integrator\nr3.h>
#include<C:\Users\windows7\projects\NS_Integrator\OdeInt.h>
#include<C:\Users\windows7\projects\NS_Integrator\stepper.h>
#include<C:\Users\windows7\projects\NS_Integrator\stepperdopr5.h>
#include <C:\Users\windows7\projects\NS_Integrator\interp_1D.h>
#include<C:\Users\windows7\projects\NS_Integrator\AP3\interpol.h>


const Doub c = 2.998e+8;
EoS ap3;

struct tov {
const Doub G = 6.674e-11;
const Doub c = 2.998e+8;

void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
dydx[0]= (-G*(ap3.rho(log10(y[0]))*pow(c,2) + y[0])*(y[1]*pow(c,2) + 4.0*M_PI*pow(x,3)*y[0]))/(pow(c,4)*x*(x - (2.0*G*y[1]/pow(c,2))));
dydx[1]= 4.0*M_PI*pow(x,2)*ap3.rho(log10(y[0]));
}
};


int main()
{    
    ap3.initialise();
    
    Doub j = 17.5;

    
    Doub rhoc = pow(10.0,j); //central mass density
    Doub pc = ap3.Pres(log10(rhoc)); // central pressure
    cout<<pc<<endl;
    const Int nvar=2;
    const Doub atol=1.0e-8, rtol=1.0e-8, h1=0.01, hmin=0.0, x1=0.01;
    VecDoub ystart(nvar);
    ystart[0]=pc;
    ystart[1]=(rhoc*4.0*M_PI*pow(x1,3)/3.0);
    Output out(-1); //Dense output at 20 points plus x1.
    tov d; //Declare d as a rhs_van object.
    Odeint<StepperDopr5<tov> > ode(ystart,x1,atol,rtol,h1,hmin,out,d);
    
    ode.integrate();
    VecDoub radius(out.count);
    VecDoub density(out.count);
    VecDoub press(out.count);
    
    for(int i = 0; i<out.count; i++){
        radius[i] = out.xsave[i]/1000.0;
        density[i] = log10(ap3.rho(log10(out.ysave[0][i])));
        press[i] = log10(out.ysave[0][i]);
        cout << radius[i]<< " "<<density[i]<<" "<<press[i]<<endl;
    }
cout<<ystart[1]/(1.989e30)<<endl;
    
ofstream file3("./AP3/densitypressprofile-ap3-0.23.dat");
for(int i = 0; i < out.count; i++)
    file3 << radius[i] << " " << density[i] <<endl;
file3.close();

}