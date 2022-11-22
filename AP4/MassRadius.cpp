#include<C:\Users\windows7\projects\NS_Integrator\nr3.h>
#include<C:\Users\windows7\projects\NS_Integrator\OdeInt.h>
#include<C:\Users\windows7\projects\NS_Integrator\stepper.h>
#include<C:\Users\windows7\projects\NS_Integrator\stepperdopr5.h>
#include <C:\Users\windows7\projects\NS_Integrator\interp_1D.h>
#include<C:\Users\windows7\projects\NS_Integrator\AP4\interpol.h>


const Doub c = 2.998e+8;
const Doub G = 6.674e-11;
const Doub c2 = pow(c,2);
EoS ap4;

struct tov {

void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {

dydx[0]= (-G*(ap4.rho(log10(y[0]))*c2 + y[0])*(y[1]*c2 + 4.0*M_PI*pow(x,3)*y[0]))/(pow(c,4)*x*(x - (2.0*G*y[1]/c2)));
dydx[1]= 4.0*M_PI*pow(x,2)*ap4.rho(log10(y[0]));

}
};


int main()
{    
    ap4.initialise();
    VecDoub Mass(500), Radius(500), cenden(500);
    
    Doub j = 17.479;

    for(int i = 0; i <500; i++){
    Doub rhoc = pow(10.0,j); //central mass density
    Doub pc = ap4.Pres(j); // central pressure
    const Int nvar=2;
    const Doub atol=1.0e-8, rtol=1.0e-8, h1=0.01, hmin=0.0, x1=0.01;
    VecDoub ystart(nvar);
    ystart[0]=pc;
    ystart[1]=(rhoc*4.0*M_PI*pow(x1,3)/3.0);
    Output out(-1); //Dense output at 20 points plus x1.
    tov d; //Declare d as a rhs_van object.
    Odeint<StepperDopr5<tov> > ode(ystart,x1,atol,rtol,h1,hmin,out,d);
    
    ode.integrate();
    
    
    Radius[i] = out.xsave[out.count -1]/1000.0;
    Mass[i] = ystart[1]/(1.989e30);
    cenden[i] = j;
    cout <<Radius[i]<< " "<<Mass[i]<<endl;

    j = j + 0.003;
} 

    
ofstream file3("./AP4/massradius_ap4.dat");
for(int i = 0; i < 500; i++)
    file3 << Mass[i] << " " << Radius[i] <<" "<<cenden[i]<< endl;
file3.close();

}