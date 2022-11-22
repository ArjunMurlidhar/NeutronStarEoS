#include<C:\Users\windows7\projects\NS_Integrator\nr3.h>
#include<C:\Users\windows7\projects\NS_Integrator\OdeInt.h>
#include<C:\Users\windows7\projects\NS_Integrator\stepper.h>
#include<C:\Users\windows7\projects\NS_Integrator\stepperdopr5.h>
#include <C:\Users\windows7\projects\NS_Integrator\interp_1D.h>
#include<C:\Users\windows7\projects\NS_Integrator\AP3\interpol.h>


const Doub c = 2.998e+8;
const Doub G = 6.674e-11;
const Doub c2 = pow(c,2);
EoS ap3;

struct tov {

void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {

dydx[0]= (-G*(ap3.rho(log10(y[0]))*c2 + y[0])*(y[1]*c2 + 4.0*M_PI*pow(x,3)*y[0]))/(pow(c,4)*x*(x - (2.0*G*y[1]/c2)));
dydx[1]= 4.0*M_PI*pow(x,2)*ap3.rho(log10(y[0]));
dydx[2]= 8.0*M_PI*pow(x,4)*ap3.rho(log10(y[0]))/3.0;

}
};


int main()
{    
    ap3.initialise();
    VecDoub MoI(42), cenden(42), mass(42);
    
    Doub j = 17.479;

    for(int i = 0; i <42; i++){
    Doub rhoc = pow(10.0,j); //central mass density
    Doub pc = ap3.Pres(j); // central pressure
    const Int nvar=3;
    const Doub atol=1.0e-8, rtol=1.0e-8, h1=0.01, hmin=0.0, x1=0.01;
    VecDoub ystart(nvar);
    ystart[0]=pc;
    ystart[1]=(rhoc*4.0*M_PI*pow(x1,3)/3.0);
    ystart[2]=rhoc*8.0*M_PI*pow(x1,5)/15.0;
    Output out(-1); //Dense output at 20 points plus x1.
    tov d; //Declare d as a rhs_van object.
    Odeint<StepperDopr5<tov> > ode(ystart,x1,atol,rtol,h1,hmin,out,d);
    
    ode.integrate();
    
    
    cenden[i] = j;
    MoI[i] = ystart[2]/(ystart[1]*pow(out.xsave[out.count -1],2));
    mass[i] = ystart[1]/(1.989e+30);
    
    cout <<cenden[i]<< " "<<MoI[i]<<endl;

    j = j + 0.02;
} 

    
ofstream file3("./AP3/moi_ap3.dat");
for(int i = 0; i < 42; i++)
    file3 << cenden[i] << " " << MoI[i] <<" "<<mass[i]<< endl;
file3.close();

}