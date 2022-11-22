#include<C:\Users\windows7\projects\NS_Integrator\nr3.h>
#include<C:\Users\windows7\projects\NS_Integrator\OdeInt.h>
#include<C:\Users\windows7\projects\NS_Integrator\stepper.h>
#include<C:\Users\windows7\projects\NS_Integrator\stepperdopr5.h>
#include <C:\Users\windows7\projects\NS_Integrator\interp_1D.h>

VecDoub pp(230), ee(230);
const Doub c = 2.998e+8;
const Doub G = 6.674e-11;
const Doub c2 = pow(c,2);
const Doub c4 = pow(c,4);
const Doub m = 1.66464;

Doub Pressure(Doub x){
    if(x>= 0.01)
    return 2.0581*pow(10.0,35)*(x*pow((1.0 + pow(x,2.0)),1.0/2.0)*((2.0/3.0)*pow(x,2.0) - 1.0) + log(x + pow((1.0 + pow(x,2.0)),1.0/2.0) ));
    else
    return 1.0976*pow(10.0,35)*(pow(x,5.0) - (5.0/14.0)*pow(x,7.0));
}

Doub MassDensity(Doub x){
    if(x>=0.01)
    return 2.2898*pow(10.0,18)*(x*pow((1.0 + pow(x,2.0)),(1.0/2.0))*(2.0*pow(x,2.0) + 1.0) - log(x + pow((1.0 + pow(x,2.0)),1.0/2.0) ));
    else 
    return 6.1059*pow(10.0,18)*(pow(x,3.0) + (3.0/10.0)*pow(x,5.0) - (3.0/56.0)*pow(x,7.0));

}


Doub rho(Doub p){
    Spline_interp obj1(pp,ee);
    return obj1.interp(p);
    

}

Doub rho_t(Doub p_t){
    Spline_interp obj1(pp,ee);
    return ((obj1.interp(pow(c,2)*1.0e+8*p_t))/1.0e+8);
    

}

Doub Pres(Doub rho){
    Spline_interp obj1(ee,pp);
    return obj1.interp(rho);
    

}

struct tov {

Doub f(const Doub x, VecDoub_I y)
{
    Doub elambda = 1.0/(1.0 - (2.0*G*y[1]/(x*c2)));
    return (2.0/x) + (elambda*G/c4)*((2.0*y[1]*c2/pow(x,2)) + 4.0*M_PI*x*(y[0] - rho(y[0])*c2));

}

Doub g(const Doub x, VecDoub_I y)
{
    Doub elambda = 1.0/(1.0 - (2.0*G*y[1]/(x*c2)));
    Doub nudash = (2.0*G/c4)*(y[1]*c2 + 4.0*M_PI*pow(x,3)*y[0])*(elambda/pow(x,2));
    return (-6.0*elambda/pow(x,2)) + (4.0*M_PI*G*elambda/c4)*(5.0*rho(y[0])*c2 + 9.0*y[0] + ((rho(y[0])*c2 + y[0])*rho(y[0])*c2/(m*y[0]))) - pow(nudash,2);
}
void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
dydx[0]= (-G*(rho(y[0])*pow(c,2) + y[0])*(y[1]*pow(c,2) + 4.0*M_PI*pow(x,3)*y[0]))/(pow(c,4)*x*(x - (2.0*G*y[1]/pow(c,2))));
dydx[1]= 4.0*M_PI*pow(x,2)*rho(y[0]);
dydx[2] = y[3];
dydx[3] = -f(x,y)*y[3] - g(x,y)*y[2];
}
};

struct tov_tilda {
const Doub g = -G*2.0e+20/pow(c,2);

void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
dydx[0]= g*(rho_t(y[0]) + y[0])*(1.0e+10*y[1] + 2.0*M_PI*pow(x,3)*1.0e-12*y[0])/(pow(x,2)*(1.0 - (2970.1181*y[1]/x)));
dydx[1]= 2.0*M_PI*1.0e-22*pow(x,2)*rho_t(y[0]);
}
};

Doub bc(Doub h, Doub h1, Doub h2, Doub R, Doub M)
    {
        Doub elambda = 1.0/(1.0 - (2.0*G*M/(R*c2)));
        Doub lambdadash = (-2.0*G*M/c2)*(elambda/pow(R,2));
        return h2 + ((2.0/R) - lambdadash)*h1 - ((6.0*elambda/pow(R,2)) + pow(lambdadash,2))*h;

    }


int main()
{    
    for (Doub i = 20.0, j = 229; i > -3.0; i = i - 0.1, j--)
    {
        pp[j] = Pressure(exp(-i));
        ee[j] = MassDensity(exp(-i));
        
    
    }
    VecDoub faval(201);
    VecDoub aval(201);
    VecDoub H(201);
    VecDoub H1(201);
    VecDoub H2(201);
    Doub j = 16.98 + 30.0*0.01;

    
  int i = 0;  
for(Doub a = -10.0; a <= 10.0; a+=0.1){
    //Doub a = 1.0e-4;
    Doub rhoc = pow(10.0,j); //central mass density
    Doub pc = Pres(rhoc); // central pressure
    const Int nvar=4;
    const Doub atol=1.0e-8, rtol=1.0e-8, h1=0.0001, hmin=0.0, x1=0.0001;
    VecDoub ystart(nvar);
    ystart[0]=pc;
    ystart[1]=(rhoc*4.0*M_PI*pow(x1,3)/3.0);
    ystart[2]= a*pow(x1,2)*(1 - (2*M_PI*G*pow(x1,2)/(7*c4))*(5.0*rhoc*c2 + 9.0*pc + (rhoc*c2 + pc)*rhoc*c2/(m*pc)));
    ystart[3]= 2*a*x1*(1 - (4*M_PI*G*pow(x1,2)/(7*c4))*(5.0*rhoc*c2 + 9.0*pc + (rhoc*c2 + pc)*rhoc*c2/(m*pc)));
    Output out(-1); //Dense output at 20 points plus x1.
    tov d; //Declare d as a rhs_van object.
    Odeint<StepperDopr5<tov> > ode(ystart,x1,atol,rtol,h1,hmin,out,d);
    ode.integrate();
    faval[i] = bc(ystart[2],ystart[3],ode.dydx[3],out.xsave[out.count -1],ystart[1]);
    aval[i] = a;
    H[i] = ystart[2];
    H1[i] = ystart[3];
    H2[i] = ode.dydx[3];
    i++;
    cout<<out.xsave[out.count -1]<<" "<<ystart[1]<<" "<<ystart[2]<<" "<<ystart[3]<<" "<<ode.dydx[3]<<endl;
    cout<<bc(ystart[2],ystart[3],ode.dydx[3],out.xsave[out.count -1],ystart[1])<<endl;
    }

ofstream file2("Hvsa.dat");
for(int i = 0; i < 201; i++)
    file2 << aval[i] << " " << H[i]<<" "<<H1[i]<<" "<<H2[i]<<" "<<faval[i]<<endl;
file2.close();

}