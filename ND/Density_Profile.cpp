#include<C:\Users\windows7\projects\NS_Integrator\nr3.h>
#include<C:\Users\windows7\projects\NS_Integrator\OdeInt.h>
#include<C:\Users\windows7\projects\NS_Integrator\stepper.h>
#include<C:\Users\windows7\projects\NS_Integrator\stepperdopr5.h>
#include <C:\Users\windows7\projects\NS_Integrator\interp_1D.h>

VecDoub pp(230), ee(230);
const Doub c = 2.998e+8;

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
const Doub G = 6.674e-11;
const Doub c = 2.998e+8;

void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
dydx[0]= (-G*(rho(y[0])*pow(c,2) + y[0])*(y[1]*pow(c,2) + 4.0*M_PI*pow(x,3)*y[0]))/(pow(c,4)*x*(x - (2.0*G*y[1]/pow(c,2))));
dydx[1]= 4.0*M_PI*pow(x,2)*rho(y[0]);
}
};

struct tov_tilda {
const Doub G = 6.674e-11;
const Doub c = 2.998e+8;
const Doub g = -G*2.0e+20/pow(c,2);

void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
dydx[0]= g*(rho_t(y[0]) + y[0])*(1.0e+10*y[1] + 2.0*M_PI*pow(x,3)*1.0e-12*y[0])/(pow(x,2)*(1.0 - (2970.1181*y[1]/x)));
dydx[1]= 2.0*M_PI*1.0e-22*pow(x,2)*rho_t(y[0]);
}
};





int main()
{    
    for (Doub i = 20.0, j = 229; i > -3.0; i = i - 0.1, j--)
    {
        pp[j] = Pressure(exp(-i));
        ee[j] = MassDensity(exp(-i));
        
    
    }
    //VecDoub radius(250);
    //VecDoub density(250);
    Doub j = 16.98;

    
    Doub rhoc = pow(10.0,j + 11.0*0.01); //central mass density
    Doub pc = Pres(rhoc); // central pressure
    //Doub pc_t = pc/(pow(c,2)*1.0e+8);
    const Int nvar=2;
    const Doub atol=1.0e-8, rtol=1.0e-8, h1=0.01, hmin=0.0, x1=0.01;
    VecDoub ystart(nvar);
    ystart[0]=pc;
    ystart[1]=(rhoc*4.0*M_PI*pow(x1,3)/3.0);//(2.0e+30);
    Output out(-1); //Dense output at 20 points plus x1.
    tov d; //Declare d as a rhs_van object.
    Odeint<StepperDopr5<tov> > ode(ystart,x1,atol,rtol,h1,hmin,out,d);
    ode.integrate();
//cout<<endl<<ystart[0]<<endl<<ystart[1]<<endl<<out.count<<endl<<ode.nok<<endl<<ode.nbad<<endl;
    VecDoub radius(out.count);
    VecDoub density(out.count);
    VecDoub press(out.count);

    for(int i = 0; i<out.count; i++){
        radius[i] = out.xsave[i]/1000.0;
        density[i] = log10(rho(out.ysave[0][i]));
        press[i] = log10(out.ysave[0][i]);
        cout << radius[i]<< " "<<density[i]<<" "<<press[i]<<endl;
    }

    /*for(int i =0; i<235; i++){
     radius[i]= out.xsave[17*i]/1000.0;
     density[i]= log10(rho(out.ysave[0][17*i]));
     cout<<radius[i]<< " "<<density[i]<<endl;
    }
    
    for(int j=235, i = 17*234 + 1; j<250; j++){
      
     radius[j]= out.xsave[i]/1000.0;
     density[j]= log10(rho(out.ysave[0][i]));
     cout<<radius[j]<< " "<<density[j]<<endl;
     i++;

    }*/

ofstream file3("densitypressprofile1.dat");
for(int i = 0; i < out.count; i++)
    file3 << radius[i] << " " << density[i] <<endl;
file3.close();

}