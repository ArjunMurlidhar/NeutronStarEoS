#include <C:\Users\windows7\projects\NS_Integrator\nr3.h>
#include <C:\Users\windows7\projects\NS_Integrator\interp_1D.h>
#include <fstream>


Doub Pressure(Doub x){
    if(x>= 0.01)
    return log10(2.0581*pow(10.0,35)*(x*pow((1.0 + pow(x,2.0)),1.0/2.0)*((2.0/3.0)*pow(x,2.0) - 1.0) + log(x + pow((1.0 + pow(x,2.0)),1.0/2.0) )));
    else
    return log10(1.0976*pow(10.0,35)*(pow(x,5.0) - (5.0/14.0)*pow(x,7.0)));
}

Doub MassDensity(Doub x){
    if(x>=0.01)
    return log10(2.2898*pow(10.0,18)*(x*pow((1.0 + pow(x,2.0)),(1.0/2.0))*(2.0*pow(x,2.0) + 1.0) - log(x + pow((1.0 + pow(x,2.0)),1.0/2.0) )));
    else 
    return log10(6.1059*pow(10.0,18)*(pow(x,3.0) + (3.0/10.0)*pow(x,5.0) - (3.0/56.0)*pow(x,7.0)));

}

Doub logPres(Doub rho, VecDoub pp, VecDoub ee){
    Spline_interp obj1(ee,pp);
    return obj1.interp(rho);
    

}

Doub logrho(Doub p, VecDoub pp, VecDoub ee){
    Spline_interp obj1(pp,ee);
    return obj1.interp(p);
    

}





int main()
{
    
    VecDoub pp(230), ee(230);
    Doub rhoc = 1.0*pow(10.0,17.0); //central mass density
    Doub pc; // central pressure


    for (Doub i = 20.0, j = 229; i > -3.0; i = i - 0.1, j--)
    {
        pp[j] = Pressure(exp(-i));
        ee[j] = MassDensity(exp(-i));
        
    
    }

    
    //for(int i = 0; i < 230; i++)
      //cout<< pp[i] << " ";
    //cout<<endl;
    //for(auto i = &ee[0]; i < &ee[ee.size()]; i++)
      //cout<< *i << " ";
    //cout<<endl;
    
    
//cout<<endl<<logPres(7.0,pp,ee)<<endl; // Pres() - pressure(massdensity)
VecDoub Rho(260);
VecDoub P(260);

for(int i=-70; i<190; i++) 
{ Rho[i+70] = 0.1*i;
  P[i+70] = logPres(Rho[i+70],pp,ee);
}

ofstream file1("logPvslogRho.dat");
for(int i = 0; i < 260; i++)
    file1 << Rho[i] << " " << P[i]<<endl;
file1.close();


    
     

    
    
return 0;
}
