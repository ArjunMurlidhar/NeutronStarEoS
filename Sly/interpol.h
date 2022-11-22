#include <fstream>

struct EoS{

const Doub c;
const Doub c2;
VecDoub pp, ee;

EoS(): c(2.998e8), c2(8.988e16), pp(99), ee(99) {}



void initialise()
{
    
    
    ifstream infile("./Sly/sly.dat");
    
    Doub baryondensity, pressure, energydensity;
    //cout.precision(10);
    for (int i = 0; i < 99; i++)
    {
       infile >> baryondensity>>pressure>>energydensity;
       pp[i] = log10(pressure*1000.0*c2);
       ee[i] = log10(energydensity*1000.0);
       
    
    }

    infile.close();

}



Doub Pres(Doub rho){
    Spline_interp obj1(ee,pp);
    return pow(10.0,obj1.interp(rho));
    

}

Doub rho(Doub p){
    Spline_interp obj1(pp,ee);
    return pow(10.0,obj1.interp(p));
    

}







};