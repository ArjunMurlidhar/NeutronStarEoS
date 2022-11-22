#include <C:\Users\windows7\projects\NS_Integrator\nr3.h>
#include <C:\Users\windows7\projects\NS_Integrator\interp_1D.h>
#include <fstream>

const Doub c = 2.998e8;
const Doub c2 = pow(c,2);

Doub logPres(Doub rho, VecDoub pp, VecDoub ee){
    Spline_interp obj1(ee,pp);
    return obj1.interp(rho);
    

}

Doub logrho(Doub p, VecDoub pp, VecDoub ee){
    Spline_interp obj1(pp,ee);
    return obj1.interp(p);
    

}

VecDoub pp(435), ee(435);



int main()
{
    
    
    ifstream infile("./alf1/alf1.dat");
    Doub baryondensity, pressure, energydensity;
    //cout.precision(10);
    for (Doub i = 0; i < 435; i++)
    {
       infile >> baryondensity>>pressure>>energydensity;
       pp[i] = log10(pressure*1000.0*c2);
       ee[i] = log10(energydensity*1000.0);
       //cout<<pp[i]<<" "<<ee[i]<<endl;
       
    }

    infile.close();
VecDoub Rho(311);
VecDoub P(311);

for(Doub i=0, j=4.13; i<311; i++, j+=0.05) 
{ Rho[i] = j;
  P[i] = logPres(Rho[i],pp,ee);
  cout<<Rho[i]<<" "<<P[i]<<endl;
}

ofstream file1("./alf1/EoS-alf1.dat");
for(int i = 0; i < 311; i++)
    file1 << Rho[i] << " " << P[i]<<endl;
file1.close();


    
     

    
    
return 0;
}
