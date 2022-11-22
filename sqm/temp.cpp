#include<C:\Users\windows7\projects\NS_Integrator\nr3.h>

int main()
{    cout<< std::scientific;
    Doub c = 2.998e8;
    Doub c2 = pow(c,2);
    
    
    VecDoub pp(176);
     VecDoub ee(176);
     ifstream infile("./sqm/sqm1.dat");
    Doub baryondensity, pressure, energydensity;
    cout.precision(14);
    for (Doub i = 0; i < 176; i++)
    {
       infile >> baryondensity>>pressure>>energydensity;
       pp[i] = log10(pressure);
       ee[i] = log10(energydensity);
       

       
    }

    infile.close();
    Doub y1 = pp[21];
    Doub y2 = pp[22];
    Doub x1 = ee[21];
    Doub x2 = ee[22];
    cout<<pow(10.0,y1)<<" "<<pow(10.0,x1)<<" "<<pow(10.0,y2)<<" "<<pow(10.0,x2)<<endl;

    Doub h = (x2 - x1)/3.0;
    Doub x,y;
    y = y1;
    x = x1 + h;
    
    while (y <= y2)
    {
        y = y1 + ((y2 - y1)/(x2 - x1))*(x - x1);

        cout << "0.0"<<" "<< pow(10.0,y) << " " << pow(10.0,x) << endl;

        x = x + h;
    }


}