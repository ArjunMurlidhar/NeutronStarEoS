struct Base_interp
//Abstract base class used by all interpolation routines in this chapter. Only the routine interp
//is called directly by the user.
{
Int n, mm, jsav, cor, dj;
const Doub *xx, *yy;
Base_interp(VecDoub_I &x, const Doub *y, Int m): n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) {
dj = MIN(1,(int)pow((Doub)n,0.25));
}              //Constructor: Set up for interpolating on a table of x’s and y’s of length m. Normally called
//by a derived class, not by the user.
Doub interp(Doub x) {
//Given a value x, return an interpolated value, using data pointed to by xx and yy.
Int jlo = cor ? hunt(x) : locate(x);
return rawinterp(jlo,x);
}
Int locate(const Doub x); //See definitions below.
Int hunt(const Doub x);
Doub virtual rawinterp(Int jlo, Doub x) = 0;
//Derived classes provide this as the actual interpolation method.
};

Int Base_interp::locate(const Doub x) 
//Given a value x, return a value j such that x is (insofar as possible) centered in the subrange
//xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either
//increasing or decreasing. The returned value is not less than 0, nor greater than n-1.
{
Int ju,jm,jl;
if (n < 2 || mm < 2 || mm > n) throw("locate size error");
Bool ascnd=(xx[n-1] >= xx[0]); //True if ascending order of table, false otherwise.
jl=0;                          //Initialize lower
ju=n-1;                        //and upper limits.
while (ju-jl > 1) {            //If we are not yet done,
jm = (ju+jl) >> 1;             //compute a midpoint,
if (x >= xx[jm] == ascnd)
jl=jm;                         //and replace either the lower limit
else
ju=jm;                         //or the upper limit, as appropriate.
}                              //Repeat until the test condition is satisfied.
cor = abs(jl-jsav) > dj ? 0 : 1; //Decide whether to use hunt or locate next time.
jsav = jl;
return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}

Int Base_interp::hunt(const Doub x)
//Given a value x, return a value j such that x is (insofar as possible) centered in the subrange
//xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either
//increasing or decreasing. The returned value is not less than 0, nor greater than n-1.
{
Int jl=jsav, jm, ju, inc=1;
if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
Bool ascnd=(xx[n-1] >= xx[0]); //True if ascending order of table, false otherwise.
if (jl < 0 || jl > n-1) {     //Input guess not useful. Go immediately to bisection.
jl=0; 
ju=n-1;
} else {
if (x >= xx[jl] == ascnd) {             //Hunt up:
for (;;) {
ju = jl + inc;
if (ju >= n-1) { ju = n-1; break;}      //Off end of table.
else if (x < xx[ju] == ascnd) break;    //Found bracket.
else {                                  //Not done, so double the increment and try again.
jl = ju;
inc += inc;
}
}
} else {                                //Hunt down:
ju = jl;
for (;;) {
jl = jl - inc;
if (jl <= 0) { jl = 0; break;}          //Off end of table.
else if (x >= xx[jl] == ascnd) break;   //Found bracket.
else {                                  //Not done, so double the increment and try again.
ju = jl;
inc += inc;
}
}
}
}
while (ju-jl > 1) {                  //Hunt is done, so begin the final bisection phase:
jm = (ju+jl) >> 1;
if (x >= xx[jm] == ascnd)
jl=jm;
else
ju=jm;
}
cor = abs(jl-jsav) > dj ? 0 : 1;     //Decide whether to use hunt or locate next time.
jsav = jl; 
return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}

struct Spline_interp : Base_interp
//Cubic spline interpolation object. Construct with x and y vectors, and (optionally) values of
//the first derivative at the endpoints, then call interp for interpolated values.
{
VecDoub y2;
Spline_interp(VecDoub_I &xv, VecDoub_I &yv, Doub yp1=1.e99, Doub ypn=1.e99)
: Base_interp(xv,&yv[0],2), y2(xv.size())
{sety2(&xv[0],&yv[0],yp1,ypn);}
Spline_interp(VecDoub_I &xv, const Doub *yv, Doub yp1=1.e99, Doub ypn=1.e99)
: Base_interp(xv,yv,2), y2(xv.size())
{sety2(&xv[0],yv,yp1,ypn);}
void sety2(const Doub *xv, const Doub *yv, Doub yp1, Doub ypn);
Doub rawinterp(Int jl, Doub xv);
};

void Spline_interp::sety2(const Doub *xv, const Doub *yv, Doub yp1, Doub ypn) 
/*This routine stores an array y2[0..n-1] with second derivatives of the interpolating function
at the tabulated points pointed to by xv, using function values pointed to by yv. If yp1 and/or
ypn are equal to 1  1099 or larger, the routine is signaled to set the corresponding boundary
condition for a natural spline, with zero second derivative on that boundary; otherwise, they are
the values of the first derivatives at the endpoints.*/
{
Int i,k;
Doub p,qn,sig,un;
Int n=y2.size();
VecDoub u(n-1);
if (yp1 > 0.99e99) //The lower boundary condition is set either to be “natural”
 y2[0]=u[0]=0.0; 
else {             //or else to have a specified first derivative.
 y2[0] = -0.5;
 u[0]=(3.0/(xv[1]-xv[0]))*((yv[1]-yv[0])/(xv[1]-xv[0])-yp1);
}
for (i=1;i<n-1;i++) {   //This is the decomposition loop of the tridiagonal algorithm.
                        //y2 and u are used for temporary storage of the decomposed factors.
 sig=(xv[i]-xv[i-1])/(xv[i+1]-xv[i-1]);
 p=sig*y2[i-1]+2.0;
 y2[i]=(sig-1.0)/p;
 u[i]=(yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1]);
 u[i]=(6.0*u[i]/(xv[i+1]-xv[i-1])-sig*u[i-1])/p;
}
if (ypn > 0.99e99) //The upper boundary condition is set either to be “natural”
 qn=un=0.0; 
else {             //or else to have a specified first derivative.
 qn=0.5;
 un=(3.0/(xv[n-1]-xv[n-2]))*(ypn-(yv[n-1]-yv[n-2])/(xv[n-1]-xv[n-2]));
}
y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
for (k=n-2;k>=0;k--) //This is the backsubstitution loop of the tridiagonal algorithm.
y2[k]=y2[k]*y2[k+1]+u[k]; 
}

Doub Spline_interp::rawinterp(Int jl, Doub x) 
/*Given a value x, and using pointers to data xx and yy, and the stored vector of second derivatives
y2, this routine returns the cubic spline interpolated value y.*/
{
Int klo=jl,khi=jl+1;
Doub y,h,b,a;
h=xx[khi]-xx[klo];
if (h == 0.0) throw("Bad input to routine splint"); //The xa’s must be distinct.
a=(xx[khi]-x)/h; 
b=(x-xx[klo])/h;                                  //Cubic spline polynomial is now evaluated.
y=a*yy[klo]+b*yy[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0;
return y;
}

