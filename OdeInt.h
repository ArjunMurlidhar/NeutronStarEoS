struct Output {
//Structure for output from ODE solver such as odeint.
Int kmax; //Current capacity of storage arrays.
Int nvar;
Doub nspace; //spacing of intervals for dense output.
bool dense; //true if dense output requested.
Int count; //Number of values actually saved.
Doub x1,x2,xout,dxout;
VecDoub xsave; //Results stored in the vector xsave[0..count-1] and the
MatDoub ysave; //matrix ysave[0..nvar-1][0..count-1].
Output() : kmax(-1),dense(false),count(0) {} //Default constructor gives no output.
Output(const Doub nspacee) : kmax(500),nspace(nspacee),count(0),xsave(kmax) {
/*Constructor provides dense output at equally spaced intervals of spacing nspace. If nsave <= 0, output
is saved only at the actual integration steps.*/
dense = nspace > 0 ? true : false;
}
void init(const Int neqn, const Doub xlo) {
//Called by Odeint constructor, which passes neqn, the number of equations, xlo, the starting
//point of the integration, and xhi, the ending point.
nvar=neqn;
if (kmax == -1) return;
ysave.resize(nvar,kmax);
if (dense) {
x1=xlo;
xout=x1;
dxout=nspace;
}
}
void resize() {
//Resize storage arrays by a factor of two, keeping saved data.
Int kold=kmax;
kmax *= 2;
VecDoub tempvec(xsave);
xsave.resize(kmax);
for (Int k=0; k<kold; k++)
xsave[k]=tempvec[k];
MatDoub tempmat(ysave);
ysave.resize(nvar,kmax);
for (Int i=0; i<nvar; i++)
for (Int k=0; k<kold; k++)
ysave[i][k]=tempmat[i][k];
}

template <class Stepper>
void save_dense(Stepper &s, const Doub xout, const Doub h) {
/*Invokes dense_out function of stepper routine to produce output at xout. Normally called
by out rather than directly. Assumes that xout is between xold and xold+h, where the
stepper must keep track of xold, the location of the previous step, and x=xold+h, the
current step.*/
if (count == kmax) resize();
for (Int i=0;i<nvar;i++)
ysave[i][count]=s.dense_out(i,xout,h);
xsave[count++]=xout;
}
void save(const Doub x, VecDoub_I &y) {
//Saves values of current x and y.
if (kmax <= 0) return;
if (count == kmax) resize();
for (Int i=0;i<nvar;i++)
ysave[i][count]=y[i];
xsave[count++]=x;
}
template <class Stepper>
void out(const Int nstp,const Doub x,VecDoub_I &y,Stepper &s,const Doub h) {
/*Typically called by Odeint to produce dense output. Input variables are nstp, the current
step number, the current values of x and y, the stepper s, and the stepsize h. A call with
nstp=-1 saves the initial values. The routine checks whether x is greater than the desired
output point xout. If so, it calls save_dense.*/
if (!dense)
throw("dense output not set in Output!");
if (nstp == -1) {
save(x,y);
xout += dxout;
} else {
while ((x-xout) > 0.0) {
save_dense(s,xout,h);
xout += dxout;
}
}
}
};

template<class Stepper>
struct Odeint {
//Driver for ODE solvers with adaptive stepsize control. The template parameter should be one
//of the derived classes of StepperBase defining a particular integration algorithm.
static const Int MAXSTP=50000; //Take at most MAXSTP steps.
Doub EPS;
Int nok;
Int nbad;
Int nvar;
Doub x1,x2,hmin;
bool dense;                 //true if dense output requested by out.
VecDoub y,dydx; 
VecDoub &ystart;
Output &out;
typename Stepper::Dtype &derivs; //Get the type of derivs from the stepper.
Stepper s; 
Int nstp;
Doub x,h;
Odeint(VecDoub_IO &ystartt,const Doub xx1,
const Doub atol,const Doub rtol,const Doub h1,
const Doub hminn,Output &outt,typename Stepper::Dtype &derivss);
void integrate(); //Does the actual integration.
};

template<class Stepper>
Odeint<Stepper>::Odeint(VecDoub_IO &ystartt, const Doub xx1, 
const Doub atol, const Doub rtol, const Doub h1, const Doub hminn, Output &outt,typename Stepper::Dtype &derivss) : nvar(ystartt.size()),
y(nvar),dydx(nvar),ystart(ystartt),x(xx1),nok(0),nbad(0),
x1(xx1),hmin(hminn),dense(outt.dense),out(outt),derivs(derivss),
s(y,dydx,x,atol,rtol,dense) {
EPS=numeric_limits<Doub>::epsilon();
h=h1;
for (Int i=0;i<nvar;i++) y[i]=ystart[i];
out.init(s.neqn,x1);
}

template<class Stepper>
void Odeint<Stepper>::integrate() {
   //cout << "Start integrate()" <<endl;
derivs(x,y,dydx);
//cout<<x<<" "<<y[0]<<" "<<y[1]<<" "<<dydx[0]<<" "<<dydx[1]<<endl;
if (dense)       //Store initial values.
out.out(-1,x,y,s,h);
else
out.save(x,y);
for (nstp=0;nstp<MAXSTP;nstp++) {
    //if(nstp == 85) break;
    //cout << "Inside for loop, step: " << nstp << endl ;
s.step(h,derivs); //Take a step.
//cout<< "took step" << nstp <<endl;
if (s.hdid == h) ++nok; else ++nbad;
if (dense)
out.out(nstp,x,y,s,s.hdid);
else
out.save(x,y);
if (y[0] >= -s.pbound && y[0] <= s.pbound) {                   //Are we done?
for (Int i=0;i<nvar;i++) ystart[i]=y[i];       //Update ystart.
if (out.kmax > 0 && abs(out.xsave[out.count-1]-x2) > 100.0*abs(x2)*EPS)
out.save(x,y);                            //Make sure last step gets saved.
return;                                    //Normal exit.
}
if (abs(s.hnext) <= hmin) throw("Step size too small in Odeint");
h=s.hnext;
//cout << "end of loop" <<endl;
//cout<<x<<" "<<y[0]<<" "<<y[1]<<" "<<dydx[0]<<" "<<dydx[1]<<endl;

}
throw("Too many steps in routine Odeint");
}