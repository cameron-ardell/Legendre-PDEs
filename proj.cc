#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
#include "stepperdopr5.h"
#include "roots.h"

using namespace std;


struct Legendre{
	Doub lambda;

	Legendre(Doub lambda_in) {lambda = lambda_in;}
	~Legendre() {;}


	void operator()(const Doub x, VecDoub_I &y, VecDoub_O &dydx){
		const Doub yy = y[0];
		const Doub kappa = y[1];
		//y[2] is lambda, but don't need to save it since it's been 
		//constant the entire time

		//at initial boundary
		if(x == 1.0){
			dydx[0] = 0.5 * lambda * yy;
			dydx[1] = (lambda * 0.5 - 1.0) * (lambda * 0.25);

		}else{
			dydx[0] = y[1];

			Doub denom = 1.0/(1.0 - x*x);

			dydx[1] = 2.0 * x * denom * kappa - lambda * denom * yy;
		}
		dydx[2] = 0.0;
	}
};

struct Shoot{
	bool even;
	string fileName;

	Shoot(Doub l) {
		//checks if even or odd function
		if( fmod(l, 2.0) == 0.0){
			even = true;
		}else{
			even = false;
		}

		//which file to write data to
		if(l == 0.0){
			fileName = "l=0.dat";
		}else if(l == 1.0){
			fileName = "l=1.dat";
		}else if(l == 2.0){
			fileName = "l=2.dat";
		}else if(l == 3.0){
			fileName = "l=3.dat";
		}else if(l == 4.0){
			fileName = "l=4.dat";
		}
	}

	~Shoot() { cout << "destructing Shoot..." << endl;}

	//evaluates function for specific lambda value as an ODE
	Doub operator()(Doub lambdaa){

		ofstream outfile;
		outfile.open(fileName.c_str());
		outfile.setf(ios::left);
		outfile << setw(16) << "# x " << setw(16) << " y " << endl;
		outfile << "#====================================================" << endl;

		//vector size for holding three variables
	    int nvar = 3;
	    //sets absolute and relative tolerances
	    //also sets minimum step size and first guess of step size
	    const Doub atol = 1.e-5, rtol = atol, h1=0.001, hmin=0.0;
	    
	    //amount of times interval should happen
	    const Doub t_in = 1.0;
	    const Doub t_final = 0.0;
    
        //total amount of steps
	    const int N_steps = 1000;
	    
	    //how much t should be incremented by each iteration
	    //is dependent on total amount of steps
	    const Doub delta_t = (t_final - t_in)/N_steps;

		//vector that holds the starting y values (will be altered as code runs)
	    VecDoub ystart(nvar);
	    
	    //need output object to control saving of intermediate values
	    Output out;

	    //initial conditions, can be reset by user
	    ystart[0] = 1.0;
	    ystart[1] = lambdaa * 0.5;
	    ystart[2] = lambdaa;

		Legendre leg(lambdaa);

		for (Doub t = t_in; t > t_final; t+= delta_t) {
	        //do Odeint, this will automatically update ystart
	        Odeint<StepperDopr5<Legendre> > ode(ystart,t,t+delta_t,atol,
	                                            rtol,h1,hmin,out,leg);
	        ode.integrate();

	        outfile << setw(16) << t << setw(16) << ystart[0] << endl;
	        if(even == true){
	        	outfile << setw(16) << -t << setw(16) << ystart[0] << endl;
	        }else{
	        	outfile << setw(16) << -t << setw(16) << -ystart[0] << endl;
	        }
    	}

    	if(even == true){
    		return ystart[1];
    	}else{
    		return ystart[0];
    	}

	}
	outfile.close();
};


int main() {

	Doub acc = 1e-16;
	Doub brack1, brack2;

	for(Doub l = 0.0; l < 5.0; l++){
		brack1 = l*(l + 1.0) - 0.5;
		brack2 = l*(l + 1.0) + 0.5;

		Shoot shot(l);

		Doub root = zbrent(shot, brack1, brack2, acc);

		cout << "l = " << l << "\nlambda should be here: " << root << endl;
	}

	return 0;
}