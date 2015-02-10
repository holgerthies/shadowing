/***
** Simulate Computations of Cray X-MP
** 
**/
#include "iRRAM.h"
#include <iostream>
#include <vector>
#include <algorithm>

using std::cout;
using std::endl;
using std::vector; 
using namespace iRRAM;
const REAL thickening_constant = power(REAL(2.0), -90);
const REAL error_constant = power(REAL(10.0), -25);

const double a=3.8;
const double p0=0.4;


INTERVAL sqrt(const INTERVAL& x){
	return INTERVAL(sqrt(x.low), sqrt(x.upp));
}

// thicken the given interval
// i.e. enlarge both sides of the interval by the 
// dyadic number given in the second parameter
 INTERVAL thicken(const INTERVAL& x, const REAL& d){
	return INTERVAL(x.low-d, x.upp+d);
}

//  computes the logistic map for a dyadic number x
REAL logistic_map(const REAL& x){
	return REAL(a)*x*(REAL(1)-x);
}


//  computes the logistic map for a dyadic number x
INTERVAL logistic_map(const INTERVAL& x){
	return REAL(a)*x*(REAL(1)-x);
}

bool mv_subset(const INTERVAL&x, const INTERVAL& y, int prec){
	return positive(y.low-x.low,prec) && positive(x.upp-y.upp, prec); 
}
// computes one solution of the inverse logistic map of the cray_interval x
// Which solution is chosen is decided by the parameter left. 
// If left is through the solution left of 0.5 is chosen.
INTERVAL inverse_logistic_map(const INTERVAL& x, bool left){
	if(left){
		return REAL(0.5)-sqrt(REAL(0.25)-x*(REAL(1)/REAL(a)));
	}
	else{
		return REAL(0.5)+sqrt(REAL(0.25)-x*(REAL(1)/REAL(a)));
	}
}


//  computes the logistic map for a dyadic number x
double logistic_map(const double x){
	return a*x*(1-x);
}


// compute the previous interval to a given interval
// using the algorithm from the paper
INTERVAL prev_interval(INTERVAL& x, bool left){

	INTERVAL x_thick = thicken(x, thickening_constant);
	INTERVAL x_prev = inverse_logistic_map(x_thick, left);
	while(!mv_subset(logistic_map(x_prev),x_thick,-10)){
		x_prev = thicken(x_prev, thickening_constant);
	}
	x_prev = thicken(x_prev, error_constant);
	return x_prev;
}

// mutlivalued max
REAL max(const REAL& a, const REAL& b, const int prec){
	if(positive(a-b,prec))
		return a;
	return b;
}

// mutlivalued max
REAL max(const REAL& a, const REAL& b, const REAL& c, const int prec){
	return max(max(a,b, prec), c, prec);
}

// computes a pseudo orbit (noisy orbit).
// Applies the logistic map N times. 
// The output is a vector of the result after each step
vector<double> pseudo_orbit(const int N){
	vector<double> ans(N+1);
	ans[0] = p0;
	for(int i=1; i<=N; i++){
		ans[i] = logistic_map(ans[i-1]);
	}
	return ans;
} 

// Computes an upper bound for the distance of the 
// shadow orbit x_n to the noisy orbit p_n for N iterates
REAL shadowing_bound(const int N){
	vector<double> orbit = pseudo_orbit(N);
	INTERVAL x = INTERVAL(static_cast<REAL>(orbit[N]), static_cast<REAL>(orbit[N]));
	REAL beta = 0; // max distance between the p_n and the shadow orbit
	for(int i=N-1; i>=0; i--){
		if(i % 10000 == 0)
			std::cout<<i<<endl;
		bool left = (orbit[i] < 0.5);
		x = prev_interval(x, left);

		beta = max(beta, abs(x.low-REAL(orbit[i])), abs(x.upp-REAL(orbit[i])), -10);
	}
	return beta;
}


void compute ()
{
	REAL x = shadowing_bound(10000000);
	rwrite(x, 30);
	iRRAM::cout <<"\n";

}

