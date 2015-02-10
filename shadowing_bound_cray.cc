/***
** Simulate Computations of Cray X-MP
** 
**/
#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/multiprecision/mpfr.hpp>  // Defines the Backend type that wraps MPFR
 #include <boost/multiprecision/cpp_bin_float.hpp> 
#include "fixed_precision_interval.h"

namespace mp = boost::multiprecision;    
const unsigned int cray_double_precision = 96;
const unsigned int cray_float_precision = 49;
typedef mp::number<mp::cpp_bin_float<cray_float_precision, mp::digit_base_2>>  cray_float;
typedef mp::number<mp::cpp_bin_float<cray_double_precision, mp::digit_base_2>>  cray_double;


using std::cout;
using std::endl;
using std::vector; 
using cray_interval = fixed_precision_interval<cray_double_precision>;

const cray_double thickening_constant = pow(cray_double(2.0), -90);
const cray_double error_constant = pow(cray_double(10.0), -25);

const cray_double a=3.8;
const cray_double p0=0.4;
const cray_float a_float=static_cast<cray_float>(a);
const cray_float p0_float=static_cast<cray_float>(p0);

// thicken the given interval
// i.e. enlarge both sides of the interval by the 
// dyadic number given in the second parameter
cray_interval thicken(const cray_interval& x, const cray_double& d){
	return cray_interval(x.left-d, x.right+d);
}

// computes the logistic map f(x) = a*x(1-x).
// Given the dyadic interval x outputs a cray_interval I  
// so that f(t) is in I for all t in x
cray_interval logistic_map(const cray_interval& x){
	cray_interval x_inv; // 1-x
	x_inv.left=1-x.right;
	x_inv.right = 1-x.left;
	return a*x*x_inv;
}

// computes one solution of the inverse logistic map of the cray_interval x
// Which solution is chosen is decided by the parameter left. 
// If left is through the solution left of 0.5 is chosen.
cray_interval inverse_logistic_map(const cray_interval& x, bool left){
	if(left){
		return cray_double(0.5)-sqrt(cray_double(0.25)-x*(cray_double(1)/a));
	}
	else{
		return cray_double(0.5)+sqrt(cray_double(0.25)-x*(cray_double(1)/a));
	}
}


//  computes the logistic map for a dyadic number x
cray_float logistic_map(const cray_float& x){
	return a_float*x*(1-x);
}


// compute the previous interval to a given interval
// using the algorithm from the paper
cray_interval prev_interval(cray_interval x, bool left){

	cray_interval x_thick = thicken(x, thickening_constant);
	cray_interval x_prev = inverse_logistic_map(x_thick, left);
	while(!logistic_map(x_prev).contains(x_thick)){
		x_prev = thicken(x_prev, thickening_constant);
	}
	x_prev = thicken(x_prev, error_constant);
	return x_prev;
}

// computes a pseudo orbit (noisy orbit).
// Applies the logistic map N times. 
// The output is a vector of the result after each step
vector<cray_float> pseudo_orbit(const int N){
	vector<cray_float> ans(N+1);
	ans[0] = p0_float;
	for(int i=1; i<=N; i++){
		ans[i] = logistic_map(ans[i-1]);
	}
	return ans;
} 

// Computes an upper bound for the distance of the 
// shadow orbit x_n to the noisy orbit p_n for N iterates
cray_double shadowing_bound(const int N){
	vector<cray_float> orbit = pseudo_orbit(N);
	cray_interval x = cray_interval(static_cast<cray_double>(orbit[N]), static_cast<cray_double>(orbit[N]));
	cray_double beta = 0; // max distance between the p_n and the shadow orbit
	for(int i=N-1; i>=0; i--){
		bool left= (orbit[i] < cray_double(0.5));
		x = prev_interval(x, left);
		beta = std::max(beta, x.max_dist(static_cast<cray_double>(orbit[i])));
	}
	return beta;
}


int main (int argc,char **argv)
{
	cray_double x = shadowing_bound(10000000);
	std::cout << x <<"\n";

}

