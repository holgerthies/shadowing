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
#include <boost/mpl/for_each.hpp>

namespace mp = boost::multiprecision;    
const unsigned int cray_double_precision = 96;
const unsigned int cray_float_precision = 49;
typedef mp::number<mp::cpp_bin_float<cray_float_precision, mp::digit_base_2>>  cray_float;
typedef mp::number<mp::cpp_bin_float<cray_double_precision, mp::digit_base_2>>  cray_double;

template <unsigned int prec>
using number = mp::number<mp::cpp_bin_float<prec, mp::digit_base_2>>; 

using std::cout;
using std::endl;
using std::vector; 
using cray_interval = fixed_precision_interval<cray_double_precision>;

// thicken the given interval
// i.e. enlarge both sides of the interval by the 
// dyadic number given in the second parameter
template <unsigned int prec>
fixed_precision_interval<prec> thicken(const fixed_precision_interval<prec>& x, const number<prec>& d){
	return fixed_precision_interval<prec>(x.left-d, x.right+d);
}

// computes the logistic map f(x) = a*x(1-x).
// Given the dyadic interval x outputs a cray_interval I  
// so that f(t) is in I for all t in x
template <unsigned int prec>
fixed_precision_interval<prec> logistic_map(const number<prec> a, const fixed_precision_interval<prec>& x){
	fixed_precision_interval<prec> x_inv; // 1-x
	x_inv.left=1-x.right;
	x_inv.right = 1-x.left;
	return a*x*x_inv;
}

// computes one solution of the inverse logistic map of the cray_interval x
// Which solution is chosen is decided by the parameter left. 
// If left is through the solution left of 0.5 is chosen.
template <unsigned int prec>
fixed_precision_interval<prec> inverse_logistic_map(const fixed_precision_interval<prec>& x, const number<prec>& a, bool left){
	if(left){
		return number<prec>(0.5)-sqrt(number<prec>(0.25)-x*(number<prec>(1)/a));
	}
	else{
		return number<prec>(0.5)+sqrt(number<prec>(0.25)-x*(number<prec>(1)/a));
	}
}


//  computes the logistic map for a dyadic number x
template <unsigned int prec>
number<prec> logistic_map(const number<prec>& a, const number<prec>& x){
	return a*x*(1-x);
}


// compute the previous interval to a given interval
// using the algorithm from the paper
template <unsigned int prec>
fixed_precision_interval<prec> prev_interval(const fixed_precision_interval<prec>& x, const number<prec>& a,  const number<prec> thickening_constant, const number<prec> error_constant, bool left){

	fixed_precision_interval<prec> x_thick = thicken(x, thickening_constant);
	fixed_precision_interval<prec> x_prev = inverse_logistic_map(x_thick, a, left);
	while(!logistic_map(a, x_prev).contains(x_thick)){
		x_prev = thicken(x_prev, thickening_constant);
	}
	x_prev = thicken(x_prev, error_constant);
	return x_prev;
}

// computes a pseudo orbit (noisy orbit).
// Applies the logistic map N times. 
// The output is a vector of the result after each step
template <unsigned int prec>
vector<number<prec>> pseudo_orbit(const number<prec>& a, const number<prec>& p0, const int N){
	vector<number<prec>> ans(N+1);
	ans[0] = p0;
	for(int i=1; i<=N; i++){
		ans[i] = logistic_map(a, ans[i-1]);
	}
	return ans;
} 

// Computes an upper bound for the distance of the 
// shadow orbit x_n to the noisy orbit p_n for N iterates
template <unsigned int float_prec, unsigned int double_prec>
number<double_prec> shadowing_bound(const int N, const number<double_prec> a, const number<double_prec> p0, const number<double_prec> thickening_constant, const number<double_prec> error_constant){
	vector<number<float_prec>> orbit = pseudo_orbit<float_prec>(static_cast<number<float_prec>>(a), static_cast<number<float_prec>>(p0), N);
	fixed_precision_interval<double_prec> x = fixed_precision_interval<double_prec>(static_cast<number<double_prec>>(orbit[N]), static_cast<number<double_prec>>(orbit[N]));
	number<double_prec> beta = 0; // max distance between the p_n and the shadow orbit
	for(int i=N-1; i>=0; i--){
		// if(i % 1000000 == 0){ 
		// 	cout<<i<<endl;
		// 	cout<<x.right<<" "<<a/4<<endl;
		// }
		bool left= (orbit[i] < number<double_prec>(0.5));
		x = prev_interval(x, a, thickening_constant, error_constant, left);
		if(x.right > a/4)
			return -1;
		beta = std::max(beta, x.max_dist(static_cast<number<double_prec>>(orbit[i])));
	}
	return beta;
}

// checks if a breakdown occurs when using the procedure to compute the  
// pseudo orbit up to the Nth iterate. A breakdown point is a point that 
// is bigger than a/4, since then the inverse function is not defined
template <unsigned int float_prec, unsigned int double_prec>
bool has_breakdown(const int N, const number<double_prec> a, const number<double_prec> p0, const number<double_prec> thickening_constant, const number<double_prec> error_constant){
	vector<number<float_prec>> orbit = pseudo_orbit<float_prec>(static_cast<number<float_prec>>(a), static_cast<number<float_prec>>(p0), N);
	fixed_precision_interval<double_prec> x = fixed_precision_interval<double_prec>(static_cast<number<double_prec>>(orbit[N]), static_cast<number<double_prec>>(orbit[N]));
	for(int i=N-1; i>=0; i--){
		
		bool left= (orbit[i] < number<double_prec>(0.5));
		x = prev_interval(x, a, thickening_constant, error_constant, left);
		if(i % 10000 == 0){ 
			cout<<i<<endl;
			cout<<x.right<<" "<<a/4<<endl;
		}
		if(x.right > a/4)
			return true;
	}
	return false;
}

// finds the first breakdown point if there is
// one between l and r, otherwise returns 0
// int binary_search_breakdown(const int l, const int r){
// 	cout << l<<" "<<r<<endl;
// 	if(r <= l){
// 		if(has_breakdown(l))
// 			return l;
// 		return 0;
// 	}
// 	int middle = (l+r)/2;
// 	if(has_breakdown(middle))
// 		return binary_search_breakdown(l, middle);
// 	return binary_search_breakdown(middle+1, r);
// }

int main (int argc,char **argv)
{
	namespace mpl = boost::mpl;
	//cout << binary_search_breakdown(1,100000000)<<endl;
	const cray_double thickening_constant = pow(cray_double(2.0), -90);
	const cray_double error_constant = pow(cray_double(10.0), -25);
	//cout << has_breakdown<35, cray_double_precision>(10000000, cray_double(3.8), cray_double(0.4), thickening_constant, error_constant)<<endl;
	vector<cray_double> as = {3.6, 3.635, 3.65, 3.7, 3.75, 3.8, 3.86, 3.91};
	vector<cray_double> p0s = {0.3, 0.4, 0.5,0.6};
	for(auto a : as){
		for(auto p0 : p0s){
			auto x = shadowing_bound<5, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 5 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<10, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 10 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<15, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 15 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<20, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 20 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<25, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 25 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<30, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 30 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<35, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 35 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<40, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 40 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<45, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 45 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<50, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 50 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<55, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 55 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<60, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 60 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<65, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 65 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<70, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 70 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<75, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 75 << " " << a << " " << p0 << " " << x << endl;
			x = shadowing_bound<80, cray_double_precision>(10000000, a, p0, thickening_constant, error_constant);
			std::cout << 80 << " " << a << " " << p0 << " " << x << endl;			
		}
	}

}

