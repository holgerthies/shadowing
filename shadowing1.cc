#include "iRRAM.h"
#include <vector>
#include <algorithm>

#include <boost/multiprecision/mpfr.hpp>  // Defines the Backend type that wraps MPFR
#include <boost/multiprecision/cpp_bin_float.hpp> 

namespace mp = boost::multiprecision;    
const unsigned int cray_double_precision = 96;
const unsigned int cray_float_precision = 49;
//typedef mp::number<mp::backends::cpp_bin_float<cray_float_precision, mp::backends::digit_base_2, void, boost::int16_t, -8192, 8191>, mp::et_off> cray_float;
//typedef double cray_float;
typedef mp::number<mp::cpp_bin_float<cray_double_precision, mp::digit_base_2>>  cray_double;
template <unsigned int prec>
using number = mp::number<mp::cpp_bin_float<prec, mp::digit_base_2>>; 
using cray_float = number<cray_float_precision>;
using namespace iRRAM;
using std::endl;
using std::vector;
using std::max;

//  computes the logistic map for a dyadic number x
template <unsigned int prec>
number<prec> logistic_map(const number<prec>& a, const number<prec>& x){
	return a*x*(1-x);
}


//  computes the logistic map for a real number x
REAL logistic_map(const REAL& a, const REAL& x){
	return a*x*(REAL(1)-x);
}

// computes one solution of the inverse logistic map of the REAL x
// Which solution is chosen is decided by the parameter left. 
// If left is through the solution left of 0.5 is chosen.
REAL inverse_logistic_map(const REAL& a, const REAL& x, bool left){
	if(left){
		return REAL(0.5)-sqrt(REAL(0.25)-x*1/a);
	}
	else{
		return REAL(0.5)+sqrt(REAL(0.25)-x*1/a);
	}
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

// shadowing distance
template <unsigned int prec>
REAL shadowing_bound(const int N, const double a, const double p0){
	vector<number<prec>> orbit = pseudo_orbit<prec>(static_cast<number<prec>>(cray_double(a)), static_cast<number<prec>>(cray_double(p0)), N);
	REAL max_dist = 0;
	REAL p = double(orbit[N]);
	for(int i=N-1;i>=0; i--){
		bool left= (double(orbit[i]) <= 0.5);
		p = inverse_logistic_map(REAL(a), p, left);
		REAL dist = abs(p-REAL(double(orbit[i])));
		if(positive(dist-max_dist,-10))
			max_dist = dist;
	}
	return max_dist;
}

void compute(){
	//const double a = 3.8;
	//const double p0 = 0.4;
	vector<double> as = {3.6, 3.635, 3.65, 3.7, 3.75, 3.8, 3.86, 3.91};
	vector<double> p0s = {0.3, 0.4, 0.5,0.6};
	for(auto a : as){
		for(auto p0 : p0s){
			for(int N=10; N <= 100000000; N *= 10){
				REAL d = shadowing_bound<10>(N, a, p0);
				cout << "10 " << a <<" "<<p0<<" "<< N << " ";
				rwrite(d, 25);
				cout << endl;
		 		d = shadowing_bound<20>(N, a, p0);
				cout << "20 " << a <<" "<<p0<<" "<< N << " ";
				rwrite(d, 25);
				cout << endl;
		 		d = shadowing_bound<30>(N, a, p0);
				cout << "30 " << a <<" "<<p0<<" "<< N << " ";
				rwrite(d, 25);
				cout << endl;
		 		d = shadowing_bound<40>(N, a, p0);
				cout << "40 " << a <<" "<<p0<<" "<< N << " ";
				rwrite(d, 25);
				cout << endl;
		 		d = shadowing_bound<50>(N, a, p0);
				cout << "50 " << a <<" "<<p0<<" "<< N << " ";
				rwrite(d, 25);
				cout << endl;
		 		d = shadowing_bound<60>(N, a, p0);
				cout << "60 " << a <<" "<<p0<<" "<< N << " ";
				rwrite(d, 25);
				cout << endl;
		 		d = shadowing_bound<70>(N, a, p0);
				cout << "70 " << a <<" "<<p0<<" "<< N << " ";
				rwrite(d, 25);
				cout << endl;
			}
		}
	}

}