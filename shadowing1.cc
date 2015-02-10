#include "iRRAM.h"
#include <vector>
#include <algorithm>

#include <boost/multiprecision/mpfr.hpp>  // Defines the Backend type that wraps MPFR
#include <boost/multiprecision/cpp_bin_float.hpp> 

namespace mp = boost::multiprecision;    
const unsigned int cray_double_precision = 96;
const unsigned int cray_float_precision = 49;
typedef mp::number<mp::backends::cpp_bin_float<cray_float_precision, mp::backends::digit_base_2, void, boost::int16_t, -8192, 8191>, mp::et_off> cray_float;
//typedef double cray_float;
typedef mp::number<mp::cpp_bin_float<cray_double_precision, mp::digit_base_2>>  cray_double;


const double a=3.8;
const cray_double ad=3.8;
const cray_double p0d=0.4;
const double p0=0.4;
using namespace iRRAM;
using std::endl;
using std::vector;
using std::max;

//  computes the logistic map for a dyadic number x
cray_float logistic_map(const cray_float x){
	cray_float ac = static_cast<cray_float>(ad);
	return ac*x*(1-x);
}

//  computes the logistic map for a dyadic number x
REAL logistic_map(const REAL& x){
	return REAL(a)*x*(REAL(1)-x);
}

// computes one solution of the inverse logistic map of the REAL x
// Which solution is chosen is decided by the parameter left. 
// If left is through the solution left of 0.5 is chosen.
REAL inverse_logistic_map(const REAL& x, bool left){
	if(left){
		return REAL(0.5)-sqrt(REAL(0.25)-x*1/REAL(a));
	}
	else{
		return REAL(0.5)+sqrt(REAL(0.25)-x*1/REAL(a));
	}
}

// computes a pseudo orbit (noisy orbit).
// Applies the logistic map N times. 
// The output is a vector of the result after each step
vector<cray_float> pseudo_orbit(const int N){
	vector<cray_float> ans(N+1);
	ans[0] = static_cast<cray_float>(p0d);
	for(int i=1; i<=N; i++){
		ans[i] = logistic_map(ans[i-1]);
	}
	return ans;
} 


void compute(){
	const int N=10000000;
	vector<cray_float> orbit = pseudo_orbit(N);
	// distance between pseudo_orbit and exact orbit
	REAL pi = p0; 
	REAL max_dist = 0;
	for(int i=1; i<=N;i++){
		REAL pi = logistic_map(REAL(p0));
		REAL dist =abs(pi-REAL(float(orbit[i])));
		if(positive(dist-max_dist,-100)){
			max_dist = dist;
		}
	}
	cout << "max distance between real orbit ans pseudo orbit : "<<max_dist<<endl;
	// max distance between pseudo orbit and closest of the real orbits currently computed
	max_dist = 0;
	REAL p = double(orbit[N]);
	//std::cout <<orbit[N]<< endl;
	cout <<p <<endl;
	for(int i=N-1;i>=0; i--){
		bool left= (double(orbit[i]) <= 0.5);
		if(double(orbit[i]) == 0.5)
			std::cout << i<<" "<<orbit[i]<<endl;
		p = inverse_logistic_map(p, left);
		REAL dist = abs(p-REAL(double(orbit[i])));
		if(positive(dist-max_dist,-10))
			max_dist = dist;
	}
	cout << max_dist <<endl;

	/*
	best_p=orbit[0];
	grid_size = power(REAL(2), -20);
	left = REAL(orbit[0])-grid_size;
	n=3;
	while(!bound(max_dist, -11)){
		cout << n << endl;
		for(int i=1; i<n; i++){
			REAL start_point_dist = 0.0;
			REAL p0=left+i*grid_size;
			REAL p=p0;
			for(int j=0;j<=N;j++){
				REAL dist = abs(p-REAL(orbit[j]));
				if(positive(dist-start_point_dist,-100))
					start_point_dist = dist;
			//cout << p << " "<<orbit[j]<<endl;
				p = logistic_map(p);
				if(positive(dist-power(REAL(2),-10),-100)){
					//cout<<"upto"<<j<<endl;
					break;
				}

			}
			if(positive(max_dist-start_point_dist,-100)){
				best_p = p0; // this is now p0
				max_dist = start_point_dist;
			}
			//return 0;
		}
		cout<<best_p<<" "<<orbit[0]<<endl;
		grid_size /= 2;
		n *= 2;
		rshow(max_dist, 500);
		cout << endl;
	}*/
}