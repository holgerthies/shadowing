// Class template for interval arithmetic.
// The interval endpoints are saved as fixed precision floating point number
// with precision given by the template parameter. 
#ifndef _FIXED_PRECISION_INTERVAL_H
#define _FIXED_PRECISION_INTERVAL_H
#include <boost/multiprecision/mpfr.hpp>  // Defines the Backend type that wraps MPFR 
#include <boost/multiprecision/cpp_bin_float.hpp> 
namespace mp = boost::multiprecision;  
template <unsigned int precision>
class fixed_precision_interval { 
	template <unsigned int p>
	using prec_float = typename mp::number<mp::cpp_bin_float<p, mp::digit_base_2>>  ;
	public:
		prec_float<precision> left;
		prec_float<precision> right;
		fixed_precision_interval() {};
		fixed_precision_interval(const prec_float<precision>& left, const prec_float<precision>& right) : left(left), right(right) {};
		// computes the sum of two intervals
		fixed_precision_interval<precision> operator+(const fixed_precision_interval& x) const;
		// computes the sum of interval and float
		fixed_precision_interval<precision> operator+(const prec_float<precision>& y) const;
		// right-side scalar plus
		template <unsigned int p>
		friend fixed_precision_interval<p> operator+(const prec_float<p>& x, const fixed_precision_interval<p>& y);
		// computes the inverse
		fixed_precision_interval<precision> operator-() const;
		// computes the difference of two intervals		
		fixed_precision_interval<precision> operator-(const fixed_precision_interval<precision>& x) const;
		// computes the difference of interval and float		
		fixed_precision_interval<precision> operator-(const prec_float<precision>& x) const;
		// right hand side scalar minus		
		template <unsigned int p>
		friend fixed_precision_interval<p> operator-(const prec_float<p>& x, const fixed_precision_interval<p>& y);
		// computes the product of two intervals
		fixed_precision_interval<precision> operator*(const fixed_precision_interval<precision>& x) const;
		// scalar product with an interval
		fixed_precision_interval<precision> operator*(const prec_float<precision>& x) const;
		// right hand side scalar product
		template <unsigned int p>
		friend fixed_precision_interval<p> operator*(const prec_float<p>& x, const fixed_precision_interval<p>& y);
		// computes the sqrt of an interval.
		template <unsigned int p>
		friend fixed_precision_interval<p> sqrt(const fixed_precision_interval<p>& x);
		// checks if the interval is contains the other
		bool contains(const fixed_precision_interval<precision>& x) const;
		// checks if point is contained in the interval
		bool contains(const prec_float<precision>& x) const;
		// returns interval length
		prec_float<precision> size() const;
		// distance of point from interval
		prec_float<precision> dist(const prec_float<precision>& x) const;
		// maximum distance of point from any point inside the interval
		prec_float<precision> max_dist(const prec_float<precision>& x) const;
		
};

// Template definitions

// [a,b]+[c,d] = [a+c, b+d]
template <unsigned int precision>
fixed_precision_interval<precision> fixed_precision_interval<precision>::operator+(const fixed_precision_interval<precision>& x) const{
	fixed_precision_interval<precision> ans;
	ans.left = this->left+x.left;
	ans.right = this->right+x.right;
	return ans;
}


// [a,b]+d = [a+d, b+d]
template <unsigned int precision>
fixed_precision_interval<precision> fixed_precision_interval<precision>::operator+(const boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>>& x) const{
	fixed_precision_interval<precision> ans;
	ans.left = this->left+x;
	ans.right = this->right+x;
	return ans;
}


template<unsigned int precision>
fixed_precision_interval<precision> operator+(const boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>>& x, const fixed_precision_interval<precision>& y){
	return y+x;
}


// [a,b]-[c,d] = [a-d, b-c]
template<unsigned int precision>
fixed_precision_interval<precision> fixed_precision_interval<precision>::operator-(const fixed_precision_interval<precision>& x) const{
	fixed_precision_interval<precision> ans;
	ans.left = this->left-x.right;
	ans.right = this->right-x.left;
	return ans;
}


// -[a,b] = [-b, -a]
template<unsigned int precision>
fixed_precision_interval<precision> fixed_precision_interval<precision>::operator-() const{
	fixed_precision_interval<precision> ans;
	ans.left = -this->right;
	ans.right = -this->left;
	return ans;
}

template<unsigned int precision>
fixed_precision_interval<precision> operator-(const boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>>& x, const fixed_precision_interval<precision>& y){
	return -y+x;
}


// [a,b]-d = [a-d, b-d]
template<unsigned int precision>
fixed_precision_interval<precision> fixed_precision_interval<precision>::operator-(const boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>>& x) const{
	fixed_precision_interval<precision> ans;
	ans.left = this->left-x;
	ans.right = this->right-x;
	return ans;
}


// [a,b]*[c,d] = [min(ac, ad, bc, bd), max(ac, ad, bc, bd)]
template<unsigned int precision>
fixed_precision_interval<precision> fixed_precision_interval<precision>::operator*(const fixed_precision_interval<precision>& x) const{
	fixed_precision_interval<precision> ans;
	// pairwise product between all interval endpoints
	std::vector<boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>>> pw_products = {this->left*x.left, this->left*x.right, this->right*x.left, this->right*x.right}; 
	ans.left = *min_element(pw_products.begin(), pw_products.end());
	ans.right = *max_element(pw_products.begin(), pw_products.end());
	return ans;
}

// [a,b]*d = [min(ad, bd), max(ad, bd)]
template<unsigned int precision>
fixed_precision_interval<precision> fixed_precision_interval<precision>::operator*(const boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>>& d) const{
	fixed_precision_interval<precision> ans;
	ans.left = std::min(d*this->left, d*this->right);
	ans.right = std::max(d*this->left, d*this->right);
	return ans;
}


template<unsigned int precision>
fixed_precision_interval<precision> operator*(const boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>>& x, const fixed_precision_interval<precision>& y){
	return y*x;
}


template<unsigned int precision>
bool fixed_precision_interval<precision>::contains(const fixed_precision_interval<precision>& x) const{
	return this->left <= x.left && this->right >= x.right;
}


template<unsigned int precision>
bool fixed_precision_interval<precision>::contains(const boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>>& x) const{
	return this->left <= x && this->right >= x;
}


template<unsigned int precision>
boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>> fixed_precision_interval<precision>::size() const{
	return this->right - this->left;
}


template<unsigned int precision>
boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>> fixed_precision_interval<precision>::dist(const boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>>& x) const{
	if(this->contains(x)) return 0;
	return std::min(abs(this->left-x), abs(this->right-x));
}


template<unsigned int precision>
boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>> fixed_precision_interval<precision>::max_dist(const boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>>& x) const{
	return std::max(abs(this->left-x), abs(this->right-x));
}


// Uses the boost multiprecision sqrt function to compute the end points
template<unsigned int precision>
fixed_precision_interval<precision> sqrt(const fixed_precision_interval<precision>& x){
	fixed_precision_interval<precision> ans;
	boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>>  p1 = boost::multiprecision::sqrt(x.left);
	boost::multiprecision::number<boost::multiprecision::cpp_bin_float<precision, boost::multiprecision::digit_base_2>>  p2 = boost::multiprecision::sqrt(x.right);
	ans.left = std::min(p1, p2);
	ans.right = std::max(p1,p2);
	return ans;
}
#endif