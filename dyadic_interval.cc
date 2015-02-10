#include "dyadic_interval.h"
#include <vector>
#include <algorithm>
using namespace iRRAM;
using std::min_element;
using std::max_element;
using std::vector;
// computes the sum of two intervals
// [a,b]+[c,d] = [a+c, b+d]
dyadic_interval dyadic_interval::operator+(const dyadic_interval& x) const{
	dyadic_interval ans;
	ans.left = this->left+x.left;
	ans.right = this->right+x.right;
	return ans;
}

// computes the sum of interval and dyadic number
// [a,b]+d = [a+d, b+d]
dyadic_interval dyadic_interval::operator+(const DYADIC& x) const{
	dyadic_interval ans;
	ans.left = this->left+x;
	ans.right = this->right+x;
	return ans;
}

// right-side scalar plus
dyadic_interval operator+(const DYADIC& x, const dyadic_interval& y){
	return y+x;
}

// computes the difference between two intervals
// [a,b]-[c,d] = [a-d, b-c]
dyadic_interval dyadic_interval::operator-(const dyadic_interval& x) const{
	dyadic_interval ans;
	ans.left = this->left-x.right;
	ans.right = this->right-x.left;
	return ans;
}

// computes the inverse
// -[a,b] = [-b, -a]
dyadic_interval dyadic_interval::operator-() const{
	dyadic_interval ans;
	ans.left = -this->right;
	ans.right = -this->left;
	return ans;
}


// right-side scalar minus
dyadic_interval operator-(const DYADIC& x, const dyadic_interval& y){
	return -y+x;
}

// computes the difference of interval and dyadic number
// [a,b]-d = [a-d, b-d]
dyadic_interval dyadic_interval::operator-(const DYADIC& x) const{
	dyadic_interval ans;
	ans.left = this->left-x;
	ans.right = this->right-x;
	return ans;
}

// computes the product between two intervals
// [a,b]*[c,d] = [min(ac, ad, bc, bd), max(ac, ad, bc, bd)]
dyadic_interval dyadic_interval::operator*(const dyadic_interval& x) const{
	dyadic_interval ans;
	// pairwise product between all interval endpoints
	vector<DYADIC> pw_products = {this->left*x.left, this->left*x.right, this->right*x.left, this->right*x.right}; 
	ans.left = *min_element(pw_products.begin(), pw_products.end());
	ans.right = *max_element(pw_products.begin(), pw_products.end());
	return ans;
}

// scalar product with an interval
// [a,b]*d = [min(ad, bd), max(ad, bd)]
dyadic_interval dyadic_interval::operator*(const DYADIC& d) const{
	dyadic_interval ans;
	ans.left = std::min(d*this->left, d*this->right);
	ans.right = std::max(d*this->left, d*this->right);
	return ans;
}

// right-side scalar product
dyadic_interval operator*(const DYADIC& x, const dyadic_interval& y){
	return y*x;
}

// computes the sqrt of an interval.
// Uses the real sqrt function to compute the end points
dyadic_interval sqrt(const dyadic_interval& x){
	dyadic_interval ans;
	DYADIC p1 = sqrt(REAL(x.left));
	DYADIC p2 = sqrt(REAL(x.right));
	ans.left = std::min(p1, p2);
	ans.right = std::max(p1,p2);
	return ans;
}

// checks if the other interval is contained in the other
bool dyadic_interval::contains(const dyadic_interval& x) const{
	return this->left <= x.left && this->right >= x.right;
}

// checks if point is contained in the other
bool dyadic_interval::contains(const DYADIC& x) const{
	return this->left <= x && this->right >= x;
}

// returns interval length
DYADIC dyadic_interval::size() const{
	return this->right - this->left;
}

// distance of dyadic point from interval
DYADIC dyadic_interval::dist(const DYADIC& x) const{
	if(this->contains(x)) return 0;
	return std::min(abs(this->left-x), abs(this->right-x));
}

// max distance of dyadic point from any point inside the interval
DYADIC dyadic_interval::max_dist(const DYADIC& x) const{
	return std::max(abs(this->left-x), abs(this->right-x));
}