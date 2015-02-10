#ifndef _DYADIC_INTERVAL_H
#define _DYADIC_INTERVAL_H
#include "iRRAM/lib.h"
#include "iRRAM/core.h"
class dyadic_interval {
	public:
		iRRAM::DYADIC left;
		iRRAM::DYADIC right;
		dyadic_interval() {};
		dyadic_interval(const iRRAM::DYADIC& left, const iRRAM::DYADIC& right) : left(left), right(right) {};
		dyadic_interval operator+(const dyadic_interval& x) const;
		dyadic_interval operator+(const iRRAM::DYADIC& y) const;
		friend dyadic_interval operator+(const iRRAM::DYADIC& x, const dyadic_interval& y);
		dyadic_interval operator-() const;
		dyadic_interval operator-(const dyadic_interval& x) const;
		dyadic_interval operator-(const iRRAM::DYADIC& x) const;
		friend dyadic_interval operator-(const iRRAM::DYADIC& x, const dyadic_interval& y);
		dyadic_interval operator*(const dyadic_interval& x) const;
		dyadic_interval operator*(const iRRAM::DYADIC& x) const;
		friend dyadic_interval operator*(const iRRAM::DYADIC& x, const dyadic_interval& y);
		//dyadic_interval operator/(const dyadic_interval& x);
		friend dyadic_interval sqrt(const dyadic_interval& x);
		bool contains(const dyadic_interval& x) const;
		bool contains(const iRRAM::DYADIC& x) const;
		iRRAM::DYADIC size() const;
		iRRAM::DYADIC dist(const iRRAM::DYADIC& x) const;
		iRRAM::DYADIC max_dist(const iRRAM::DYADIC& x) const;
		
};
#endif