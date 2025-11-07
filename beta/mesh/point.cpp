#include "stdafx.h"

#include "point.hpp"
#include "node.hpp"
#include <cmath>

double Point::epsilon = 1e-5;


Point::Point(Node &node) : xi(x), eta(y), zeta(z),x1(x), x2(y), x3(z)
{
	x = node.x;
	y = node.y;
	z = node.z;
}

Point Point::operator+(Point &p)
{
	Point tmp;
	tmp.x = x + p.x;
	tmp.y = y + p.y;
	tmp.z = z + p.z;
	return tmp;   //returns copy of object, not pointer to object
}
Point &Point::operator+=(const Point &p1)
{
	x+=p1.x;
	y+=p1.y;
	z+=p1.z;
	return (*this);
}
Point &Point::operator-=(const Point &p1)
{
	x-=p1.x;
	y-=p1.y;
	z-=p1.z;
	return (*this);
}
bool Point::operator==(const Point &p) const
{
	if(fabs(x-p.x)<epsilon && fabs(y-p.y)<epsilon && fabs(z-p.z)<epsilon) 
		return true;
	return false;
}
bool Point::operator<(const Point &p) const
{
	if(x<p.x) return true;
	if(fabs(x-p.x)<epsilon) {
		if(y<p.y) return true;
		if(fabs(y-p.y)<epsilon)
			if(z<p.z) return true;
	}
	return false;
}

ostream & operator <<(ostream& o, const Point &p)
{
	o << '(' << p.x << ", " << p.y << ", " << 
			p.z << ')';
	return o;
} 
istream & operator >>(istream& i, Point &p)
{
	i >> p.x >> p.y >> p.z;
	if(i.fail()) cerr << "Bad Point data." << endl;
	return i;
}

Point Point::rotate(int axis, double angle)
{
	Point tmp;
	double c=cos(angle *M_PI/180);
	double s=sin(angle *M_PI/180);

	switch(axis){
	case 1:	
			tmp.x=x;
			tmp.y=c*y + s*z;
			tmp.z=-s*y + c*z;
			break;
	case 2:	
			tmp.x=c*x -s*z;
			tmp.y=y;
			tmp.z=s*x + c*z;
			break;
	case 3:	
			tmp.x=c*x + s*y;
			tmp.y=-s*x + c*y;
			tmp.z=z;
			break;	
	}
	return tmp;   //returns copy of object, not pointer to object
}
