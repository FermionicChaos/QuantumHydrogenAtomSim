#ifndef VECPATHR2_H
#define VECPATHR2_H
#pragma once

#include "R1.h"

#include "R2.h"
#include "VecFieldR2.h"

class VecPathR2 {
private:
	R1 x, y; //Parameter dependent e.g. time, legnth, whatever.
public:
	VecPathR2();
	~VecPathR2();
	
	VecPathR2(const R1& a,const R1& b);

};
#endif
