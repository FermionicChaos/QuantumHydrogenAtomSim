#ifndef GeodesicFuncs_H
#define GeodesicFuncs_H
#pragma once

#include "R1.h"
#include "VecR4.h"

void EvalLE(VecR4& pos, VecR4& vel);

double Rdot(int i, double r);

double PHIdot(double r);

void PHI(const R1& r, const R1& phidot, R1& phi);

double Tdot(double r);

void T(R1& tdot, R1& t);

#endif // GeodesicFuncs_H
