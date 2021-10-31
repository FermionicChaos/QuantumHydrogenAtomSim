#ifndef FUNCTIONSR2_H
#define FUNCTIONSR2_H
#pragma once

//#include "stdafx.h"
#include "R2.h"
#include "cR2.h"
//#include "C2.h"

cR2 Plm(int l, int m, cR2& arg);
cR2 Ylm(int l, int m, cR2& arg);
cR2 Lka(int k, double a, cR2& arg);
cR2 Rnl(double n, double l, cR2& arg);
cR2 Psi_nlm(int n, int l, int m, cR2& xx, cR2& yy);





R2 Plm(int l, int m, const R2& arg);
R2 Ylm(int l, int m, const R2& arg);
R2 Lka(int k, double a, const R2& arg);
R2 Rnl(double n, double l, const R2& arg);
R2 Psi_nlm(int n, int l, int m, R2& xx, R2& yy);

R2 exp(const R2& arg);
R2 sqrt(const R2& arg);
R2 pow(const R2& arg, double exponent);

R2 atan2(const R2& y, const R2& x);
#endif