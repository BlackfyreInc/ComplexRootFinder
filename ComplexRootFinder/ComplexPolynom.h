#pragma once

#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;
#include <cstring>
#include <cstdlib>
#include <complex>

#ifndef COMPLEXPOLYNOM_H
#define COMPLEXPOLYNOM_H

class ComplexPolynom
{
public:
    int degree;
    complex <double>* koef;

    ComplexPolynom();
    ComplexPolynom(int d, complex <double>* k);
    ~ComplexPolynom();
    void deletePoly();

    void SetPoly(int d, complex <double>* k);
    void copyPoly(ComplexPolynom& ob);
    ComplexPolynom getcopyPoly();
    ComplexPolynom operator/(ComplexPolynom& ob);
    ComplexPolynom operator*(ComplexPolynom& ob);
    ComplexPolynom operator*(complex <double> c);
    ComplexPolynom operator-(ComplexPolynom& ob);
    ComplexPolynom operator+(ComplexPolynom& ob);
    ComplexPolynom operator%(ComplexPolynom& ob);
    ComplexPolynom& operator=(const ComplexPolynom& ob);
    void printPoly();
    void deleteBy_a0();
    ComplexPolynom derivative();
    complex <double> value_complex(complex <double> x);
    void removeMultipleRoots(void);
    ComplexPolynom TransformToCauchyPoly();
    complex <double> SolveNewton(ComplexPolynom cp, complex <double> left, complex <double> right);
    complex <double> get_s();
    ComplexPolynom NomrPoly();
    complex <double> findOneRoot_Newton_CP(ComplexPolynom poly, complex <double> left, complex <double> right);
};
#endif