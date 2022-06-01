#pragma once

#include "Polynom.h"
#include "ComplexPolynom.h"

#define NOROOT 11111111111
#define EPSILON 0.00000001

class RootSeeker
{
private:
    double* roots;
    complex <double>* complex_roots;
    Polynom p;
    ComplexPolynom cp;
    double a, b;
    complex <double> a_c, b_c;
    int cr;
    bool correctPoly;

public:
    RootSeeker() {};

    void SetPoly(int d, double* coef);
    void SetPoly(int d, complex <double>* coef);
    Polynom SetPoly_(int d, double* coef);
    ComplexPolynom SetPoly_(int d, complex <double>* coefs);
    void SetPoly(Polynom poly);
    void printPoly();
    void printCPoly();
    void SetInterval(double a, double b);
    void countOfRoots();
    double* getRoots();
    void printRoots();

    static double findOneRoot_FalsePosition(Polynom& poly, double left, double right);
    static double findOneRoot_Newton(Polynom& poly, double left, double right);
    static double findOneRoot_Secant(Polynom& poly, double x0, double x1);
    static double findOneRoot_Muller(Polynom& poly, double left, double right);
    static double findOneRoot_Halley(Polynom& poly, double left, double right);

    static complex <double> findOneRoot_Halley_Complex(Polynom& poly, complex <double> left, complex <double> right);
    static complex <double> findOneRoot_Halley_Complex(ComplexPolynom& poly, complex <double> left, complex <double> right);

    static complex <double> findOneRoot_Newton_Complex(Polynom& poly, complex <double> left, complex <double> right);
    static complex <double> findOneRoot_Newton_Complex(ComplexPolynom& poly, complex <double> left, complex <double> right);

    static complex <double> findOneRoot_Secant_Complex(Polynom& poly, complex <double> x0, complex <double> x1);
    static complex <double> findOneRoot_Secant_Complex(ComplexPolynom& poly, complex <double> x0, complex <double> x1);

    static complex <double> findOneRoot_FalsePosition_Complex(Polynom& poly, complex <double> left, complex <double> right);
    static complex <double> findOneRoot_FalsePosition_Complex(ComplexPolynom& poly, complex <double> left, complex <double> right);

    static complex <double> findOneRoot_Muller_Complex(Polynom& poly, complex <double> left, complex <double> right);
    static complex <double> findOneRoot_Muller_Complex(ComplexPolynom& poly, complex <double> left, complex <double> right);

    static complex <double>* DurandKerner(Polynom& poly, int deg);
    static complex <double>* DurandKerner(ComplexPolynom& poly, int deg);

    complex <double> findOneRootJenkinsTraub(ComplexPolynom c_poly);

    void findRoots(double (*methodName)(Polynom&, double, double));
};