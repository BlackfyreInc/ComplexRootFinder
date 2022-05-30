#include "Polynom.h"
#include <complex>

#define NOROOT 11111111111
#define EPSILON 0.00000001

class RootSeeker
{
private:
    double* roots;
    complex <double>* croots;
    Polynom p;
    double a, b;
    int cr;
    bool correctPoly;

public:
    RootSeeker() {};

    void SetPoly(int d, double* coef);
    Polynom SetPoly_(int d, double* coef);
    void SetPoly(Polynom poly);
    void printPoly();
    void SetInterval(double a, double b);
    void countOfRoots();
    double* getRoots();
    void printRoots();

    static double findOneRoot_FalsePosition(Polynom& poly, double left, double right);
    static double findOneRoot_Newton(Polynom& poly, double left, double right);
    static double findOneRoot_Secant(Polynom& poly, double x0, double x1);
    static double findOneRoot_Muller(Polynom& poly, double left, double right);
    static double findOneRoot_Halley(Polynom& poly, double left, double right);
    static complex <double> findOneRoot_Secant_Complex(Polynom& poly, complex <double> x0, complex <double> x1);
    static complex <double> findOneRoot_Newton_Complex(Polynom poly, complex <double> left, complex <double> right);
    static complex <double> findOneRoot_FalsePosition_Complex(Polynom& poly, complex <double> left, complex <double> right);
    static complex <double>* DurandKerner(Polynom& poly, int deg);

    void findRoots(double (*methodName)(Polynom&, double, double));
};