#include "RootSeeker.h"

void RootSeeker::SetPoly(int d, double* coefs)
{
    p.deletePoly();

    int degr = d;
    for (int i = degr; i > 0; i--)
    {
        if (coefs[i] == 0)
        {
            d--;
        }
        else
        {
            break;
        }
    }
    if (d == 0)
    {
        cout << "Your polynom is incorrect\n";
        correctPoly = false;
        p.SetPoly(degr, coefs);
    }
    else
    {
        correctPoly = true;
        p.SetPoly(d, coefs);
    }
}

Polynom RootSeeker::SetPoly_(int d, double* coefs)
{
    p.deletePoly();

    int degr = d;
    for (int i = degr; i > 0; i--)
    {
        if (coefs[i] == 0)
        {
            d--;
        }
        else
        {
            break;
        }
    }
    if (d == 0)
    {
        cout << "Your polynom is incorrect\n";
        correctPoly = false;
        p.SetPoly(degr, coefs);
    }
    else
    {
        correctPoly = true;
        p.SetPoly(d, coefs);
        return p;
    }
}

void RootSeeker::SetPoly(Polynom poly)
{
    p.deletePoly();

    int d = poly.degree;
    int degr = d;
    for (int i = degr; i > 0; i--)
    {
        if (poly.koef[i] == 0)
        {
            d--;
        }
        else
        {
            break;
        }
    }
    if (d == 0)
    {
        cout << "Your polynom is incorrect\n";
        correctPoly = false;
        p.SetPoly(degr, poly.koef);
    }
    else
    {
        correctPoly = true;
        p = poly;
    }
}

void RootSeeker::printPoly()
{
    p.printPoly();
}

void RootSeeker::SetInterval(double a, double b)
{
    this->a = a;
    this->b = b;
}

void RootSeeker::countOfRoots()
{
    if (!correctPoly)
    {
        return;
    }

    // Sturm's theorem

    Polynom* hi;
    hi = new Polynom[p.degree + 1];
    hi[0] = p;
    hi[1] = p.derivative();

    for (int i = 2; i <= p.degree; i++)
    {
        hi[i] = hi[i - 2] % hi[i - 1];
        for (int j = 0; j <= hi[i].degree; j++)
        {
            hi[i].koef[j] *= -1;
        }
    }

    int c1, c2;
    double y1, y2, y3, y4;
    c1 = 0;
    c2 = 0;
    y1 = hi[0].value(a);
    y3 = hi[0].value(b);

    for (int i = 1; i <= p.degree; i++)
    {

        y2 = hi[i].value(a);
        y4 = hi[i].value(b);

        if (!((y1 > 0 && y2 > 0) || (y1 < 0 && y2 < 0)))
        {
            c1++;
        }
        if (!((y3 > 0 && y4 > 0) || (y3 < 0 && y4 < 0)))
        {
            c2++;
        }

        y1 = y2;
        y3 = y4;
    }

    this->cr = c1 - c2;
    cout << "Counts of the real roots on interval [" << a << ", " << b << "]: " << this->cr << "\n";
    roots = new double[cr];
}
/// <summary>
/// False Position method for finding one root of a polynomial.
/// </summary>
/// <param name="poly"></param>
/// <param name="left"></param>
/// <param name="right"></param>
/// <returns>
/// Returns <see cref="double a"/> in case of successful cycle.
/// Returns <see cref="NOROOT"/> and prints "This method doesn't work for this polynom." in case of unsuccessful cycle.
/// </returns>
double RootSeeker::findOneRoot_FalsePosition(Polynom& poly, double left, double right)
{
    double a = left;
    double b = right;
    double c = (a * poly.value(b) - b * poly.value(a)) / (poly.value(b) - poly.value(a));

    int i = 0;
    while (fabs(poly.value(a)) > EPSILON)
    {
        if (poly.value(c) < EPSILON)
        {
            a = c;
        }
        else
        {
            b = c;
        }

        c = (a * poly.value(b) - b * poly.value(a)) / (poly.value(b) - poly.value(a));


        i++;
        if (i > 10000)
        {
            cout << "This method doesn't work for this polynom." << endl;
            return NOROOT;
        }
    }

    return a;
}
/// <summary>
/// Newton method for finding one real root of a polynomial.
/// </summary>
/// <param name="poly"></param>
/// <param name="left"></param>
/// <param name="right"></param>
/// <returns>
/// Returns <see cref="double xn"/> in case of successful cycle.
/// </returns>
double RootSeeker::findOneRoot_Newton(Polynom& poly, double left, double right)
{
    Polynom pd = poly.derivative();

    double x0 = (left + right) / 2;
    if (x0 == 0)
    {
        x0 = (left + right) / 2;
        if (x0 == 0)
        {
            x0 -= 5;
        }
    }

    double xn = x0 - poly.value(x0) / pd.value(x0);

    while (fabs(xn - x0) > EPSILON)
    {
        x0 = xn;
        xn = xn - poly.value(xn) / pd.value(xn);
    }

    return xn;
}
/// <summary>
/// Secant method for finding one real root of a polynomial.
/// </summary>
/// <param name="poly"></param>
/// <param name="x0"></param>
/// <param name="x1"></param>
/// <returns>
/// Returns <see cref="double xn"/> and prints "Method diverges" in case of successful cycle.
/// </returns>
double RootSeeker::findOneRoot_Secant(Polynom& poly, double x0, double x1)
{
    double xn, f0, f1, fn = 1, e = 1;

    while (e > 0.005 && fn != 0)
    {
        f0 = poly.value(x0);
        f1 = poly.value(x1);
        xn = (x1 - (f1 * (x1 - x0) / (f1 - f0)));
        fn = poly.value(xn);
        e = fabs((x1 - x0) / x1);
        x0 = x1;
        x1 = xn;
    }

    return xn;
}
/// <summary>
/// Muller method for finding one real root of a polynomial.
/// </summary>
/// <param name="poly"></param>
/// <param name="x0"></param>
/// <param name="x1"></param>
/// <param name="x2"></param>
/// <returns>
/// Returns <see cref="double res"/> in case of successful cycle.
/// Returns <see cref="-1"/> and prints "Method diverges. More iterations may be needed." in case of unsuccessful cycle.
/// </returns>
double RootSeeker::findOneRoot_Muller(Polynom& poly, double left, double right)
{
    if (left == 0 || right == 0 || left == right) 
    {
        left += 5;
        right += 10;
    }
    double x2 = left + right + 10;
    int i = 1, j = 0;
    double res, t, h4, f1, f2, f3, d1, d2, h1, h2, a0, a1, a2, x, y;

    do
    {
        f1 = poly.value(left);
        f2 = poly.value(right);
        f3 = poly.value(x2);
        d1 = f1 - f3;
        d2 = f2 - f3;
        h1 = left - x2;
        h2 = right - x2;
        a0 = f3;
        a1 = ((d2 * h1 * h1 - d1 * h2 * h2) / (h1 * h2 * (h1 - h2)));
        a2 = ((d1 * h2 - d2 * h1) / (h1 * h2 * (h1 - h2)));
        x = -((2 * a0) / (a1 + sqrt(fabs(a1 * a1 - (4 * a2 * a0)))));
        y = -((2 * a0) / (a1 - sqrt(fabs(a1 * a1 - (4 * a2 * a0)))));

        if (fabs((a1 + sqrt(fabs(a1 * a1 - (4 * a2 * a0))))) > fabs((a1 - sqrt(fabs(a1 * a1 - (4 * a2 * a0))))))
            h4 = x;
        else
            h4 = y;
        res = x2 + h4;
        if (fabs(poly.value(res)) < 0.0001)
        {
            i = 0;
        }
        else
        {
            t = res;
            left = right;
            right = x2;
            x2 = res;
        }
        if (j > 10000) 
        {
            break;
        }
        if (left != left) 
        {
            break;
        }
        j++;
    } while (i != 0);

    return res;
}
/// <summary>
/// Halley method for finding one real root of a polynomial.
/// </summary>
/// <param name="poly"></param>
/// <param name="left"></param>
/// <param name="right"></param>
/// <returns>
/// Returns <see cref="double xn"/> in case of success.
/// </returns>
double RootSeeker::findOneRoot_Halley(Polynom& poly, double left, double right)
{
    Polynom pd = poly.derivative();
    Polynom pdd = pd.derivative();

    double x0 = (left + right) / 2;
    if (x0 == 0)
    {
        x0 = (left + right) / 2;
        if (x0 == 0)
        {
            x0 -= 5;
        }
    }

    double xn = x0 - (2*poly.value(x0)*pd.value(x0)) / (2*(pow(pd.value(x0), 2)) - poly.value(x0) * pdd.value(x0));

    while (fabs(xn - x0) > EPSILON)
    {
        x0 = xn;
        xn = xn - (2 * poly.value(xn) * pd.value(xn)) / (2 * (pow(pd.value(xn), 2)) - poly.value(xn) * pdd.value(xn));
    }

    return xn;
}
/// <summary>
/// Halley method for finding one complex root of a polynomial.
/// </summary>
/// <param name="poly"></param>
/// <param name="left"></param>
/// <param name="right"></param>
/// <returns>
/// Returns <see cref="complex <double> xn"/> in case of success.
/// </returns>
complex <double> RootSeeker::findOneRoot_Halley_Complex(Polynom& poly, complex <double> left, complex <double> right)
{
    Polynom pd = poly.derivative();
    Polynom pdd = pd.derivative();

    complex <double> x0 = (left + right) / 2.;
    if (x0 == complex <double>(0,0))
    {
        x0 = (left + right) / 2.;
        if (x0 == 0.0)
        {
            x0 -= 5;
        }
    }

    complex <double> xn = x0 - (2. * poly.value_complex(x0) * pd.value_complex(x0)) / (2. * (pow(pd.value_complex(x0), 2)) - poly.value_complex(x0) * pdd.value_complex(x0));

    while (abs(xn - x0) > EPSILON)
    {
        x0 = xn;
        xn = xn - (2. * poly.value_complex(xn) * pd.value_complex(xn)) / (2. * (pow(pd.value_complex(xn), 2)) - poly.value_complex(xn) * pdd.value_complex(xn));
    }
    if (fabs(xn.real()) < EPSILON)
    {
        xn.real (0);
    }
    else if (fabs(xn.imag()) < EPSILON)
    {
        xn.imag (0);
    }
    return xn;
}
/// <summary>
/// Newton–Raphson method for finding one complex root of polynomial.
/// </summary>
/// <param name="poly"></param>
/// <param name="x0"></param>
/// <returns>
/// Returns <see cref="complex <double> xn"/> in case of successful cycle.
/// Returns <see cref="-1"/> and prints "Method diverges" in case of unsuccessful test.
/// </returns>
complex <double> RootSeeker::findOneRoot_Newton_Complex(Polynom poly, complex <double> left, complex <double> right)
{
    Polynom pd = poly.derivative();

    complex <double> x0 = (left + right) / 2.0;
    
    if (pd.value_complex(x0) == 0.0) 
    {
        cout << "Method diverges" << endl;
        return -1;
    }

    complex <double> xn = x0 - poly.value_complex(x0) / pd.value_complex(x0);

    while (abs(xn - x0) > EPSILON)
    {
        x0 = xn;
        xn = xn - poly.value_complex(xn) / pd.value_complex(xn);
    }
    if (fabs(xn.real()) < EPSILON)
    {
        xn.real (0);
    }
    else if (fabs(xn.imag()) < EPSILON)
    {
        xn.imag (0);
    }
    return xn;
}
/// <summary>
/// Secant method for finding one complex root of polynomial.
/// </summary>
/// <param name="poly"></param>
/// <param name="x0"></param>
/// <param name="x1"></param>
/// <returns>
/// Returns <see cref="complex <double> xn"/> in case of successful cycle.
/// Returns <see cref="-1"/> and prints "Method diverges" in case of unsuccessful test.
/// </returns>
complex <double> RootSeeker::findOneRoot_Secant_Complex(Polynom& poly, complex <double> x0, complex <double> x1)
{
    complex <double> xn, f0, f1, fn = 1;
    double e = 1;

    while (e > 0.005 && fn != 0.0)
    {
        f0 = poly.value_complex(x0);
        f1 = poly.value_complex(x1);
        xn = (x1 - (f1 * (x1 - x0) / (f1 - f0)));
        fn = poly.value_complex(xn);
        e = abs((x1 - x0) / x1);
        x0 = x1;
        x1 = xn;
    }
    if (fabs(xn.real()) < EPSILON)
    {
        xn.real(0);
    }
    else if (fabs(xn.imag()) < EPSILON)
    {
        xn.imag(0);
    }
    return xn;
}
/// <summary>
/// FalsePosition method for finding one complex root of polynomial.
/// </summary>
/// <param name="poly"></param>
/// <param name="left"></param>
/// <param name="right"></param>
/// <returns>
/// Returns <see cref="complex <double> a"/> in case of successful check.
/// Returns <see cref="NOROOT"/> and prints "This method doesn't work for this polynom." in case of unsuccessful check.
/// </returns>
complex <double> RootSeeker::findOneRoot_FalsePosition_Complex(Polynom& poly, complex <double> left, complex <double> right)
{
    complex <double> a = left;
    complex <double> b = right;
    complex <double> c = (a * poly.value_complex(b) - b * poly.value_complex(a)) / (poly.value_complex(b) - poly.value_complex(a));

    int i = 0;
    while (abs(poly.value_complex(a)) > EPSILON)
    {
        if (norm(poly.value_complex(c)) < EPSILON)
        {
            a = c;
        }
        else
        {
            b = c;
        }

        c = (a * poly.value_complex(b) - b * poly.value_complex(a)) / (poly.value_complex(b) - poly.value_complex(a));


        i++;
        if (i > 10000)
        {
            cout << "This method doesn't work for this polynom." << endl;
            return NOROOT;
        }
    }
    if (fabs(a.real()) < EPSILON) 
    {
        a.real (0);
    }
    else if (fabs(a.imag()) < EPSILON)
    {
        a.imag (0);
    }
    return a;
}

complex <double> RootSeeker::findOneRoot_Muller_Complex(Polynom& poly, complex <double> left, complex <double> right)
{
    complex <double> x0 = left;
    complex <double> x2 = right;
    complex <double> x1 = (left + right) / 2.;
    int i;
    complex <double> res;

    for (i = 0;; ++i)
    {
        complex <double> f1 = poly.value_complex(x0);
        complex <double> f2 = poly.value_complex(x1);
        complex <double> f3 = poly.value_complex(x2);
        complex <double> d1 = f1 - f3;
        complex <double> d2 = f2 - f3;
        complex <double> h1 = x0 - x2;
        complex <double> h2 = x1 - x2;
        complex <double> a0 = f3;
        complex <double> a1 = (((d2 * pow(h1, 2)) - (d1 * pow(h2, 2))) / ((h1 * h2) * (h1 - h2)));
        complex <double> a2 = (((d1 * h2) - (d2 * h1)) / ((h1 * h2) * (h1 - h2)));
        complex <double> x = ((-2. * a0) / (a1 + abs(sqrt(a1 * a1 - 4. * a0 * a2))));
        complex <double> y = ((-2. * a0) / (a1 - abs(sqrt(a1 * a1 - 4. * a0 * a2))));
        double modx = norm(x);
        double mody = norm(y);
        if (modx >= mody)
            res = x + x2;
        else
            res = y + x2;

        complex <double> m = res * 100.;
        complex <double> n = x2 * 100.;
        double mm = floor(norm(m));
        double nn = floor(norm(n));
        if (mm == nn)
            break;
        x0 = x1;
        x1 = x2;
        x2 = res;
        if (i > 10000)
        {
            cout << "Method diverges. More iterations may be needed." << endl;
            break;
        }
    }
    if (fabs(res.real()) < EPSILON)
    {
        res.real(0);
    }
    else if (fabs(res.imag()) < EPSILON)
    {
        res.imag(0);
    }
    return res;
}
/// <summary>
/// One of numerical methods to solve and find roots of polynomials. Finds both real and complex roots.
/// </summary>
/// <param name="poly">Polynom class poly</param>
/// <returns>
/// Returns <see cref="rootApprox"/> list of found roots.
/// </returns>
complex <double>* RootSeeker::DurandKerner(Polynom& poly, int deg)
{
    complex<double> constant(0.4, 0.9), reset(1, 0);
    complex<double>* rootApprox = NULL;
    complex<double>* firstDenom = new complex<double>[deg];
    complex<double>* secondDenom = new complex<double>[deg];
    rootApprox = new complex<double>[deg];
    fill_n(firstDenom, deg, 1);
    fill_n(secondDenom, deg, 1);

    for (int i = 0; i < deg; i++)
    {
        rootApprox[i] = pow(constant, i);
    }

    for (int i = 0; i <= 10000; i++)
    {
        for (int xn = 0; xn < deg; xn++)
        {
            for (int firstHalf = 1; firstHalf <= deg - 1; firstHalf++)
            {
                if (firstHalf == 1)
                {
                    firstDenom[xn] = reset;
                }

                if (xn < firstHalf)
                {
                    firstDenom[xn] *= (rootApprox[xn] - rootApprox[firstHalf]);
                }
            }
            complex<double> y = poly.value_complex(rootApprox[xn]);
            if (xn > 0)
            {
                for (int secondHalf = 0; secondHalf <= deg - 1; secondHalf++)
                {
                    if (secondHalf == 0)
                    {
                        secondDenom[xn] = reset;
                    }

                    if (xn > secondHalf)
                    {
                        secondDenom[xn] *= (rootApprox[xn] - rootApprox[secondHalf]);
                    }
                }
            }
            rootApprox[xn] = (rootApprox[xn] - ((y) / (firstDenom[xn] * secondDenom[xn])));
        }
    }
    return rootApprox;
}

void RootSeeker::findRoots(double (*methodName)(Polynom&, double, double))
{
    if (!correctPoly || cr == 0)
    {
        return;
    }

    Polynom dvochlen;
    dvochlen.degree = 1;
    dvochlen.koef = new double[2];
    dvochlen.koef[1] = 1;

    Polynom cop;
    p.copyPoly(cop);
    roots[0] = methodName(cop, a, b);

    if (fabs(roots[0]) < EPSILON)
    {
        roots[0] = 0;
    }

    if (roots[0] == NOROOT)
    {
        return;
    }

    dvochlen.koef[0] = -roots[0];

    for (int i = 1; i < cr; i++)
    {
        cop = cop / dvochlen;

        roots[i] = methodName(cop, a, b);

        if (roots[i] == NOROOT)
        {
            return;
        }

        if (fabs(roots[i]) < EPSILON)
        {
            roots[i] = 0;
        }

        dvochlen.koef[0] = -roots[i];
    }

}

void RootSeeker::printRoots()
{
    if (correctPoly && cr != 0 && roots[0] != NOROOT)
    {
        cout << "Roots: ";
        for (int i = 0; i < cr; i++)
        {
            cout << roots[i] << " ";
        }
        cout << "\n";
    }
    else
    {
        cout << "No roots\n";
    }
}

double* RootSeeker::getRoots()
{
    return roots;
}


