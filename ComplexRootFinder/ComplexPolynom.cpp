#include "ComplexPolynom.h"

ComplexPolynom::ComplexPolynom()
{
    degree = 0;
    koef = new complex <double>[degree + 1];
    koef[0] = 0;
}

ComplexPolynom::ComplexPolynom(int d, complex <double>* k)
{
    degree = d;
    koef = new complex <double>[degree + 1];
    for (int i = 0; i <= degree; i++)
    {
        koef[i] = k[i];
    }
}

ComplexPolynom::~ComplexPolynom()
{
    degree = 0;
    //delete[] koef;
}

void ComplexPolynom::deletePoly()
{
    degree = 0;
    delete[] koef;
}

void ComplexPolynom::deleteBy_a0()
{
    for (int i = 0; i <= degree; i++)
    {
        koef[i] /= koef[degree];
    }
}

ComplexPolynom ComplexPolynom::derivative()
{
    ComplexPolynom cp;
    int d = degree - 1;
    complex <double>* coef;
    coef = new complex <double>[d + 1];
    for (int i = 0; i <= d; i++)
    {
        coef[i] = koef[i + 1] * (double)(i + 1);
    }
    cp.SetPoly(d, coef);
    return cp;
}

void ComplexPolynom::SetPoly(int d, complex <double>* k)
{
    degree = d;
    koef = new complex <double>[degree + 1];
    for (int i = 0; i <= degree; i++)
    {
        koef[i] = k[i];
    }
}

ComplexPolynom& ComplexPolynom::operator=(const ComplexPolynom& ob)
{
    degree = ob.degree;
    koef = new complex <double>[degree + 1];
    for (int i = 0; i <= degree; i++)
    {
        koef[i] = ob.koef[i];
    }

    return *this;
}

void ComplexPolynom::copyPoly(ComplexPolynom& ob)
{
    ob.degree = degree;
    ob.koef = new complex <double>[ob.degree + 1];
    for (int i = 0; i <= ob.degree; i++)
    {
        ob.koef[i] = koef[i];
    }
}

ComplexPolynom ComplexPolynom::getcopyPoly()
{
    ComplexPolynom res;
    res.degree = degree;
    res.koef = new complex <double>[res.degree + 1];
    for (int i = 0; i <= res.degree; i++)
    {
        res.koef[i] = koef[i];
    }
    return res;
}

void ComplexPolynom::printPoly()
{
    cout << "(" << koef[degree].real() << " + " << koef[degree].imag() << " * i" << ")" << "* x^" << degree;
    for (int i = degree - 1; i >= 0; i--)
    {
        cout << " + " << "(" << koef[i].real() << " + " << koef[i].imag() << " * i" << ")" << "*" << "x^" << i;
    }
    cout << "\n";
}

ComplexPolynom ComplexPolynom::operator /(ComplexPolynom& ob)
{

    bool inAlgoritm = true;

    ComplexPolynom temp;
    ComplexPolynom ob_1;
    ComplexPolynom ob_2;
    ComplexPolynom ob_4;

    temp.degree = degree - ob.degree;
    temp.koef = new complex <double>[temp.degree + 1];
    memset(temp.koef, 0, (temp.degree + 1) * sizeof(complex <double>));

    ob_1.degree = degree;
    ob_1.koef = new complex <double>[degree + 1];
    for (int i = degree; i >= 0; i--)
    {
        ob_1.koef[i] = koef[i];
    }

    ob_2.degree = ob.degree;
    ob_2.koef = new complex <double>[ob.degree + 1];
    for (int i = ob.degree; i >= 0; i--)
    {
        ob_2.koef[i] = ob.koef[i];
    }

    ob_4.degree = ob_1.degree;
    ob_4.koef = new complex <double>[ob_1.degree + 1];

    complex <double> mnojnik;
    int k = 0;
    int i, j;
    while (inAlgoritm)
    {
        for (int i = ob.degree; i >= 0; i--)
        {
            ob_4.koef[i] = ob.koef[i];
        }

        if (ob_2.degree < ob_1.degree)
        {
            for (i = ob_1.degree, j = ob_2.degree; i >= 0; i--, j--)
            {
                if (j < 0)
                {
                    ob_4.koef[i] = 0;
                }
                else
                {
                    ob_4.koef[i] = ob_2.koef[j];
                }
            }
        }

        mnojnik = ob_1.koef[ob_1.degree] / ob_4.koef[ob_1.degree];

        temp.koef[temp.degree - k] = mnojnik;
        k++;

        for (int i = 0; i <= ob_1.degree; i++)
        {
            ob_4.koef[i] *= mnojnik;
        }

        for (int i = 0; i <= ob_1.degree; i++)
        {
            ob_1.koef[i] -= ob_4.koef[i];
        }

        ob_1.degree--;
        if (ob_2.degree > ob_1.degree)
        {
            inAlgoritm = false;
        }
    }

    return temp;
}

ComplexPolynom ComplexPolynom::operator %(ComplexPolynom& ob)
{

    bool inAlgoritm = true;

    ComplexPolynom temp;
    ComplexPolynom ob_1;
    ComplexPolynom ob_2;
    ComplexPolynom ob_4;

    temp.degree = degree - ob.degree;
    temp.koef = new complex <double>[temp.degree + 1];
    memset(temp.koef, 0, (temp.degree + 1) * sizeof(complex <double>));

    ob_1.degree = degree;
    ob_1.koef = new complex <double>[degree + 1];
    for (int i = degree; i >= 0; i--)
    {
        ob_1.koef[i] = koef[i];
    }

    ob_2.degree = ob.degree;
    ob_2.koef = new complex <double>[ob.degree + 1];
    for (int i = ob.degree; i >= 0; i--)
    {
        ob_2.koef[i] = ob.koef[i];
    }

    ob_4.degree = ob_1.degree;
    ob_4.koef = new complex <double>[ob_1.degree + 1];

    complex <double> mnojnik;
    int k = 0;
    int i, j;
    while (inAlgoritm)
    {
        for (int i = ob.degree; i >= 0; i--)
        {
            ob_4.koef[i] = ob.koef[i];
        }

        if (ob_2.degree < ob_1.degree)
        {
            for (i = ob_1.degree, j = ob_2.degree; i >= 0; i--, j--)
            {
                if (j < 0)
                {
                    ob_4.koef[i] = 0;
                }
                else
                {
                    ob_4.koef[i] = ob_2.koef[j];
                }
            }
        }

        mnojnik = ob_1.koef[ob_1.degree] / ob_4.koef[ob_1.degree];

        temp.koef[temp.degree - k] = mnojnik;
        k++;

        for (int i = 0; i <= ob_1.degree; i++)
        {
            ob_4.koef[i] *= mnojnik;
        }

        for (int i = 0; i <= ob_1.degree; i++)
        {
            ob_1.koef[i] -= ob_4.koef[i];
        }

        ob_1.degree--;
        if (ob_2.degree > ob_1.degree)
        {
            inAlgoritm = false;
        }
    }

    return ob_1;
}

complex <double> ComplexPolynom::value_complex(complex <double> x)
{
    complex <double> y(0, 0);
    for (int i = 0; i <= degree; i++)
    {
        y += koef[i] * pow(x, i);
    }
    return y;
}

void ComplexPolynom::removeMultipleRoots()
{
    ComplexPolynom p;
    this->copyPoly(p);

    ComplexPolynom r1, r2, r3;
    p.copyPoly(r1);
    r2 = p.derivative();

    while (true)
    {
        r3 = r1 % r2;

        if (r3.degree <= 0)
        {
            break;
        }

        r1.deletePoly();
        r2.copyPoly(r1);
        r2.deletePoly();
        r3.copyPoly(r2);
        r3.deletePoly();
    }

    if (r3.koef[0] == 0.)
    {
        if (r2.degree > 0)
        {
            this->deletePoly();
            *this = p / r2;
        }
    }
    else
    {
        if (r3.degree > 0)
        {
            this->deletePoly();
            *this = p / r3;
        }
    }
}

ComplexPolynom ComplexPolynom::operator *(ComplexPolynom& ob)
{
    ComplexPolynom res;
    res.degree = degree + ob.degree;
    res.koef = new complex <double>[res.degree + 1];

    memset(res.koef, 0, (res.degree + 1) * sizeof(complex <double>));

    for (int i = 0; i <= ob.degree; i++)
    {
        for (int j = 0; j <= degree; j++)
        {
            res.koef[i + j] += ob.koef[i] * koef[j];
        }
    }

    return res;
}

ComplexPolynom ComplexPolynom::operator *(complex <double> c)
{
    ComplexPolynom res;
    res.degree = degree;
    res.koef = new complex <double>[res.degree + 1];

    memset(res.koef, 0, (res.degree + 1) * sizeof(complex <double>));

    for (int i = 0; i <= degree; i++)
    {
        res.koef[i] += c * koef[i];
    }

    return res;
}

ComplexPolynom ComplexPolynom::operator -(ComplexPolynom& ob)
{
    ComplexPolynom res;
    res.degree = degree;
    if (degree - ob.degree > 0)
    {
        res.degree = degree;
    }
    else if (degree - ob.degree < 0)
    {
        res.degree = ob.degree;
    }
    else
    {
        res.degree = degree;
        for (int i = 0; i <= degree; i++)
        {
            if (koef[i] == ob.koef[i])
            {
                res.degree -= 1;
            }
        }
    }
    res.koef = new complex <double>[res.degree + 1];

    memset(res.koef, 0, (res.degree + 1) * sizeof(complex <double>));

    for (int i = 0; i <= degree; i++)
    {
        for (int j = 0; j <= ob.degree; j++)
        {
            res.koef[i + j] += koef[i] - ob.koef[j];
        }
    }

    return res;
}

ComplexPolynom ComplexPolynom::operator +(ComplexPolynom& ob)
{
    ComplexPolynom res;
    res.degree = degree;
    if (degree > ob.degree)
    {
        res.degree = degree;
    }
    else if (degree < ob.degree)
    {
        res.degree = ob.degree;
    }
    else
    {
        res.degree = degree;
        for (int i = 0; i <= degree; i++)
        {
            if (koef[i] == -ob.koef[i])
            {
                res.degree -= 1;
            }
        }
    }
    res.koef = new complex <double>[res.degree + 1];

    memset(res.koef, 0, (res.degree + 1) * sizeof(complex <double>));

    for (int i = 0; i <= degree; i++)
    {
        for (int j = 0; j <= ob.degree; j++)
        {
            res.koef[i + j] += koef[i] + ob.koef[j];
        }
    }

    return res;
}

ComplexPolynom ComplexPolynom::TransformToCauchyPoly() 
{
    ComplexPolynom cauchy_poly;
    cauchy_poly.degree = degree;
    bool normalize = false;
    complex <double> norm_c = 0., val;
    for (int i = 0; i < degree + 1; i++) 
    {
        if (i == 0)
        {
            val = koef[i];
            if (abs(val) != 1)
            {
                normalize = true;
                norm_c = val;
            }
            cauchy_poly.koef[i] = 1;
        }
        else if (i == degree) 
        {
            val = koef[i];
            if (normalize) 
            {
                val /= 1. * norm_c;
            }
            val = -val;
            cauchy_poly.koef[i] = val;
        }
        else 
        {
            val = koef[i];
            if (normalize)
            {
                val /= 1. * norm_c;
            }
            cauchy_poly.koef[i] = val;
        }
    }
    return cauchy_poly;
}

complex <double> ComplexPolynom::SolveNewton(ComplexPolynom cp, complex <double> left, complex <double> right)
{
    ComplexPolynom pd = cp.derivative();

    complex <double> x0 = (left + right) / 2.0;

    if (pd.value_complex(x0) == 0.0)
    {
        cout << "Method diverges" << endl;
        return -1;
    }

    complex <double> xn = x0 - cp.value_complex(x0) / pd.value_complex(x0);

    while (abs(xn - x0) > 0.00001)
    {
        x0 = xn;
        xn = xn - cp.value_complex(xn) / pd.value_complex(xn);
    }
    if (fabs(xn.real()) < 0.00001)
    {
        xn.real(0);
    }
    else if (fabs(xn.imag()) < 0.00001)
    {
        xn.imag(0);
    }
    return xn;
}

complex <double> ComplexPolynom::get_s()
{
    ComplexPolynom cauchy = TransformToCauchyPoly();
    complex <double> c(0.3, 0.5), c1(0.7, 0.9), c2(0.0, 1.0);
    complex <double> beta = SolveNewton(cauchy, c, c1);
    double r = ((double)rand() / (RAND_MAX)) * 2 * 3.14;
    return beta * exp(r * c2);
}

ComplexPolynom ComplexPolynom::NomrPoly() 
{
    ComplexPolynom res;
    res.degree = degree;
    for (int i = 0; i < degree + 1; i++) 
    {
        res.koef[i] = (koef[i]) / (koef[0]);
    }
    return res;
}