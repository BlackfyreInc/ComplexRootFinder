#include "Polynom.h"

Polynom::Polynom()
{
    degree = 0;
    koef = new double[degree + 1];
    koef[0] = 0;
}

Polynom::Polynom(int d, double* k)
{
    degree = d;
    koef = new double[degree + 1];
    for (int i = 0; i <= degree; i++)
    {
        koef[i] = k[i];
    }
}

Polynom::~Polynom()
{
    degree = 0;
    //delete[] koef;
}

void Polynom::deletePoly()
{
    degree = 0;
    delete[] koef;
}

void Polynom::deleteBy_a0()
{
    for (int i = 0; i <= degree; i++)
    {
        koef[i] /= koef[degree];
    }
}

Polynom Polynom::derivative()
{
    Polynom p;
    int d = degree - 1;
    double* coef;
    coef = new double[d + 1];
    for (int i = 0; i <= d; i++)
    {
        coef[i] = koef[i + 1] * (i + 1);
    }
    p.SetPoly(d, coef);
    return p;
}

void Polynom::SetPoly(int d, double* k)
{
    degree = d;
    koef = new double[degree + 1];
    for (int i = 0; i <= degree; i++)
    {
        koef[i] = k[i];
    }
}

Polynom& Polynom::operator=(const Polynom& ob)
{
    degree = ob.degree;
    koef = new double[degree + 1];
    for (int i = 0; i <= degree; i++)
    {
        koef[i] = ob.koef[i];
    }

    return *this;
}

void Polynom::copyPoly(Polynom& ob)
{
    ob.degree = degree;
    ob.koef = new double[ob.degree + 1];
    for (int i = 0; i <= ob.degree; i++)
    {
        ob.koef[i] = koef[i];
    }
}

void Polynom::printPoly()
{
    cout << koef[degree] << "x^" << degree;
    for (int i = degree - 1; i >= 0; i--)
    {
        cout << " + " << koef[i] << "x^" << i;
    }
    cout << "\n";
}

Polynom Polynom::operator /(Polynom& ob)
{

    bool inAlgoritm = true;

    Polynom temp;
    Polynom ob_1;
    Polynom ob_2;
    Polynom ob_4;

    temp.degree = degree - ob.degree;
    temp.koef = new double[temp.degree + 1];
    memset(temp.koef, 0, (temp.degree + 1) * sizeof(double));

    ob_1.degree = degree;
    ob_1.koef = new double[degree + 1];
    for (int i = degree; i >= 0; i--)
    {
        ob_1.koef[i] = koef[i];
    }

    ob_2.degree = ob.degree;
    ob_2.koef = new double[ob.degree + 1];
    for (int i = ob.degree; i >= 0; i--)
    {
        ob_2.koef[i] = ob.koef[i];
    }

    ob_4.degree = ob_1.degree;
    ob_4.koef = new double[ob_1.degree + 1];

    double mnojnik;
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

Polynom Polynom::operator %(Polynom& ob)
{

    bool inAlgoritm = true;

    Polynom temp;
    Polynom ob_1;
    Polynom ob_2;
    Polynom ob_4;

    temp.degree = degree - ob.degree;
    temp.koef = new double[temp.degree + 1];
    memset(temp.koef, 0, (temp.degree + 1) * sizeof(double));

    ob_1.degree = degree;
    ob_1.koef = new double[degree + 1];
    for (int i = degree; i >= 0; i--)
    {
        ob_1.koef[i] = koef[i];
    }

    ob_2.degree = ob.degree;
    ob_2.koef = new double[ob.degree + 1];
    for (int i = ob.degree; i >= 0; i--)
    {
        ob_2.koef[i] = ob.koef[i];
    }

    ob_4.degree = ob_1.degree;
    ob_4.koef = new double[ob_1.degree + 1];

    double mnojnik;
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

double Polynom::value(double x)
{
    double y = 0;
    for (int i = 0; i <= degree; i++)
    {
        y += koef[i] * powl(x, i);
    }
    return y;
}

complex <double> Polynom::value_complex(complex <double> x)
{
    complex <double> y(0, 0);
    for (int i = 0; i <= degree; i++)
    {
        y += koef[i] * pow(x, i);
    }
    return y;
}

void Polynom::removeMultipleRoots()
{
    Polynom p;
    this->copyPoly(p);

    Polynom r1, r2, r3;
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

    if (r3.koef[0] == 0)
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

Polynom Polynom::operator *(Polynom& ob)
{
    Polynom res;
    res.degree = degree + ob.degree;
    res.koef = new double[res.degree + 1];

    memset(res.koef, 0, (res.degree + 1) * sizeof(double));

    for (int i = 0; i <= ob.degree; i++)
    {
        for (int j = 0; j <= degree; j++)
        {
            res.koef[i + j] += ob.koef[i] * koef[j];
        }
    }

    return res;
}

Polynom Polynom::operator -(Polynom& ob)
{
    Polynom res;
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
    res.koef = new double[res.degree + 1];

    memset(res.koef, 0, (res.degree + 1) * sizeof(double));

    for (int i = 0; i <= degree; i++)
    {
        for (int j = 0; j <= ob.degree; j++)
        {
            res.koef[i + j] += koef[i] - ob.koef[j];
        }
    }

    return res;
}


