#include <ctime>
#include <random>
#include "RootSeeker.h"


#define RAND_MAX 10000

double randomDouble(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double* randomArray(int n)
{
    double* randArray;
    randArray = new double[n];
    for (int i = 0; i < n; i++)
    {
        randArray[i] = randomDouble(-100, 100);
    }
    return randArray;
}

Polynom randomPoly(int maxDegree)
{
    double* roots;
    int n;

    while (true)
    {
        n = rand() % (maxDegree + 1);
        if (n != 0)
        {
            break;
        }
    }

    roots = randomArray(n);
    cout << "Roots: ";
    for (int i = 0; i < n; i++)
    {
        cout << roots[i] << " ";
    }
    cout << "\n";

    Polynom dvochlen, p;
    dvochlen.degree = 1;
    dvochlen.koef = new double[2];
    dvochlen.koef[1] = 1;

    dvochlen.koef[0] = -roots[0];
    p = dvochlen;
    for (int i = 1; i < n; i++)
    {
        dvochlen.koef[0] = -roots[i];
        p = p * dvochlen;
    }

    return p;

}


// Input Polynomial: {1, 2, 3, 4} -> 4x^3 + 3x^2 + 2x + 1

int main()
{
    //   Test 1. Generate random the polynomials by their roots
    cout << "\nTEST 1\n\n";

    Polynom p;
    srand(time(NULL));

    int tests = 10;
    int maxDegree = 5;

    cout << "-------------------------------------------------------------------------\n";
    for (int i = 0; i < tests; i++)
    {
        p = randomPoly(maxDegree);

        cout << "Polynom: ";
        p.printPoly();
        cout << "\n";

        RootSeeker rs1;
        rs1.SetInterval(-600, 600);

        rs1.SetPoly(p);

        rs1.countOfRoots();

        cout << "(Newton Method) " << endl;
        rs1.findRoots(rs1.findOneRoot_Newton);
        rs1.printRoots();

        cout << "(False Position Method) " << endl;
        rs1.findRoots(rs1.findOneRoot_FalsePosition);
        rs1.printRoots();

        cout << "(Secant Method) " << endl;
        rs1.findRoots(rs1.findOneRoot_Secant);
        rs1.printRoots();

        cout << "(Muller Method) " << endl;
        rs1.findRoots(rs1.findOneRoot_Muller);
        rs1.printRoots();

        cout << "(Halley Method) " << endl;
        rs1.findRoots(rs1.findOneRoot_Halley);
        rs1.printRoots();
        cout << "-------------------------------------------------------------------------\n";
    }





    //   Test 2. Find roots of the polynomials (1, 2 and 3 orders)

    cout << "\n\n\n\n\n\n\n\n\n TEST 2\n\n";

    RootSeeker rs;
    rs.SetInterval(-100, 100);

    cout << "1st order\n";
    double test_poly[][2] = { {5,1}, {-1,2}, {4,2}, {0,2}, {2,0} ,{0,0} };
    double roots_poly[6] = { -5, 0.5, -2, 0, NOROOT, NOROOT };
    for (int i = 0; i < 6; ++i)
    {
        rs.SetPoly(1, test_poly[i]);
        rs.printPoly();

        cout << "Root: ";
        if (roots_poly[i] == NOROOT)
        {
            cout << "No root\n";
        }
        else
        {
            cout << roots_poly[i] << "\n";
        }

        rs.countOfRoots();

        cout << "(Newton Method) ";
        rs.findRoots(rs.findOneRoot_Newton);
        rs.printRoots();

        cout << "(False Position Method) ";
        rs.findRoots(rs.findOneRoot_FalsePosition);
        rs.printRoots();

        cout << "(Secant Method) " << endl;
        rs.findRoots(rs.findOneRoot_Secant);
        rs.printRoots();

        cout << "(Muller Method) " << endl;
        rs.findRoots(rs.findOneRoot_Muller);
        rs.printRoots();

        cout << "(Halley Method) " << endl;
        rs.findRoots(rs.findOneRoot_Halley);
        rs.printRoots();
        cout << "\n";

    }
    cout << "-------------------------------------------------------------------------\n";


    cout << "2nd order\n";
    double test_poly2[][3] = { {1.,-3.,2.0}, {1.,-2.,1}, {2.,4, -6.}, {1,1,1}, {2,0,-4} ,{3,0,9},{0.,-3.,2.0}, {0.,-2.,1}, {0.,0., 4.}, {0,2,4}, {0,0,0} ,{3,0,0}, {-1, 0, 1} };
    double roots_poly2[][2] = { {1, 0.5}, {1, 1}, {1, -1.0 / 3}, {NOROOT, NOROOT}, {sqrt(0.5), -sqrt(0.5)}, {NOROOT, NOROOT}, {0, 1.5}, {0, 2},{0, 0}, {0, -0.5}, {NOROOT, NOROOT}, {NOROOT, NOROOT}, {-1, 1} };
    for (int i = 0; i < 13; ++i)
    {
        rs.SetPoly(2, test_poly2[i]);
        rs.printPoly();

        cout << "Roots: ";
        if (roots_poly2[i][0] == NOROOT && roots_poly2[i][1] == NOROOT)
        {
            cout << "No roots\n";
        }
        else
        {
            cout << roots_poly2[i][0] << " " << roots_poly2[i][1] << "\n";
        }

        rs.countOfRoots();

        cout << "(Newton Method) ";
        rs.findRoots(rs.findOneRoot_Newton);
        rs.printRoots();

        cout << "(False Position Method) ";
        rs.findRoots(rs.findOneRoot_FalsePosition);
        rs.printRoots();

        cout << "(Secant Method) " << endl;
        rs.findRoots(rs.findOneRoot_Secant);
        rs.printRoots();

        cout << "(Muller Method) " << endl;
        rs.findRoots(rs.findOneRoot_Muller);
        rs.printRoots();

        cout << "(Halley Method) " << endl;
        rs.findRoots(rs.findOneRoot_Halley);
        rs.printRoots();
        cout << "\n";


    }
    cout << "-------------------------------------------------------------------------\n";


    cout << "3d order\n";
    double test_poly3[][4] = { {0, -3, 2, 1}, {1, 2, 1, 0}, {56, -34, 1, 1}, {0, 1, 0, 1}, {1, 0, 0, 1}, {-6, 2, -1, -4} };
    double roots_poly3[][3] = { {0, 1, -3}, {-1, -1, NOROOT}, {-7, 4, 2}, {0, NOROOT, NOROOT}, {NOROOT, NOROOT, -1}, {NOROOT, NOROOT, -1.38835} };
    for (int i = 0; i < 6; ++i)
    {
        rs.SetPoly(3, test_poly3[i]);
        rs.printPoly();

        cout << "Roots: ";
        if (roots_poly3[i][0] == NOROOT && roots_poly3[i][1] == NOROOT && roots_poly3[i][2] == NOROOT)
        {
            cout << "No roots\n";
        }
        else
        {
            cout << roots_poly3[i][0] << " " << roots_poly3[i][1] << " " << roots_poly3[i][2] << "\n";
        }

        rs.countOfRoots();

        cout << "(Newton Method) ";
        rs.findRoots(rs.findOneRoot_Newton);
        rs.printRoots();

        cout << "(False Position Method) ";
        rs.findRoots(rs.findOneRoot_FalsePosition);
        rs.printRoots();

        cout << "(Secant Method) " << endl;
        rs.findRoots(rs.findOneRoot_Secant);
        rs.printRoots();

        cout << "(Muller Method) " << endl;
        rs.findRoots(rs.findOneRoot_Muller);
        rs.printRoots();

        cout << "(Halley Method) " << endl;
        rs.findRoots(rs.findOneRoot_Halley);
        rs.printRoots();
        cout << "\n";

    }



    //   Test 3. Generate random polynomials and calculate time work of the methods

    cout << "\n\n\n\n\n\n\n\n\n TEST 3\n\n";

    double t1 = 0;
    double t2 = 0;
    double t3 = 0;
    double t4 = 0;
    double t5 = 0;
    double* ar;
    cout.setstate(ios_base::failbit);
    for (int i = 0; i < 10; i++)
    {
        int n = rand() % 15;
        if (n == 0)
        {
            continue;
        }

        ar = new double[n];
        ar = randomArray(n);
        rs.SetPoly(n - 1, ar);

        rs.countOfRoots();

        clock_t StartT1 = clock();
        rs.findRoots(rs.findOneRoot_Newton);
        clock_t EndT1 = clock();
        t1 += (EndT1 - StartT1) / (double)CLOCKS_PER_SEC;

        clock_t StartT2 = clock();
        rs.findRoots(rs.findOneRoot_Halley);
        clock_t EndT2 = clock();
        t2 += (EndT2 - StartT2) / (double)CLOCKS_PER_SEC;

        clock_t StartT3 = clock();
        rs.findRoots(rs.findOneRoot_FalsePosition);
        clock_t EndT3 = clock();
        t3 += (EndT3 - StartT3) / (double)CLOCKS_PER_SEC;

        clock_t StartT4 = clock();
        rs.findRoots(rs.findOneRoot_Muller);
        clock_t EndT4 = clock();
        t4 += (EndT4 - StartT4) / (double)CLOCKS_PER_SEC;

        clock_t StartT5 = clock();
        rs.findRoots(rs.findOneRoot_Secant);
        clock_t EndT5 = clock();
        t5 += (EndT5 - StartT5) / (double)CLOCKS_PER_SEC;

        delete[] ar;

    }
    cout.clear();
    cout << "\n";
    cout << "Newton: " << t1 << " seconds\n";

    cout << "Halley: " << t2 << " seconds\n";

    cout << "False Position: " << t3 << " seconds\n";

    cout << "Muller: " << t4 << " seconds\n";

    cout << "Secant: " << t5 << " seconds\n";


    //   Test 4. Find complex roots of the polynomials and count time

    cout << "\n\n\n\n\n\n\n\n\n TEST 4\n\n";
    double time1 = 0;
    double time2 = 0;
    double time3 = 0;
    double time4 = 0;
    double time5 = 0;
    double time6 = 0;
    double test_poly4[][4] = { {3, 1, 3, 1}, { 6, 3, 2, 1 }, { -3, -1, -1, 2 }, { -10, 3, -8, 1 } };

    for (int i = 0; i < 4; ++i) 
    {
        Polynom p1 = rs.SetPoly_(3, test_poly4[i]);
        rs.printPoly();

        rs.countOfRoots();

        cout << "(Newton Method) ";
        complex <double> c1(1.0, 3.0);
        complex <double> c2(4.0, 5.0);
        clock_t StartNC = clock();
        complex <double> a = rs.findOneRoot_Newton_Complex(p1, c1, c2);
        clock_t EndNC = clock();
        time1 += (EndNC - StartNC) / (double)CLOCKS_PER_SEC;
        cout << a << endl;

        cout << "(False Position Method) ";
        clock_t StartFP = clock();
        complex <double> b = rs.findOneRoot_FalsePosition_Complex(p1, c1, c2);
        clock_t EndFP = clock();
        time2 += (EndFP - StartFP) / (double)CLOCKS_PER_SEC;
        cout << b << endl;

        cout << "(Secant Method) ";
        clock_t StartSc = clock();
        complex <double> c = rs.findOneRoot_Secant_Complex(p1, c1, c2);
        clock_t EndSc = clock();
        time3 += (EndSc - StartSc) / (double)CLOCKS_PER_SEC;
        cout << c << endl;

        cout << "(Halley Method) ";
        clock_t StartHl = clock();
        complex <double> d = rs.findOneRoot_Halley_Complex(p1, c1, c2);
        clock_t EndHl = clock();
        time4 += (EndHl - StartHl) / (double)CLOCKS_PER_SEC;
        cout << d << endl;

        cout << "(Muller Method) ";
        clock_t StartMl = clock();
        complex <double> e = rs.findOneRoot_Halley_Complex(p1, c1, c2);
        clock_t EndMl = clock();
        time5 += (EndMl - StartMl) / (double)CLOCKS_PER_SEC;
        cout << e << endl;

        cout << "(Durand-Kerner Method) " << endl;
        clock_t StartDK = clock();
        complex <double>* f = rs.DurandKerner(p1, p1.degree);
        clock_t EndDK = clock();
        time6 += (EndDK - StartDK) / (double)CLOCKS_PER_SEC;
        for (int i = 0; i < p1.degree; i++) 
        {
            cout << "     " << f[i] << endl;
        }  
        cout << "\n";

        cout << "Newton: " << time1 << " seconds\n";

        cout << "False Position: " << time2 << " seconds\n";

        cout << "Secant: " << time3 << " seconds\n";

        cout << "Halley: " << time4 << " seconds\n";

        cout << "Muller: " << time5 << " seconds\n";

        cout << "Durand-Kerner: " << time6 << " seconds\n";

        cout << "\n";
    }
    
    /* May work slowly, but feel free to check
    cout << "\n\n\n\n\n\n\n\n\n TEST 4.5\n\n";

    double test_poly5[][21] = { {-40,52,4,2,1,-40,52,4,2,1,-40,52,4,2,1,-40,52,4,2,1,20} };

    for (int i = 0; i < 1; i++)
    {
        Polynom p2 = rs.SetPoly_(20, test_poly5[i]);
        rs.printPoly();

        rs.countOfRoots();

        cout << "(Newton Method) ";
        complex <double> c1(0.2, 0.5);
        complex <double> c2(0.6, 0.9);
        complex <double> a = rs.findOneRoot_Newton_Complex(p2, c1, c2);
        cout << a << endl;

        cout << "(False Position Method) ";
        complex <double> b = rs.findOneRoot_FalsePosition_Complex(p2, c1, c2);
        cout << b << endl;

        cout << "(Secant Method) ";
        complex <double> c = rs.findOneRoot_Secant_Complex(p2, c1, c2);
        cout << c << endl;

        cout << "(Durand-Kerner Method) ";
        complex <double>* d = rs.DurandKerner(p2, p2.degree+1);
        for (int i = 0; i < p2.degree; i++)
        {
            cout << d[i] << endl;
        }
        cout << "\n";
    }
    */
    // Test 6. Complex polynomials, Jenking-Traub method WIP
    cout << "\n\n\n\n\n\n\n\n\n TEST 5\n\n";

    complex <double> c(1.0, 0.0), cc(0.0, 2.0), cc1(0.0, 1.0), cc2(1.0, 0.0);
    complex <double> carr[][4] = { {c,cc,cc1,cc2} };
    ComplexPolynom cp = rs.SetPoly_(3,carr[0]);

    rs.printCPoly();

    cout << "(Newton Method) ";
    complex <double> c4(1.0, 3.0);
    complex <double> c5(4.0, 5.0);
    complex <double> a = rs.findOneRoot_Newton_Complex(cp, c4, c5);
    cout << a << endl;

    cout << "(False Position Method) ";
    complex <double> b = rs.findOneRoot_FalsePosition_Complex(cp, c4, c5);
    cout << b << endl;

    cout << "(Secant Method) ";
    complex <double> sc = rs.findOneRoot_Secant_Complex(cp, c4, c5);
    cout << sc << endl;

    cout << "(Durand-Kerner Method) " << endl;
    complex <double>* d = rs.DurandKerner(cp, cp.degree);
    for (int i = 0; i < cp.degree; i++)
    {
        cout << "     " << d[i] << endl;
    }


    cout << "(Muller Method) ";
    complex <double> e = rs.findOneRoot_Muller_Complex(cp, c4, c5);
    cout << e << endl;

    cout << "(Halley Method) ";
    complex <double> f = rs.findOneRoot_Halley_Complex(cp, c4, c5);
    cout << f << endl;

    //cout << "(JT Method) ";
    //complex <double> g = rs.findOneRootJenkinsTraub(cp);
    //cout << g << endl;
    cout << "\n";

}
