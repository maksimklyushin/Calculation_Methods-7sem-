#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

vector<vector<double>> A()
{
    vector<vector<double>> a;
    a.resize(3);
    for(int i=0; i<3; i++)
        a[i].resize(3);
    a[0][0] = 7.1056*pow(10, -6);
    a[0][1] = 6.5056*pow(10, -3);
    a[0][2] = 0.93056;
    a[1][0] = -6.8888*pow(10, -3);
    a[1][1] = -0.68888;
    a[1][2] = 0.91112;
    a[2][0] = 4.41168;
    a[2][1] = 1.64168;
    a[2][2] = 1.09168;
    return a;
}

vector<double> B()
{
    vector<double> b;
    b.resize(3);
    b[0] = 4.09877;
    b[1] = 1.83706;
    b[2] = 2.20580;
    return b;
}

vector<vector<double>> AB()
{
    vector<vector<double>> a;
    a.resize(4);
    for(int i=0; i<4; i++)
        a[i].resize(3);
    a[0][0] = 7.1056*pow(10, -6);
    a[0][1] = 6.5056*pow(10, -3);
    a[0][2] = 0.93056;
    a[1][0] = -6.8888*pow(10, -3);
    a[1][1] = -0.68888;
    a[1][2] = 0.91112;
    a[2][0] = 4.41168;
    a[2][1] = 1.64168;
    a[2][2] = 1.09168;
    a[3][0] = 4.09877;
    a[3][1] = 1.83706;
    a[3][2] = 2.20580;
    return a;
}

void print_nevyazka(vector<double> x)
{
    vector<vector<double>> a = AB();
    double err = 0;
    cout << endl;
    cout << "Error:"<< endl;
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            err = err + a[j][i]*x[j];
        }
        if(i==1)
        {
            cout << "|A*X - B| =   " << setw(10) << abs(err - a[3][i]) << endl;
        }
        else
            {
            cout << "              "<< setw(10) << abs(err - a[3][i]) << endl;
            }
        err = 0;
    }
}

void printmatrix3(vector<vector<double>> a)
{
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<4; j++)
        {
            cout << setw(10) << a[j][i] << "     ";
        }
        cout << endl;
    }

}

void printsolution(vector<double> a)
{
    cout << endl;
    cout << "Solution:" << endl;
    cout << "    ";
    printf("%.15f    ", a[0]);
    cout << endl;
    cout << "X = ";
    printf("%.15f    ", a[1]);
    cout << endl;
    cout << "    ";
    printf("%.15f    ", a[2]);
    cout << endl;
}

void det()
{
    vector<vector<double>> a = A();
    double d;
    d = a[0][0]*a[1][1]*a[2][2] + a[1][0]*a[2][1]*a[0][2] + a[2][0]*a[0][1]*a[1][2] - a[0][2]*a[1][1]*a[2][0] - a[1][2]*a[2][1]*a[0][0] - a[2][2]*a[0][1]*a[1][0];
    cout << endl;
    cout << "Determinant: det(A) = " << d << endl;
}

vector<vector<double>> swap1 (vector<vector<double>> &a, int m, int n)
{
    double p;
    for(int i=0; i<4; i++)
    {
        p = a[i][m];
        a[i][m] = a[i][n];
        a[i][n] = p;
    }
    return a;
}

vector<vector<double>> swap2 (vector<vector<double>> &a, int m, int n)
{
    double p;
    for(int i=0; i<6; i++)
    {
        p = a[i][m];
        a[i][m] = a[i][n];
        a[i][n] = p;
    }
    return a;
}

vector<vector<double>> inverse_matrix(vector<vector<double>> a)
{
    vector<vector<double>> aa;
    vector<vector<double>> bb;
    double m, tmp;
    int q;
    aa.resize(6);
    for(int i=0; i<6; i++)
        aa[i].resize(3);
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
        {
            aa[j][i] = a[j][i];
            aa[i+3][i] = 1;
        }
    for(int k=0; k<3; k++)
    {
        q = k;
        m = abs(aa[k][k]);
        for(int p=k+1; p<3; p++)
        {
            if(abs(aa[k][p])>m)
            {
                m = abs(aa[k][p]);
                q = p;
            }
        }
        if(q!=k)
        {
        aa = swap2(aa, q, k);
        }
        tmp = aa[k][k];
        for(int j=k; j<6; j++)
        {
         aa[j][k] = aa[j][k]/tmp;
        }
        for(int i=k+1; i<3; i++)
        {
            tmp = aa[k][i];
            for(int j=k; j<6; j++)
            {
                aa[j][i] = aa[j][i] - aa[j][k]*tmp;
            }
        }
    }
    for(int i=2; i>0; i--)
    {
        for(int j=i-1; j>=0; j--)
        {
            tmp = aa[i][j];
            for(int k=i; k<6; k++)
            {
                aa[k][j] = aa[k][j] - aa[k][i]*tmp;
            }
        }
    }
    bb.resize(3);
    for(int i=0; i<3; i++)
        bb[i].resize(3);
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
            bb[j][i] = aa[j+3][i];
    }
    cout << "Inverse matrix, A^-1:" << endl;
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
            cout << setw(10) << bb[j][i] << "       ";
        cout << endl;
    }
    cout << endl;
    return bb;
}

void checking2(vector<vector<double>> a, vector<vector<double>> aa)
{
    vector<vector<double>> bb;
    double s=0;
    bb.resize(3);
    for(int i=0; i<3; i++)
        bb[i].resize(3);
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {   for(int k=0; k<3; k++)
                s = s + a[k][i]*aa[j][k];
            bb[j][i] = s;
            s = 0;
        }

    }
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
            cout << setw(10) << bb[j][i] << "       ";
        cout << endl;
    }
}

void Gauss1(double epsilon)
{
    cout << "I.GAUSS METHOD" << endl;
    cout << endl;
    vector<vector<double>> ab = AB();
    vector<double> x;
    x.resize(3);
    double tmp;
    cout << endl;
    cout << "Initial matrix:" << endl;
    printmatrix3(ab);
    cout << endl;
    cout << endl;
    cout << "Pryamoi hod:" << endl;
    for(int k=0; k<3; k++)
    {
        cout << endl;
        cout << "Step " << k+1 << ":" << endl;
        tmp = ab[k][k];
        if(tmp <= epsilon)
        {
            cout << "Warning! The first element <= " << epsilon << " !" << endl;
        }
        for(int j=k; j<4; j++)
        {
         ab[j][k] = ab[j][k]/tmp;
        }
        for(int i=k+1; i<3; i++)
        {
            tmp = ab[k][i];
            for(int j=k; j<4; j++)
            {
                ab[j][i] = ab[j][i] - ab[j][k]*tmp;
            }

        }
         printmatrix3(ab);
         cout << endl;
    }
    ab[3][1] = ab[3][1] - ab[2][1]*ab[3][2];
    ab[2][1] = 0;
    ab[3][0] = ab[3][0] - ab[2][0]*ab[3][2] - ab[1][0]*ab[3][1];
    ab[2][0] = 0;
    ab[1][0] = 0;
    cout << endl;
    cout << "Obratnyi hod:" << endl;
    printmatrix3(ab);
    cout << endl;
    for(int i=0; i<3; i++)
        x[i] = ab[3][i];
    printsolution(x);
    cout << endl;
    print_nevyazka(x);
    cout << endl;
    det();
    cout << endl;
}

void Gauss2()
{
    vector<vector<double>> ab = AB();
    vector<vector<double>> a = A();
    vector<vector<double>> aa;
    vector<double> x;
    double m;
    double tmp;
    int q;
    x.resize(3);
    cout << endl;
    cout << "II. GAUSS with choosing the column element" << endl;
    cout << endl;
    cout << "Initial matrix:" << endl;
    printmatrix3(ab);
    cout << endl;
    cout << endl;
    cout << "Pryamoi hod:" << endl;
    cout << endl;
    for(int k=0; k<3; k++)
    {
        q = k;
        m = abs(ab[k][k]);
        for(int p=k+1; p<3; p++)
        {
            if(abs(ab[k][p])>m)
            {
                m = abs(ab[k][p]);
                q = p;
            }
        }
        if(q!=k)
        {
        cout << "Swap:" << endl;
        ab = swap1(ab, q, k);
        printmatrix3(ab);
        cout << endl;
        }
        cout << endl;
        cout << "Step " << k+1 << ":" << endl;
        tmp = ab[k][k];
        for(int j=k; j<4; j++)
        {
         ab[j][k] = ab[j][k]/tmp;
        }
        for(int i=k+1; i<3; i++)
        {
            tmp = ab[k][i];
            for(int j=k; j<4; j++)
            {
                ab[j][i] = ab[j][i] - ab[j][k]*tmp;
            }

        }
         printmatrix3(ab);
         cout << endl;
    }
    ab[3][1] = ab[3][1] - ab[2][1]*ab[3][2];
    ab[2][1] = 0;
    ab[3][0] = ab[3][0] - ab[2][0]*ab[3][2] - ab[1][0]*ab[3][1];
    ab[2][0] = 0;
    ab[1][0] = 0;
    cout << endl;
    cout << endl;
    cout << "Obratnyi hod:" << endl;
    printmatrix3(ab);
    cout << endl;
    for(int i=0; i<3; i++)
        x[i] = ab[3][i];
    printsolution(x);
    cout << endl;
    print_nevyazka(x);
    cout << endl;
    det();
    cout << endl;
    cout << endl;
    aa = inverse_matrix(a);
    cout << endl;
    cout << endl;
    cout << "Checking: A*A^-1 = " << endl;
    checking2(a, aa);
    cout << endl;
}


void LU()
{
    vector<vector<double>> a = A();
    vector<vector<double>> L, U, LB, UB;
    vector<double> y, x, b;
    y.resize(3);
    x.resize(3);
    double s,ss, tmp;
    L.resize(3);
    U.resize(3);
    LB.resize(4);
    UB.resize(4);
    for(int i=0; i<4; i++)
    {
        LB[i].resize(3);
        UB[i].resize(3);
    }
    for(int i=0; i<3; i++)
    {
        L[i].resize(3);
        U[i].resize(3);
    }
    L[0][0] = a[0][0];
    for(int i=0; i<3; i++)
    {
        U[i][i] = 1;
        for(int j=i; j<3; j++)
        {
            s = 0;
            ss = 0;
            for(int k=0; k<=i-1; k++)
            {
                s = s + L[k][j]*U[i][k];
                ss = ss + L[k][i]*U[j][k];
            }
            L[i][j] = a[i][j] - s;
            U[j][i] = (a[j][i] - ss)/L[i][i];
        }
    }
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "III. LU method." << endl;
    cout << endl;
    cout << endl;
    cout << "Matrix L:" << endl;
    cout << endl;
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
            cout << setw(10) << L[j][i] << "       ";
        cout << endl;
    }
    cout << endl;
    cout << endl;
    cout << "Matrix U:" << endl;
    cout << endl;
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
            cout << setw(10) << U[j][i] << "       ";
        cout << endl;
    }
    cout << endl;
    cout << endl;
    b = B();
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            LB[j][i] = L[j][i];
            UB[j][i] = U[j][i];
        }
    }
    for(int i=0; i<3; i++)
    {
        LB[3][i] = b[i];
    }
    for(int k=0; k<3; k++)
    {
        tmp = LB[k][k];
        for(int j=k; j<4; j++)
        {
         LB[j][k] = LB[j][k]/tmp;
        }
        for(int i=k+1; i<3; i++)
        {
            tmp = LB[k][i];
            for(int j=k; j<4; j++)
            {
                LB[j][i] = LB[j][i] - LB[j][k]*tmp;
            }

        }
    }
    cout << "Result of Gauss for system Ly = b:" << endl;
    printmatrix3(LB);
    cout << endl;
    cout << endl;
    for(int i=0; i<3; i++)
        y[i] = LB[3][i];

    for(int i=0; i<3; i++)
    {
        UB[3][i] = y[i];
    }
    for(int i=2; i>0; i--)
    {
        for(int j=i-1; j>=0; j--)
        {
            tmp = UB[i][j];
            for(int k=i; k<4; k++)
            {
                UB[k][j] = UB[k][j] - UB[k][i]*tmp;
            }
        }
    }
    for(int i=0; i<3; i++)
        x[i] = UB[3][i];
    cout << endl;
    cout << "Result of Gauss for system Ux = y:" << endl;
    printmatrix3(UB);
    cout << endl;
    printsolution(x);
    cout << endl;
    print_nevyazka(x);
}

void LU2()
{
    vector<vector<double>> a = A();
    vector<vector<double>> L, U, LB, UB;
    vector<vector<double>> ab = AB();
    vector<double> y, x, b;
    y.resize(3);
    x.resize(3);
    double s,ss, tmp, m;
    int q;
    L.resize(3);
    U.resize(3);
    LB.resize(4);
    UB.resize(4);
    for(int i=0; i<4; i++)
    {
        LB[i].resize(3);
        UB[i].resize(3);
    }
    for(int i=0; i<3; i++)
    {
        L[i].resize(3);
        U[i].resize(3);
    }
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "IV. LU method (improved)." << endl;
    cout << endl;
    cout << endl;
    cout << "Initial matrix:" << endl;
    cout << endl;
    printmatrix3(ab);
    cout << endl;
    for(int i=0; i<3; i++)
    {
        q = i;
        m = abs(ab[i][i]);
        for(int p=i+1; p<3; p++)
        {
            if(abs(ab[i][p])>m)
            {
                m = abs(ab[i][p]);
                q = p;
            }
        }
        if(q!=i)
        {
        cout << "Swap:" << endl;
        ab = swap1(ab, q, i);
        printmatrix3(ab);
        cout << endl;
        }
        L[0][0] = ab[0][0];
        U[i][i] = 1;
        for(int j=i; j<3; j++)
        {
            s = 0;
            ss = 0;
            for(int k=0; k<=i-1; k++)
            {
                s = s + L[k][j]*U[i][k];
                ss = ss + L[k][i]*U[j][k];
            }
            L[i][j] = ab[i][j] - s;
            U[j][i] = (ab[j][i] - ss)/L[i][i];
        }
    }
    cout << endl;
    cout << endl;
    cout << "Matrix L:" << endl;
    cout << endl;
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
            cout << setw(10) << L[j][i] << "       ";
        cout << endl;
    }
    cout << endl;
    cout << endl;
    cout << "Matrix U:" << endl;
    cout << endl;
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
            cout << setw(10) << U[j][i] << "       ";
        cout << endl;
    }
    cout << endl;
    cout << endl;
    b = B();
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            LB[j][i] = L[j][i];
            UB[j][i] = U[j][i];
        }
    }
    for(int i=0; i<3; i++)
    {
        LB[3][i] = ab[3][i];
    }
    for(int k=0; k<3; k++)
    {
        tmp = LB[k][k];
        for(int j=k; j<4; j++)
        {
         LB[j][k] = LB[j][k]/tmp;
        }
        for(int i=k+1; i<3; i++)
        {
            tmp = LB[k][i];
            for(int j=k; j<4; j++)
            {
                LB[j][i] = LB[j][i] - LB[j][k]*tmp;
            }

        }
    }
    cout << "Result of Gauss for system Ly = b:" << endl;
    cout << endl;
    printmatrix3(LB);
    cout << endl;
    cout << endl;
    for(int i=0; i<3; i++)
        y[i] = LB[3][i];

    for(int i=0; i<3; i++)
    {
        UB[3][i] = y[i];
    }
    for(int i=2; i>0; i--)
    {
        for(int j=i-1; j>=0; j--)
        {
            tmp = UB[i][j];
            for(int k=i; k<4; k++)
            {
                UB[k][j] = UB[k][j] - UB[k][i]*tmp;
            }
        }
    }
    for(int i=0; i<3; i++)
        x[i] = UB[3][i];
    cout << endl;
    cout << "Result of Gauss for system Ux = y:" << endl;
    cout << endl;
    printmatrix3(UB);
    cout << endl;
    printsolution(x);
    cout << endl;
    print_nevyazka(x);
}


int main()
{
    system("color f0");
    Gauss1(0.0001);
    cout << endl;
    Gauss2();
    cout << endl;
    LU();
    cout << endl;
    LU2();
    return 0;
}
