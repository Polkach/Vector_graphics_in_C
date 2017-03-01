//Кубическое приближение двумерной функции

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "main.h"
#include "metod1.h"
#include "functions.h"

extern double pogr1;
static int n;
static double a;
static double b;

static int m;
static double c;
static double d;

static double *xarr;
static double *yarr;
static double **f_xy;

static double **D_x;
static double **D_yt;
static double **G_x;
static double **F;
static double **Gamma;


int Init_1(int n, double a, double b, int m, double c, double d);

void Finalize_1(void);

void Input_1(void);

void Delta_1(int number1, int number2, double delta);

void Calc_1(void);

double Polinom1(double x, double y);

void Coeff_1(int n, int m);

double Value_1(double x, double y, int n, int m, double *xarr, double *yarr);

double Pogreshnost1(double x, double y);

double FinePogreshnost1(double x, double y);

double *allocation(int size);

void ShowM(double **a, int n, int m);

//Выделение память под массивы, необходимые для аппроксимации
int Init_1(int n_, double a_, double b_, int m_, double c_, double d_){
	int i;

	n = n_;
	a = a_;
	b = b_;

	m = m_;
	c = c_;
	d = d_;

    xarr=(double*)malloc(n*sizeof(double));
    yarr=(double*)malloc(m*sizeof(double));
    f_xy = malloc(sizeof(double*)*(n));
    for(i=0;i<n;i++) f_xy[i] = malloc(sizeof(double)*(m));


    D_x = malloc(sizeof(double*)*(n));
    for(i=0;i<n;i++) D_x[i] = malloc(sizeof(double)*(m));
    D_yt = malloc(sizeof(double*)*(n));
    for(i=0;i<n;i++)D_yt[i] = malloc(sizeof(double)*(m));
    G_x = malloc(sizeof(double*)*(n));
    for(i=0;i<n;i++) G_x[i] = malloc(sizeof(double)*(n));
    F = malloc(sizeof(double*)*(n));
    for(i=0;i<n;i++) F[i] = malloc(sizeof(double)*(m));
    Gamma = malloc(sizeof(double*)*((n-1)*4));
    for(i=0;i<(n-1)*4;i++) Gamma[i] = malloc(sizeof(double)*((m-1)*4));


	if (!(xarr && yarr && f_xy && D_x && D_yt && G_x && F && Gamma))
        return 0;

	return 1;
}

//Очистка памяти
void Finalize_1(void){
    int i;
	free(xarr);
	free(yarr);
	for(i=0;i<n;i++)free(F[i]);
    free(F);
    for(i=0;i<(n-1)*4;i++)free(Gamma[i]);
    free(Gamma);
    for(i=0;i<n;i++)free(D_x[i]);
    free(D_x);
    for(i=0;i<n;i++)free(D_yt[i]);
    free(D_yt);
    for(i=0;i<n;i++)free(G_x[i]);
    free(G_x);
}

//Равномерная расстановка узлов и вычисление значений функции в этих точках
void Input_1(void){
	int i, j;
	double h;
    int p=10;

        printf("Ravnomern1\n");

        h=(b-a)/(n-1.0);

        for(i=0; i<n; i++){
            xarr[i]=a+h*i;
        }

        h=(d-c)/(m-1.0);

        for(i=0; i<m; i++){
            yarr[i]=c+h*i;
        }
        if(n<p)
            p=n;

        printf("xarr: ");
        for(i=0; i<p; i++){
            printf("%-10.6lf ",xarr[i]);
        }
        printf("\n");
        printf("yarr: ");

        if(m<p) p=m;

        for(i=0; i<p; i++){
            printf("%-10.6lf ",yarr[i]);
        }
        printf("\n");



    for(i=0; i<n; i++)
        for(j=0; j<m; j++){
            f_xy[i][j]=f(xarr[i],yarr[j]);
        }
    printf("f_xy\n");
}

//Вызов вычисления коэффициентов
void Calc_1(void){

	Coeff_1(n, m);
}

//Вычисление коэффициентов для аппркосимирующей функции
void Coeff_1(int n, int m){

    int i,j,k,p,q;
    double hx = xarr[1]-xarr[0];
    double hy = yarr[1]-yarr[0];
    double **A;
    double **B;
    double **C;
    A = malloc(sizeof(double*)*(4));
    for(i=0;i<4;i++) A[i] = malloc(sizeof(double)*(4));
    B = malloc(sizeof(double*)*(4));
    for(i=0;i<4;i++) B[i] = malloc(sizeof(double)*(4));
    C = malloc(sizeof(double*)*(4));
    for(i=0;i<4;i++) C[i] = malloc(sizeof(double)*(4));

    for(i=1;i<n-1;i++)
        for(j=0;j<m;j++)
            D_x[i][j] = (f(xarr[i+1],yarr[j])-f(xarr[i-1],yarr[j]))/2.0/hx;
    for(j=0;j<m;j++)
    {
        D_x[0][j] = (3.0*(f(xarr[1],yarr[j])-f(xarr[0],yarr[j]))/hx-D_x[1][j])/2.0;
        D_x[n-1][j] = (3.0*(f(xarr[n-1],yarr[j])-f(xarr[n-2],yarr[j]))/hx-D_x[n-2][j])/2.0;
    }

    printf("D_x\n");

    for(j=1;j<m-1;j++)
        for(i=0;i<n;i++)
            D_yt[i][j] = (f(xarr[i],yarr[j+1])-f(xarr[i],yarr[j-1]))/hy/2.0;
    for(i=0;i<n;i++)
    {
        D_yt[i][0] = (3.0*(f(xarr[i],yarr[1])-f(xarr[i],yarr[0]))/hy-D_yt[i][1])/2.0;
        D_yt[i][m-1] = (3.0*(f(xarr[i],yarr[m-1])-f(xarr[i],yarr[m-2]))/hy-D_yt[i][m-2])/2.0;
    }

    printf("D_yt\n");

    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            G_x[i][j] = 0.0;
    for(i=1;i<n-1;i++)
    {
        G_x[i][i-1] = -1.0/2.0/hx;
        G_x[i][i+1] = 1.0/2.0/hx;
    }
    G_x[0][0] = -5.0/4.0/hx;
    G_x[0][1] = 3.0/2.0/hx;
    G_x[0][2] = -1.0/4.0/hx;
    G_x[n-1][n-1] = 5.0/4.0/hx;
    G_x[n-1][n-2] = -3.0/2.0/hx;
    G_x[n-1][n-3] = 1.0/4.0/hx;


    printf("G_x\n");

    for(i=0;i<n;i++)
        for(j=0;j<m;j++)
        {
            F[i][j] = 0.0;
            for(k=0;k<n;k++)
                F[i][j] += G_x[i][k]*D_yt[k][j];
        }

    printf("F\n");

    for(i=0;i<n-1;i++)
        for(j=0;j<m-1;j++)
        {
            Gamma[i*4+0][j*4+0] = f_xy[i][j];
            Gamma[i*4+0][j*4+1] = D_yt[i][j];
            Gamma[i*4+0][j*4+2] = f_xy[i][j+1];
            Gamma[i*4+0][j*4+3] = D_yt[i][j+1];
            Gamma[i*4+1][j*4+0] = D_x[i][j];
            Gamma[i*4+1][j*4+1] = F[i][j];
            Gamma[i*4+1][j*4+2] = D_x[i][j+1];
            Gamma[i*4+1][j*4+3] = F[i][j+1];
            Gamma[i*4+2][j*4+0] = f_xy[i+1][j];
            Gamma[i*4+2][j*4+1] = D_yt[i+1][j];
            Gamma[i*4+2][j*4+2] = f_xy[i+1][j+1];
            Gamma[i*4+2][j*4+3] = D_yt[i+1][j+1];
            Gamma[i*4+3][j*4+0] = D_x[i+1][j];
            Gamma[i*4+3][j*4+1] = F[i+1][j];
            Gamma[i*4+3][j*4+2] = D_x[i+1][j+1];
            Gamma[i*4+3][j*4+3] = F[i+1][j+1];
        }

    printf("Gamma\n");

    A[0][0] = 1.0;
    A[0][1] = 0.0;
    A[0][2] = 0.0;
    A[0][3] = 0.0;
    A[1][0] = 0.0;
    A[1][1] = 1.0;
    A[1][2] = 0.0;
    A[1][3] = 0.0;
    A[2][0] = -3.0/hx/hx;
    A[2][1] = -2.0/hx;
    A[2][2] = 3.0/hx/hx;
    A[2][3] = -1.0/hx;
    A[3][0] = 2.0/hx/hx/hx;
    A[3][1] = 1.0/hx/hx;
    A[3][2] = -2.0/hx/hx/hx;
    A[3][3] = 1.0/hx/hx;

    B[0][0] = 1.0;
    B[1][0] = 0.0;
    B[2][0] = 0.0;
    B[3][0] = 0.0;
    B[0][1] = 0.0;
    B[1][1] = 1.0;
    B[2][1] = 0.0;
    B[3][1] = 0.0;
    B[0][2] = -3.0/hy/hy;
    B[1][2] = -2.0/hy;
    B[2][2] = 3.0/hy/hy;
    B[3][2] = -1.0/hy;
    B[0][3] = 2.0/hy/hy/hy;
    B[1][3] = 1.0/hy/hy;
    B[2][3] = -2.0/hy/hy/hy;
    B[3][3] = 1.0/hy/hy;
    printf("A\n");
    printf("B\n");

    for(i=0;i<n-1;i++)
        for(j=0;j<m-1;j++)
        {
            for(p=0;p<4;p++)
                for(q=0;q<4;q++)
                {
                    C[p][q] = 0.0;
                    for(k=0;k<4;k++)
                        C[p][q] += A[p][k]*Gamma[i*4+k][j*4+q];
                }

            for(p=0;p<4;p++)
                for(q=0;q<4;q++)
                {
                    Gamma[i*4+p][j*4+q] = 0.0;
                    for(k=0;k<4;k++)
                        Gamma[i*4+p][j*4+q] += C[p][k]*B[k][q];
                }
        }

    printf("Gamma\n");

    for(i=0;i<4;i++)free(A[i]);
    free(A);
    for(i=0;i<4;i++)free(B[i]);
    free(B);
    for(i=0;i<4;i++)free(C[i]);
    free(C);
}

//Аппроксимирующий полином
double Polinom1(double x, double y){

	return Value_1(x, y, n, m, xarr, yarr);
}

//Вычисление значений аппроксимирующей функции с помощью вычисленных коэффициентов
double Value_1(double x, double y, int n, int m, double *xarr, double *yarr){

	int i,j,k,l;
	double t;
	for(i=1;i<n;i++)
		for(j=1;j<m;j++)
		{
			if(x<xarr[i] && y<yarr[j])
			{
				t = 0.0;
				for(k=0;k<4;k++)
					for(l=0;l<4;l++)
					{
						t += Gamma[(i-1)*4+k][(j-1)*4+l]*pow(x-xarr[i-1],k)*pow(y-yarr[j-1],l);
					}
				return t;
			}
		}
	for(i=1;i<n;i++)
		if(x<xarr[i])
		{
			t = 0.0;
			for(k=0;k<4;k++)
				for(l=0;l<4;l++)
				{
					t += Gamma[(i-1)*4+k][(m-2)*4+l]*pow(x-xarr[i-1],k)*pow(y-yarr[m-2],l);

				}
			return t;

		}

	for(j=1;j<m;j++)
		if(y<yarr[j])
		{
			t = 0.0;
			for(k=0;k<4;k++)
				for(l=0;l<4;l++)
				{
					t += Gamma[(n-2)*4+k][(j-1)*4+l]*pow(x-xarr[n-2],k)*pow(y-yarr[j-1],l);
				}
			return t;

		}
	return f(xarr[n-1],yarr[m-1]);
}

//Дельта функция
void Delta_1(int number1, int number2, double delta){

	f_xy[number1][number2] += delta;
}

//Погрешность
double Pogreshnost1(double x, double y){

	return f(x, y) - Polinom1(x, y);
}

//Нормированная погрешность
double FinePogreshnost1(double x, double y){

    if(  (x < xarr[0]) || (x > xarr[n-1]) || (y < yarr[0]) || (y > yarr[m-1]) )return 0;

    return (1.0/pogr1)*Pogreshnost1(x,y);

}

//Выводит на экран пересечение первых n строк и m столбцов
void ShowM(double **a, int n, int m)
{
    int i, j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<m;j++)
            printf("%f    ",a[i][j]);
        printf("\n");
    }
}


