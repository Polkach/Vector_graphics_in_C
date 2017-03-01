#include <stdio.h>
#include "math.h"
#include "plot_x11.h"

/*static*/ void DrawWindowContent(void);
/*static*/ int KeyPressFunction(int nKeySym);

double XX(int i);
double XX2(int i, int N);
double F1(double x);
double F2(double x);
double F3(double x);
double F4(double x);
double F5(double x);
void STR(char *str, int what);
int estimation(void);
double Pf(double x);
double F(double x);
double dF(double x);
double oF(double x);
double ksi(int i);
int solve(void);
int estimation2(void);
double Pf2(double x);
double Deviation (double x);
void Show(double (*fun)(double));
void Show2(double (*fun)(double), double (*fun2)(double));

int image(int file2, int n2, int k2, double a2, double b2, char* name2);
