#include <math.h>
#include <stdio.h>

double f(double x, double y);

//Функция, которую аппроксимируем
double f(double x, double y){

	return cos(sqrt(x*x+y*y));
    //return cos(x*x+y*y);
    //return x*x/2.0-y*y/2.0;
}

