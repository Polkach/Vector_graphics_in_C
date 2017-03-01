#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "math.h"

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include "plot_x11.h"

static int number = 1;//номер комбинации графиков, которая отображается
static int number2 = 4;//номер оригинального графика
static int dev = 1;//индикатор девиации
static int file = 0;//1 - если график из файла
static int n = 32;//число узлов
static int k = 16;//смещаемый узел
static double add = 0;//смещение выбранного узла
static double a = -2;//левая граница графика
static double b = 2;//правая граница графика

static double *Mf;//значения в узлах
static double *coeffs;//коэффициенты первого метода аппроксимации
static double *g;//коэффициенты первого метода аппроксимации
static double *xxx;//координаты узлов
static double *cR;//коэффициенты второго метода аппроксимации
static double *c, *cL, *r, *v;//коэффициенты второго метода аппроксимации

static double coeff_y = 1;//относительные размеры графика по оси у
static double coeff_x = 1;//относительные размеры графика по оси х

static double max = 0;
static double min = 0;
static double mx = 1e-16;
static double df1;//коэффициент для приближенного вычисления производной во втором методе
static double df2;//коэффициент для приближенного вычисления производной во втором методе

//Узлы Чебышева
double XX(int i)
{
	if(file == 1) return (double)(xxx[i]);
	if(number > 3) i = n-1-i;
	return (a+b)/2. + cos(M_PI*(2.*((double)(i)+1.)-1.)/(2.*(double)(n)))*(b-a)/2.;
}

//Равномерно распределенные узлы
double XX2(int i, int N)
{
	return a+(double)(i)*(b-a)/(double)(N-1);
}

//Стандартные функции
double F1(double x) {x = x; return 1.;}
double F2(double x) {return x;}
double F3(double x) {return x*x/2.;}
double F4(double x) {return sin(x);}
double F5(double x) {x = - x * x / 2; return exp(x);}

//Преобразователь числа в строку
void STR(char *str, int what)
{
	str[4] = str[5] = ' ';
	if(what<10) {str[0] = (what) + '0'; str[1] = ' '; str[2] = ' '; str[3] = ' ';}
	else if(what<100) {str[0] = (what/10) + '0'; str[1] = (what%10) + '0'; str[2] = ' '; str[3] = ' ';}
	else if(what<1000) {str[0] = (what/100) + '0'; str[1] = ((what/10)%10) + '0'; str[2] = (what%10) + '0'; str[3] = ' ';}
	else if(what<10000) {str[0] = (what/1000) + '0'; str[1] = ((what/100)%10) + '0'; str[2] = ((what/10)%10) + '0'; str[3] = (n%10) + '0';}
	else {str[0] = str[1] = str[2] = str[3] = '9'; str[4] = '+';}
}

//Вычисление коэффициентов первого метода аппроксимации
int estimation(void)
{
    int i, j;
    double z, g1 = 0., g2 = 0.;
    for(i=0; i<n; i++) coeffs[i] = g[i] = 0.;

    for(j=1; j<=n; j++)
    {
        z = 2.*cos(M_PI*(2.*(double)(j)-1.)/(2.*(double)(n)));
        for(i=0; i<n; i++)
        {
            if(i==0)
            {
                g[i] = Mf[j-1];
                g2 = g[i];
            }
            else
            {
                if(i==1)
                {
                    g[i] = z*Mf[j-1]/2;
                    g1 = g[i];
                }
                else
                {
                    g[i] = z*g1 - g2;
                    g2 = g1;
                    g1 = g[i];
                }
            }
        }
        for(i=0; i<n; i++) coeffs[i] += g[i];
    }

    coeffs[0] = coeffs[0]/(double)(n);
    for(i=1; i<n; i++) coeffs[i] = 2.*coeffs[i]/(double)(n);
    return 0;
}

//Аппрокисмирующая функция первым методом
double Pf(double x)
{
    double z;
    double P = 0.;
    double T0 = 1., T1 = 0., T = 0.;
    int i;
    z = 2.*(2.*x-(b+a))/(b-a); if(file==1) z = 2.*(2.*x-(xxx[n-1]+xxx[0]))/(xxx[n-1]-xxx[0]);
    T1 = z/2.;
    P = coeffs[0]*T0 + coeffs[1]*T1;
    for(i=2; i<n; i++)
    {
        T = z*T1 - T0;
        P += coeffs[i]*T;
        T0 = T1;
        T1 = T;
    }
    return P;
}

//Оригинальный график
double F(double x)
{
	if(number2 == 1) return F1(x);
	if(number2 == 2) return F2(x);
	if(number2 == 3) return F3(x);
	if(number2 == 4) return F4(x);
	return F5(x);
}

//Производная оригинального графика
double dF(double x)
{
	if(number2 == 1) return 0.;
	if(number2 == 2) return F1(x);
	if(number2 == 3) return F2(x);
	if(number2 == 4) return cos(x);
	return (-x)*F5(x);
}

//Кусочно-линейная функция с заданными в точках значениями
double oF(double x)
{
	int i;
	for(i=0; i<n-1; i++)
	{
		if(x < xxx[i] || x > xxx[i+1]) continue;
		else return Mf[i] + ((x-xxx[i])/(xxx[i+1]-xxx[i]))*(Mf[i+1]-Mf[i]);
	}
	printf("fault!\n");
	return 0;
}

//Середины отрезков между узлами Чебвышева
double ksi(int i)
{
    if(i==0) return a - (ksi(1)-a);
    else if(i==n) return b + (b-ksi(n-1));
    else return (XX(i-1) + XX(i))/2;
}

//Вспомогательная процедура для второго метода аппрокисмации
int solve(void)
{
    int i;
    cR[0] /= c[0]; r[0] /= c[0]; c[0] = 1.;
    for(i=1; i<=n; i++)
    {
        c[i] -= cL[i]*cR[i-1];
        r[i] -= cL[i]*r[i-1];
        //cL[i] = 0.;
        cR[i] /= c[i];
        r[i] /= c[i];
        //c[i] = 1.;
    }
    v[n] = r[n];
    for(i=n-1; i>=0; i--) v[i] = r[i] - v[i+1]*cR[i];
    return 0;
}

//Вычисление коэффициентов первого метода аппроксимации
int estimation2(void)
{
    int i;
    if(file != 1) {df1 = dF(XX(0)); df2 = dF(XX(n-1));}
    r[0] = df1 - Mf[0]*( (1/(XX(0)-ksi(0))) - (1/(ksi(1)-XX(0))) );
    c[0] = (1/(ksi(1)-ksi(0))) - (1/(XX(0)-ksi(0)));
    cR[0] = (1/(ksi(1)-XX(0))) - (1/(ksi(1)-ksi(0)));
    r[n] = df2 - Mf[n-1]*( (1/(XX(n-1)-ksi(n-1))) - (1/(ksi(n)-XX(n-1))) );
    cL[n] = (1/(ksi(n)-ksi(n-1))) - (1/(XX(n-1)-ksi(n-1)));
    c[n] = (1/(ksi(n)-XX(n-1))) - (1/(ksi(n)-ksi(n-1)));
    for(i=1; i<n; i++)
    {
        cL[i] = (1/(XX(i-1)-ksi(i-1))) - (1/(ksi(i)-ksi(i-1)));
         c[i] = (1/(ksi(i)-XX(i-1))) + (1/(ksi(i)-ksi(i-1))) + (1/(XX(i)-ksi(i))) + (1/(ksi(i+1)-ksi(i)));
        cR[i] = (1/(ksi(i+1)-XX(i))) - (1/(ksi(i+1)-ksi(i)));
         r[i] = Mf[i-1]*( (1/(XX(i-1)-ksi(i-1))) + (1/(ksi(i)-XX(i-1))) ) + Mf[i]*( (1/(XX(i)-ksi(i))) + (1/(ksi(i+1)-XX(i))) );
    }
    solve();
    for(i=0; i<n; i++)
    {
        cL[i] = v[i];
        cR[i] = (1/(ksi(i+1)-ksi(i))) * ( ((v[i+1]-Mf[i])/(ksi(i+1)-XX(i))) - ((Mf[i]-v[i])/(XX(i)-ksi(i))) );
         c[i] = (Mf[i] - v[i])/(XX(i)-ksi(i)) - (XX(i)-ksi(i))*cR[i];
    }
    return 0;
}

//Аппрокисмирующая функция первым методом
double Pf2(double x)
{
    int i;
    for(i=0; i<n; i++)
    {
        if(x<ksi(i) || x>ksi(i+1)) continue;
        else
        {
            return cL[i] + c[i]*(x-ksi(i)) + cR[i]*(x-ksi(i))*(x-ksi(i));
        }
    }
	return 0;
}

//Функция девиации
double Deviation(double x)
{
	if(file == 1 && number < 4) return fabs(oF(x)-Pf(x))/mx;
	if(file == 1 && number > 3) return fabs(oF(x)-Pf2(x))/mx;
	if(number < 4) return fabs(F(x)-Pf(x))/mx;
	else return fabs(F(x)-Pf2(x))/mx;
}

//Отображение одного графика на экране
void Show(double (*fun)(double))
{
  	int w = nWWidth, h = nWHeight, i;
	int x_start, y_start, x_end, y_end;
	int flag;
	double flag2;
	flag2 = fabs(fun(XX2(0, w)));
	max = min = fun(XX2(0, w));
	for(i = 0; i < w; i++)
	{
		if(max < fun(XX2(i, w))) max = fun(XX2(i, w));
		if(min > fun(XX2(i, w))) min = fun(XX2(i, w));
	} if(max-min<1e-4) {min--; max++;}
  	for (i = 0; i < w-1; i++)
 	{
   		x_start = (int)( coeff_x *w*(XX2(i, w) - a)/ (b-a) - w*(coeff_x-1)/2 );
   		y_start = (int)( coeff_y *h*(fabs(fun(XX2(i, w))-max))/ (fabs(max-min)) - h*(coeff_y-1)/2 );
  		x_end = (int)( coeff_x *w*(XX2(i+1, w) - a)/ (b-a) - w*(coeff_x-1)/2  );
  		y_end = (int)( coeff_y *h*(fabs(fun(XX2(i+1, w))-max))/ (fabs(max-min)) - h*(coeff_y-1)/2 );
   		WDrawLine (x_start, y_start, x_end, y_end);
		if(fabs(fun(XX2(i, w))) < flag2)
		{
			flag = y_start; flag2 = fabs(fun(XX2(i, w)));
		}
  	}
	WSetColor (BLACK);
	if(fabs(fabs(a)-fabs(b))<1e-15 && a<b && coeff_x >= 1./2. && coeff_x <= 2. && coeff_y >= 1./2. && coeff_y <= 2. )
	{
		WDrawLine(w/2, 0, w/2, h);
		WDrawLine(0, flag, w, flag);
	}
	//printf("show1()a=%lf,b=%lf\n",a,b);
}

//Отображение двух графиков на экране
void Show2(double (*fun)(double), double (*fun2)(double))
{
  	int w = nWWidth, h = nWHeight, i;
	int x_start, y_start, x_end, y_end;
	int flag;
	double flag2;
	flag2 = fabs(fun(XX2(0, w)));
	max = min = fun(XX2(0, w));

	for(i = 0; i < w; i++)
	{
		if(max < fun(XX2(i, w))) max = fun(XX2(i, w));
		if(max < fun2(XX2(i, w))) max = fun2(XX2(i, w));
		if(min > fun(XX2(i, w))) min = fun(XX2(i, w));
		if(min > fun2(XX2(i, w))) min = fun2(XX2(i, w));
	} if(max-min<1e-4) {min--; max++;}
  	for (i = 0; i < w-1; i++)
 	{
		WSetColor (GREEN);
   		x_start = (int)( coeff_x *w*(XX2(i, w) - a)/ (b-a) - w*(coeff_x-1)/2 );
   		y_start = (int)( coeff_y *h*(fabs(fun(XX2(i, w))-max))/ (fabs(max-min)) - h*(coeff_y-1)/2 );
  		x_end = (int)( coeff_x *w*(XX2(i+1, w) - a)/ (b-a) - w*(coeff_x-1)/2  );
  		y_end = (int)( coeff_y *h*(fabs(fun(XX2(i+1, w))-max))/ (fabs(max-min)) - h*(coeff_y-1)/2 );
   		WDrawLine (x_start, y_start, x_end, y_end);
		if(fabs(fun(XX2(i, w))) < flag2)
		{
			flag = y_start; flag2 = fabs(fun(XX2(i, w)));
		}
		WSetColor (RED);
		x_start = (int)( coeff_x *w*(XX2(i, w) - a)/ (b-a) - w*(coeff_x-1)/2 );
   		y_start = (int)( coeff_y *h*(fabs(fun2(XX2(i, w))-max))/ (fabs(max-min)) - h*(coeff_y-1)/2 );
  		x_end = (int)( coeff_x *w*(XX2(i+1, w) - a)/ (b-a) - w*(coeff_x-1)/2  );
  		y_end = (int)( coeff_y *h*(fabs(fun2(XX2(i+1, w))-max))/ (fabs(max-min)) - h*(coeff_y-1)/2 );
   		WDrawLine (x_start, y_start, x_end, y_end);
		if(fabs(fun2(XX2(i, w))) < flag2)
		{
			flag = y_start; flag2 = fabs(fun2(XX2(i, w)));
		}
  	}
	WSetColor (BLACK);
	if(fabs(fabs(a)-fabs(b))<1e-15 && a<b && coeff_x >= 1./2. && coeff_x <= 2. && coeff_y >= 1./2. && coeff_y <= 2. )
	{
		WDrawLine(w/2, 0, w/2, h);
		WDrawLine(0, flag, w, flag);
	}
}

//Окно программы
void DrawWindowContent(void)
{
  double x;
  int w = nWWidth;
  int h = nWHeight;
  int i;
  char str[3], stri[8]="00000000";
	double A, B;

//Титл окна и цвета оформления
  if(number2 == 1) WSetTitle("f(x) = 1");
  if(number2 == 2) WSetTitle("f(x) = x");
  if(number2 == 3) WSetTitle("f(x) = x*x/2");
  if(number2 == 4) WSetTitle("f(x) = sin(x)");
  if(number2 == 5) WSetTitle("f(x) = e^(-x*x)");
  if(file == 1) WSetTitle("FILE");
  WSetColor (DARKGRAY);
  WFillRectangle (0, 0, w, h);

  WSetColor (BLACK);

//Пояснительные надписи
  WDrawString ("Q=quit, F1..F4 -- change scale, F5/F6 -- change node count, G -- change graphic, D -- show deviation, F11/F12 -/+ k", 10, 20);
  WDrawString ("N= ", 10, 60); STR(str, n); WDrawString (str, 22, 60);

  WDrawString ("x= ", 10, 80);
  if(coeff_x>1) {STR(str, (int)(coeff_x)); WDrawString (str, 22, 80);}
  else {WDrawString ("1/", 22, 80); STR(str, (int)(1/coeff_x)); WDrawString (str, 34, 80);}
  WDrawString ("y= ", 10+50, 80);
  if(coeff_y>1) {STR(str, (int)(coeff_y)); WDrawString (str, 22+50, 80);}
  else {WDrawString ("1/", 22+50, 80); STR(str, (int)(1/coeff_y)); WDrawString (str, 34+50, 80);}

	WSetColor (BLACK);
	if(number == 1) WDrawString("Original function", 10, 40);
	if(number == 2) WDrawString("First method", 10, 40);
	if(number == 3) WDrawString("Original + first method", 10, 40);
	if(number == 4) WDrawString("Second method", 10, 40);
	if(number == 5) WDrawString("Original + Second method", 10, 40);


  if(file == 1)
  {
	number2 = 2;

	if(number == 1) {WSetColor (GREEN); Show(oF);}
	else {
	if(number < 4)
	{
    		estimation();
		if(number == 2) {WSetColor (RED); Show(Pf);}
		if(number == 3) {Show2(oF, Pf);}
		for(i = 0; i < n; i++) if( mx < fabs( Mf[i]-Pf(xxx[i]) ) ) mx = fabs( Mf[i]-Pf(xxx[i]) );
	}
	else
	{
		estimation2();
		if(number == 4) {WSetColor (RED); Show(Pf2);}
		if(number == 5) {Show2(oF, Pf2);}
	      for(i = 0; i < w; i++) if( mx < fabs( oF(XX2(i, w))-Pf2(XX2(i, w)) ) ) mx = fabs( oF(XX2(i, w))-Pf2(XX2(i, w)));
	}
	if(dev==2) {WSetColor (BROWN); Show(Deviation);}
	printf("mx == %e\n", mx);
	for(i=0; i<8; i++)
	{
		stri[i] = ' ';
		mx = mx*10;
		stri[i] = ((int)(mx)%10) + '0';
	}
	WSetColor (BLACK); WDrawString ("deviation= 0.", 22+50, 60); WDrawString (stri, 22+130, 60);
	mx = 1e-16;}
  } else

  if(number == 1) {WSetColor (GREEN); Show(F);}
  else if(number == 2 || number == 3)
  {
	Mf = (double *)malloc(n * sizeof(double));
    	for(i=0; i<n; i++) {Mf[i] = F(XX(i)); if(i==k) Mf[i] += add;}
     	g = (double *)malloc(n * sizeof(double));
    	coeffs = (double *)malloc(n * sizeof(double));
    	estimation();
	if(number == 3) Show2(F, Pf);
	if(number == 2) {WSetColor (RED); Show(Pf);}
	for(i = 0; i < w; i++) if( mx < fabs( F(XX2(i, w))-Pf(XX2(i, w)) ) ) mx = fabs( F(XX2(i, w))-Pf(XX2(i, w)) );
	if(dev==2) {WSetColor (BROWN); Show(Deviation);}
	for(i=0; i<8; i++)
	{
		stri[i] = ' ';
		mx = mx*10;
		stri[i] = ((int)(mx)%10) + '0';
	}
	WSetColor (BLACK); WDrawString ("deviation= 0.", 22+50, 60); WDrawString (stri, 22+130, 60);

	mx = 1e-16;
	free(Mf); free(g); free(coeffs);
	printf("free mf,g,coeffs\n");
  }
  else
  {
	Mf = (double *)malloc(n * sizeof(double));
    	for(i=0; i<n; i++) {Mf[i] = F(XX(i)); if(i==k) Mf[i] += add;}
	cR = (double *)malloc((n) * sizeof(double));
    	c = (double *)malloc((n+1) * sizeof(double));
    	cL = (double *)malloc((n+1) * sizeof(double));
    	r = (double *)malloc((n+1) * sizeof(double));
    	v = (double *)malloc((n+1) * sizeof(double));
	estimation2();
	if(number == 5) Show2(F, Pf2);
	if(number == 4) {WSetColor(RED); Show(Pf2);}

	mx = 1e-16;
	for(i = 0; i < w; i++) if( mx < fabs( F(XX2(i, w))-Pf2(XX2(i, w)) ) ) {mx = fabs( F(XX2(i, w))-Pf2(XX2(i, w)) );}
	if(dev==2) {WSetColor (BROWN); Show(Deviation);}
	printf("mx == %e\n", mx);
	for(i=0; i<8; i++)
	{
		stri[i] = ' ';
		mx = mx*10;
		stri[i] = ((int)(mx)%10) + '0';
	}
	WSetColor (BLACK); WDrawString ("deviation= 0.", 22+50, 60); WDrawString (stri, 22+130, 60);
	mx = 1e-16;
	free(Mf);
	free(cR); free(c); free(cL); free(r); free(v);
  }

}

//Реакции на различные нажатия клавиш
int KeyPressFunction(int nKeySym)
{
  switch (nKeySym)
  {
  case XK_Q:
  case XK_q:
  case XK_Escape:
    free(Mf); free(xxx); free(coeffs); free(g); free(cR); free(c); free(cL); free(r); free(v);
    return KEY_PRESS_QUIT;

  case XK_Return:
    if(file != 1) n = 8;
    k = n/2;
    add = 0.;
    coeff_x = 1.0;
    coeff_y = 1.0;
    break;

//Меняет горизонтальные размеры графика
  case XK_F1:
    if(coeff_x > 1./64.) coeff_x *= 0.5;
    break;
  case XK_F2:
    if(coeff_x < 1024.) coeff_x *= 2.0;
    break;

//Меняет вертикальные размеры графика
  case XK_F3:
    if(coeff_y > 1./64.) coeff_y *= 0.5;
    break;
  case XK_F4:
    if(coeff_y < 1024.) coeff_y *= 2.0;
    break;

//Меняет количество узлов аппроксимации
  case XK_F6:
    if(file != 1) n *=2;
    break;
  case XK_F5:
    if(file != 1) if (n > 4) n /= 2;
    if(n%2==1) n++;
    break;

//Меняет оригинальный график
  case XK_G:
  case XK_g:
    number2 = (number2+1)%5; if(number2 == 0) number2 = 5;
    break;

//Показывает девиацию (нормированный модуль разности)
  case XK_D:
  case XK_d:
    dev = (dev+1)%2; if(dev == 0) dev = 2;
    break;

//Оригинальный график
  case XK_1:
    number = 1;
    break;

//График, построенный первым методом
  case XK_2:
    number = 2;
    break;

//Оригинальный график и график, построенный первым методом
  case XK_3:
    number = 3;
    break;

//График, построенный втормы методом
  case XK_4:
    number = 4;
    break;

//Оригинальный график и график, построенный вторым методом
  case XK_5:
    number = 5;
    break;

//Смещения выбранного узла
  case XK_F12:
    add = add+1;
    break;
  case XK_F11:
    add = add-1;
    break;

  default:
    return KEY_PRESS_NOTHING;
  }

  return KEY_PRESS_EXPOSE;
}

//Заполняет глобальные переменные и создает окно программы
int image(int file2, int n2, int k2, double a2, double b2, char* name2)
{
    file = file2;
	if(file2==1)
	{
        FILE * f=fopen(name2,"r");
        if (f == NULL)
        {
            printf("Cannot find file");
            return -1;
        }
        else
        {
            fscanf(f, "%d", &n);
            xxx = (double *)malloc(n * sizeof(double));
            Mf = (double *)malloc(n * sizeof(double));
            int i;
            for(i=0; i<n; i++)
            {
                fscanf(f, "%lf", &xxx[i]);
                fscanf(f, "%lf", &Mf[i]);
            }
            a=xxx[0]; b=xxx[n-1];
            fscanf(f, "%lf", &df1);
            fscanf(f, "%lf", &df2);
            fclose(f);
            k=0;
        }
	}
	else
	{
        n = n2;
        if(n%2==1) n++;
        k = k2;
        a = a2;
        b = b2;
        xxx = (double *)malloc(n * sizeof(double));
		Mf = (double *)malloc(n * sizeof(double));
	}

    g = (double *)malloc(n * sizeof(double));
	coeffs = (double *)malloc(n * sizeof(double));
	cR = (double *)malloc(n * sizeof(double));
	c = (double *)malloc((n+1) * sizeof(double));
	cL = (double *)malloc((n+1) * sizeof(double));
	r = (double *)malloc((n+1) * sizeof(double));
	v = (double *)malloc((n+1) * sizeof(double));

	printf("Retcode: %d", DrawWindow(DrawWindowContent, KeyPressFunction));

	return 1;
}

