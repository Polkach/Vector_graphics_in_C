#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "main.h"
#include "xlib.h"

#include "metod1.h"
#include "functions.h"

static int n;
static double a;
static double b;

static int m;
static double c;
static double d ;

static double max, min;

static double x_left, x_right;
static double y_up, y_down;

static int variant = 1;

static int num_steps = 80;
static double phi = 30.0 * 3.1075926/180.0;
static double psi = 10.0 * 3.1075926/180.0;

double pogr1=0.0;
double pogr2=0.0;

int YMax[1000];
int YMin[1000];

#define NO_VALUE  7777

//Отрисовка линии (проектирование на поверхность экрана)
void drawLine(Point *p1, Point *p2)
{
	int dx = fabs(p2->x - p1->x);
	int dy = fabs(p2->y - p1->y);
	int sx = p2->x >= p1->x ? 1 : -1;
	int sy = p2->y >= p1->y ? 1 : -1;
	int d;
	int d1;
	int d2;
	int x;
	int y;
	int i;
	int m1;
	int m2;

	if (dy <= dx)
	{
		d = -dx;
		d1 = dy << 1;
		d2 = (dy - dx) << 1;

		for (x = p1->x, y = p1->y, i = 0; i <= dx; i++, x += sx)
		{
			if (YMin[x] == NO_VALUE)
			{
				xlibDrawPoint(x, y);
				YMin[x] = YMax[x] = y;
			}
			else if (y < YMin[x])
			{
				xlibDrawPoint(x, y);
				YMin[x] = y;
			}
			else if (y > YMax[x])
			{
				xlibDrawPoint(x, y);
				YMax[x] = y;
			}

			if (d > 0)
			{
				d += d2;
				y += sy;
			}
			else d += d1;
		}
	}
	else
	{
		d = -dy;
		d1 = dx << 1;
		d2 = (dx - dy) << 1;
		m1 = YMin[p1->x];
		m2 = YMax[p1->x];

		for (x = p1->x, y = p1->y, i = 0; i <= dy; i++, y += sy)
		{
			if (YMin[x]  == NO_VALUE)
			{
				xlibDrawPoint(x, y);
				YMin[x] = YMax[x] = y;
			}
			else if (y < m1)
			{
				xlibDrawPoint(x, y);
				if (y < YMin[x]) YMin[x] = y;
			}
			else if (y > m2)
			{
				xlibDrawPoint(x, y);
				if (y > YMax[x]) YMax[x] = y;
			}

			if (d > 0)
			{
				d += d2;
				x += sx;
				m1 = YMin[x];
				m2 = YMax[x];
			}
			else d += d1;
		}
	}
}

//Отрисовка функции (проектирование трехмерного объекта на поверхность экрана)
void DrawFunction(int width, int height, double (*func)(double x, double y))
{
	Point* curLine = (Point*)malloc(num_steps * sizeof(Point));
	Point* nextLine = (Point*)malloc(num_steps * sizeof(Point));
	double sphi = sin(phi);
	double cphi = cos(phi);
	double spsi = sin(psi);
	double cpsi = cos(psi);
	double e1[3] = {cphi, sphi, 0.0};
	double e2[3] = {sphi * spsi, -spsi * cphi, cpsi};
	double x, y, z;
	double hx = (x_right - x_left)/(num_steps - 1);
	double hy = (y_up - y_down)/(num_steps - 1);
	double xMin = (e1[0] >= 0.0 ? x_left : x_right) * e1[0] + (e1[1] >= 0.0 ? y_down : y_up) * e1[1];
	double xMax = (e1[0] >= 0.0 ? x_right : x_left) * e1[0] + (e1[1] >= 0.0 ? y_up : y_down) * e1[1];
	double yMin = (e2[0] >= 0.0 ? x_left : x_right) * e2[0] + (e2[1] >= 0.0 ? y_down : y_up) * e2[1];
	double yMax = (e2[0] >= 0.0 ? x_right : x_left) * e2[0] + (e2[1] >= 0.0 ? y_up : y_down) * e2[1];
	double ax;
	double bx;
	double ay;
	double by;
	int i, j;
    double freemy;

    int p1_x,p1_y,p2_x,p2_y;

    freemy=0.0*(width+height);


	if (e2[2] >= 0.0)
	{
		yMin += min * e2[2];
		yMax += max * e2[2];
	}
	else
	{
		yMin += max * e2[2];
		yMax += min * e2[2];
	}

	ax = -750 * xMin/(xMax - xMin);
	bx = 750/(xMax - xMin);
	ay = -550 * yMin/(yMax - yMin);
	by = 550/(yMax - yMin);

	for (i = 0; i < (int)(sizeof(YMax)/sizeof(int)); i++)
		YMin[i] = YMax[i] = NO_VALUE;

	for (i = 0; i < num_steps; i++)
	{
		x = x_left + i * hx;
		y = y_down + (num_steps - 1) * hy;
		z = func(x, y);
		curLine[i].x = (int)(ax + bx * (x * e1[0] + y * e1[1]));
		curLine[i].y = height - (int)(ay + by * (x * e2[0] + y * e2[1] + z * e2[2]));
	}
	for (i = num_steps - 1; i > -1; i--)
	{
		for (j = 0; j < num_steps - 1; j++) drawLine(&curLine[j], &curLine[j + 1]);
		if (i > 0)
			for (j = 0; j < num_steps; j++)
			{
				x = x_left + j * hx;
				y = y_down + (i - 1) * hy;
				z = func(x, y);
				nextLine[j].x = (int)(ax + bx * (x * e1[0] + y * e1[1]));
				nextLine[j].y = height - (int)(ay + by * (x * e2[0] + y * e2[1] + z * e2[2]));
				drawLine(&curLine[j], &nextLine[j]);
				curLine[j] = nextLine[j];
			}
	}

	free(curLine);
	free(nextLine);
	width = freemy/height;
	xlibSetColor(0.0, 0.0, 0.0);
	xlibDrawString(10, 10,"1 - change graph");
	xlibDrawString(10, 20,"2 - zoom in");
	xlibDrawString(10, 30,"3 - zoom out");
	xlibDrawString(10, 40,"4 - points /=2");
	xlibDrawString(10, 50,"5 - points *=2");
	xlibDrawString(10, 60,"6 - delta -=1");
	xlibDrawString(10, 70,"7 - delta +=1");
	xlibDrawString(10, 80,"8 - turn left");
	xlibDrawString(10, 90,"9 - turn right");
	xlibDrawString(10, 100,"- - high detalization");
	xlibDrawString(10, 110,"= - low detalization");
	xlibDrawString(10, 120,"q - exit");

    p1_x = (int)(ax+bx*((x_left+x_right)/2*e1[0]+y_down*e1[1]));
	p1_y = height-(int)(ay+by*((x_left+x_right)/2*e2[0]+y_down*e2[1]));
	p2_x = (int)(ax+bx*((x_left+x_right)/2*e1[0]+y_up*e1[1]));
	p2_y = height-(int)(ay+by*((x_left+x_right)/2*e2[0]+y_up*e2[1]));
	xlibDrawLine(p1_x,p1_y,p2_x,p2_y);
    xlibDrawString(p2_x+10, p2_y+10,"y");
    p1_x = (int)(ax+bx*(x_left*e1[0]+(y_down+y_up)/2*e1[1]));
	p1_y = height-(int)(ay+by*(x_left*e2[0]+(y_down+y_up)/2*e2[1]));
	p2_x = (int)(ax+bx*(x_right*e1[0]+(y_down+y_up)/2*e1[1]));
	p2_y = height-(int)(ay+by*(x_right*e2[0]+(y_down+y_up)/2*e2[1]));
	xlibDrawLine(p1_x,p1_y,p2_x,p2_y);
    xlibDrawString(p2_x+10, p2_y+10,"x");
    	if (variant != 3){
        xlibDrawString(10,130,"max of function = %.2f  min of function = %.2f\n", max, min);
        p1_x = (int)(ax+bx*((x_left+x_right)/2*e1[0]+(y_down+y_up)/2*e1[1]));
        p1_y = height-(int)(ay+by*((x_left+x_right)/2*e2[0]+(y_down+y_up)/2*e2[1]+(int)(min)*e2[2]));
        p2_x = (int)(ax+bx*((x_left+x_right)/2*e1[0]+(y_down+y_up)/2*e1[1]));
        p2_y = height-(int)(ay+by*((x_left+x_right)/2*e2[0]+(y_down+y_up)/2*e2[1]+(int)(max)*e2[2]));
        xlibDrawLine(p1_x,p1_y,p2_x,p2_y);
        xlibDrawString(p2_x+10, p2_y+10,"z");
	}
}

//Вычисление максимума и минимум оображаемого сегмента
void CalcMaxAndMin(int width, int height, double (*func)(double x, double y), double* maxim, double* minim)
{
	int i, j;
	double x, y, z;
	double hx, hy;
    double freemy;

	hx = (x_right - x_left)/num_steps;//(width - 1);
	hy = (y_up - y_down)/num_steps;//(height - 1);

	*maxim = *minim = func(x_left, y_down);
	for (i = 0; i < num_steps; i++)
		for (j = 0; j < num_steps; j++)
		{
			x = x_left + i * hx;
			y = y_down + j * hy;
			z = func(x, y);
			if (z > *maxim) *maxim = z;
			if (z < *minim) *minim = z;
		}
    freemy=0.0*(width+height);
	width = freemy/height;
}



void render(int width, int height){
	double r1, r2;

	xlibSetColor(1.0, 1.0, 1.0);
	xlibFillRectangle(0, 0, width, height);

    printf("\nn: %d  \n", n);
    printf("m: %d  \n", m);

//Титл окна
	switch (variant)
	{
	case 1:
		printf("Real Function\n");
		xlibSetTitle("Original");
		break;
	case 2:
		printf("Interplation method1\n");
		xlibSetTitle("Interpolation");
		break;
	case 3:
		printf("Pogreshnost1\n");
		xlibSetTitle("Error");
		break;
	}

	xlibSetColor(0.0, 0.0, 0.0);

//Информация о видимом сегменте
	if (variant != 3){

		CalcMaxAndMin(width, height, f, &max, &min);
		printf("max of function = %.2f  min of function = %.2f\n", max, min);
	}
    if(variant==2 || variant==3){
        CalcMaxAndMin(width, height, Pogreshnost1, &r1, &r2);
        if (fabs(r1) < fabs(r2)) r1 = r2;
        pogr1 = fabs(r1);
        printf("Pogreshnost1 = %le\n", pogr1);
    }

//Отрисовка самого графика
	switch (variant){
        case 1:
            xlibSetColor(1.0, 0.0, 0.0);
            DrawFunction(width, height, f);
            break;
        case 2:
            xlibSetColor(0.0, 0.0, 1.0);
            DrawFunction(width, height, Polinom1);
            break;
        case 3:
            xlibSetColor(0.0, 1.0, 0.0);
            DrawFunction(width, height, FinePogreshnost1);
            break;
    }
}

//Отвечает за действия, привязанные к клавишам
void keyboard(int key){
	double tmp;
	int number1;
	int number2;
	double delta;

	switch (key)
	{
	//Выход
	case 'q':
	case 'Q':
	case KEY_ESC:
		xlibPostExit();
		return;

    //Смена графика
	case '1':
        printf("\nXK_1\n");
		if (variant == 1) variant = 2;
		else
		{if (variant == 2) variant = 3; else{
		if (variant == 3) variant = 1;}}
		if (!variant) break;
		if (variant == 2 || variant == 3)
		{
			Input_1();
			Calc_1();
		}
		break;

    //Приближение и отдаление
	case '2':
	case '3':
		tmp = (x_right - x_left);
        if (key == '2'){
            printf("\nXK_2 down size:\n");
            tmp /= 4.0;
        }
        else {
            printf("\nXK_3 up sixe:\n");
            tmp /= -2.0;
        }
		x_left += tmp;
		x_right -= tmp;
		tmp = (y_up - y_down);
		if (key == '2') tmp /= 4.0;
		else tmp /= -2.0;
		y_down += tmp;
		y_up -= tmp;
		break;

    //Число узлов аппроксимации
	case '4':
	case '5':
            printf("Before finalizw\n");
		Finalize();
            printf("Afterfinalize\n");

		if (key == '4'  && n > 10 && m > 10){
			n /= 2;
			m /= 2;
		}
		else{
			n *= 2;
			m *= 2;
		}

		if (!Init_1(n, a, b, m, c, d)){

			printf("Not enough memory.\n");
			Finalize();
			xlibPostExit();
			return;
		}

		if (variant == 2 || variant == 3){

			Input_1();
			Calc_1();
		}
		break;

    //Смещение одного значения
	case '6':
	case '7':
		number1 = n/3;
		number2 = m/3;
        if (key == '6'){
            printf("\nXK_6 delta = -1.0:\n");
            delta = -1.0;
		}
        else{
            printf("\nXK_7 delta = +1.0:\n");
            delta = +1.0;
		}

        if (variant == 2 || variant == 3){

			Delta_1(number1, number2, delta);
			Calc_1();
		}
		break;

    //Поворот графика
	case '8':
        printf("\nXK_8 phi-= :\n");
		phi -= 12*(3.1415926/180.0)/4.0;
		break;
	case '9':
        printf("\nXK_9 phi+= :\n");
		phi += 12*(3.1415926/180.0)/4.0;
		break;

    //Детализация
	case '-':
        printf("\nXK_- num_steps *= 2; :\n");
		num_steps *= 2;
		break;
	case '=':
        printf("\nXK_= num_steps /= 2 :\n");
		if (num_steps > 5) num_steps /= 2;
		break;

    //Смещение вдоль осей
    case KEY_LEFT:
		x_left -= 1.0;
		x_right -= 1.0;
		break;
	case KEY_RIGHT:
		x_left += 1.0;
		x_right += 1.0;
		break;
	case KEY_DOWN:
		y_down -= 1.0;
		y_up -= 1.0;
		break;
	case KEY_UP:
		y_down += 1.0;
		y_up += 1.0;
		break;
	default:
		return;
	}

	xlibPostRedisplay();
}

int Init(void){
    int zn;

//Считываем с клавиатуры требуемые диапазоны
    printf("Input from: 1=input.txt  and  0=function\n");
    if(!(scanf("%d",&zn))) return -1;
    if(zn==0){
        printf("Please, enter a and b: ");
        if(!(scanf("%lg %lg", &a, &b))) return -1;
        //a=-3; b=3;

        printf("Now, please, enter n: ");
        if(!(scanf("%d", &n))) return -1;
        //n=4;

        printf("Please, enter c and d: ");
        if(!(scanf("%lg %lg", &c, &d))) return -1;
        //c=-1.5; d=1.5;

        printf("Now, please, enter m: ");
        if(!(scanf("%d", &m))) return -1;
        //m=4;

        if (b <= a || n < 3 || d <= c || m < 3){
            printf("Some input error.\n");
            return 0;
        }

//Задаем размеры
        x_left = a;
        x_right = b;
        y_down = c;
        y_up = d;

//Аппроксимируем
        if (!Init_1(n, a, b, m, c, d)){
            printf("Not enough memory.\n");
            Finalize();
            return 0;
        }

    }
    if(zn==1){
        printf("input.txt\n");


    }

	return 1;
}

void Finalize(void){

	Finalize_1();
}

int main(void){
	int width = 800;
	int height = 600;

//Инициируем окно программы
	xlibInitPosition(0, 0);
	xlibInitWindowSize(width, height);
	xlibRenderFunc(render);
	xlibKeyboardFunc(keyboard);

	if (!Init()) return -1;

	xlibMainLoop("Interpolation");

//Чистим память
	Finalize();

	return 0;
}
