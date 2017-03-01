
double f(double x, double y);

double dxf(double x, double y);

double dyf(double x, double y);

double ddf(double x, double y);

double func1(double x, double y);

int Init(void);

void Finalize(void);

typedef struct Point{
	int x;
	int y;
} Point;

void drawLine(Point *p1, Point *p2);

void DrawFunction(int width, int height, double (*func)(double x, double y));

void CalcMaxAndMin(int width, int height, double (*func)(double x, double y),double* maxim, double* minim);

void render(int width, int height);

void keyboard(int key);