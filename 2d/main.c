#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "plot_x11.h"
#include "image.h"

int main(void)
{
/* Выбираем, как будет задана функция: формульно или значениями в точках (из файла) */

    int file;
    printf("file? (1 = yes): ");
	scanf("%d", &file);

    char name[100];
    int n,k;
    double a,b;

    if(file != 1)
	{
		printf("Enter n: ");
		scanf("%d", &n);
		printf("Enter k: ");
		scanf("%d", &k);
		printf("Enter segment:\na = ");
		scanf("%lf", &a);
		printf("b = ");
		scanf("%lf", &b);
		if((n<k) || (0>k))
		{
            printf("Incorrect value of k (must be lower then n)");
            return -1;
		}
		if(b<=a)
		{
            printf("Incorrect value of a (must be lower then b)");
            return -1;
		}
	}
	else
	{
        printf("Enter filename: ");
        if(!scanf("%s",name))
        {
            printf("Cannot read this");
            return -1;
        }
        printf("%s",name);
	}

/* Вызываем отрисовку */

	image(file,n,k,a,b,name);

	return 1;
}








