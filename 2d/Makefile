all:
	gcc -c -O3 -ffast-math -W -Wall -I. -I/usr/X11R6/include plot_x11.c
	gcc -c image.c
	gcc -c main.c
	gcc plot_x11.o image.o main.o -L. -L/usr/X11R6/lib -lX11 -lm -o go
