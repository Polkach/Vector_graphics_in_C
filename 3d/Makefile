NAME            = plot

CC              = gcc
LINK            = gcc

DEBUG           = 
OPTIMAZE_COMMON   = -O3 -ffast-math
OPTIMAZE_SPECIFIC =
OPTIMAZE        = $(OPTIMAZE_COMMON) $(OPTIMAZE_SPECIFIC)

WARNINGS        = -W -Wall
INCLUDE_DIR     = -I. -I/usr/X11R6/include
CFLAGS          = -c $(DEBUG) $(OPTIMAZE) $(WARNINGS) $(INCLUDE_DIR)

LIB_DIR         = -L. -L/usr/X11R6/lib
LIB             = -lX11
LDFLAGS         = $(DEBUG)

OBJS            = main.o metod1.o xlib.o functions.o	

all:    $(NAME)

$(NAME)         : $(OBJS)
	$(LINK) $(LDFLAGS) $(OBJS) $(LIB_DIR) $(LIB) -lm -o $(NAME)

clean:
	rm -f $(OBJS) $(NAME)

.c.o:
	$(CC) $(CFLAGS) $<

main.o      : main.c main.h metod1.h xlib.h functions.h

metod1.o          : metod1.c metod1.h functions.h

xlib.o		: xlib.c xlib.h

functions.o	: functions.h functions.c
