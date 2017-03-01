#ifndef _WIN32

/* linux
*/

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include "xlib.h"

static int wnd_x = 0;
static int wnd_y = 0;
static int wnd_w = 800;
static int wnd_h = 600;
static int done = 0;
static char wnd_title[1024];

static Display *display = NULL;
static int ScrNum;
static Window window;
static Drawable draw;
static GC prGC;
static Font font;

typedef void (* PFRENDER) (int width, int height);
typedef void (* PFKEYBOARD) (int key);

static PFRENDER render = NULL;
static PFKEYBOARD keyboard = NULL;

void xlibInitPosition(int x, int y)
{
	wnd_x = x;
	wnd_y = y;
}

void xlibInitWindowSize(int width, int height)
{
	wnd_w = width;
	wnd_h = height;
}

void xlibRenderFunc(void (*func)(int width, int height))
{
	render = func;
}

void xlibKeyboardFunc(void (*func)(int key))
{
	keyboard = func;
}

void xlibPostRedisplay(void)
{
	XEvent event;

	event.type = Expose;
	event.xexpose.count = 0;
	XSendEvent(display, window, True, 0, &event);
}

void xlibPostExit(void)
{
	done = 1;
}

static void Exit(void)
{
	if (draw) XFreePixmap(display, draw);
	if (font) XUnloadFont(display, font);
	if (prGC) XFreeGC(display, prGC);
	if (window) XDestroyWindow(display, window);
	XFlush(display);
	if (display) XCloseDisplay(display);
}

static int translateKey(int key)
{
	int ret;

	if (key & 0xff00)
	{
		switch (key)
		{
		case XK_Escape:		ret = KEY_ESC;		break;
		case XK_Tab:		ret = KEY_TAB;		break;
		case XK_Return:		ret = KEY_RETURN;	break;
		case XK_BackSpace:	ret = KEY_BACKSPACE;	break;
		case XK_Delete:		ret = KEY_DELETE;	break;
		case XK_Home:		ret = KEY_HOME;		break;
		case XK_End:		ret = KEY_END;		break;
		case XK_Page_Up:	ret = KEY_PGUP;		break;
		case XK_Page_Down:	ret = KEY_PGDOWN;	break;
		case XK_Left:		ret = KEY_LEFT;		break;
		case XK_Right:		ret = KEY_RIGHT;	break;
		case XK_Up:		ret = KEY_UP;		break;
		case XK_Down:		ret = KEY_DOWN;		break;
		case XK_Shift_L:	ret = KEY_SHIFT;	break;
		case XK_Shift_R:	ret = KEY_SHIFT;	break;
		case XK_Control_L:	ret = KEY_CTRL;		break;
		case XK_Control_R:	ret = KEY_CTRL;		break;
		case XK_Alt_L:		ret = KEY_ALT;		break;
		case XK_Alt_R:		ret = KEY_ALT;		break;
		case XK_F1:		ret = KEY_F1;		break;
		case XK_F2:		ret = KEY_F2;		break;
		case XK_F3:		ret = KEY_F3;		break;
		case XK_F4:		ret = KEY_F4;		break;
		case XK_F5:		ret = KEY_F5;		break;
		case XK_F6:		ret = KEY_F6;		break;
		case XK_F7:		ret = KEY_F7;		break;
		case XK_F8:		ret = KEY_F8;		break;
		case XK_F9:		ret = KEY_F9;		break;
		case XK_F10:		ret = KEY_F10;		break;
		case XK_F11:		ret = KEY_F11;		break;
		case XK_F12:		ret = KEY_F12;		break;
		default: ret = (key < 0x06ff) ? (key & 0x00ff) : 0;
		}
	}
	else ret = key;

	if (ret >= 512) ret = 0;

	return ret;
}

void xlibMainLoop(const char *title)
{
	int wnd_depth;
	XEvent event;
	XWindowAttributes attr;
	KeySym key;
	Atom WM_DELETE_WINDOW;

	if ((display = XOpenDisplay(NULL)) == NULL)
	{
		printf("Cannot connect to the X server.\n");
		return;
	}

	ScrNum = DefaultScreen(display);

	window = XCreateSimpleWindow(display, RootWindow(display, ScrNum),
				    wnd_x, wnd_y, wnd_w, wnd_h, 5,
				    BlackPixel(display, ScrNum), WhitePixel(display, ScrNum));

	WM_DELETE_WINDOW = XInternAtom(display, "WM_DELETE_WINDOW", False);
	XSetWMProtocols(display, window, &WM_DELETE_WINDOW, 1);

	XSelectInput(display, window, ExposureMask | KeyPressMask);
	XFlush(display);

	strcpy(wnd_title, title);
	XStoreName(display, window, wnd_title);
	XSetIconName(display, window, wnd_title);

	XMapWindow(display, window);

	prGC = XCreateGC(display, window, 0, NULL);
	if ((font = XLoadFont(display, "*fixed*120*koi8-r"))) XSetFont(display, prGC, font);

	if (!XGetWindowAttributes(display, window, &attr))
	{
		printf("Bad window.\n");
		Exit();
		return;
	}

	wnd_w = attr.width;
	wnd_h = attr.height;
	wnd_depth = attr.depth;

	draw = XCreatePixmap(display, window, wnd_w, wnd_h, wnd_depth);
	if (!draw)
	{
		printf("Bad draw.\n");
		Exit();
		return;
	}

	while (!done)
	{
		XNextEvent(display, &event);

		switch (event.type)
		{
		case ClientMessage:
			if (event.xclient.format == 32 && event.xclient.data.l[0] == (long)WM_DELETE_WINDOW)
			{
				Exit();
				return;
			}
			break;
		case Expose:
			if (event.xexpose.count != 0) break;

			if (!XGetWindowAttributes(display, window, &attr))
			{
				printf("Bad window.\n");
				Exit();
				return;
			}

			if (wnd_w != attr.width || wnd_h != attr.height || wnd_depth != attr.depth)
			{
				XFreePixmap(display, draw);
				wnd_w = attr.width;
				wnd_h = attr.height;
				wnd_depth = attr.depth;
				draw = XCreatePixmap(display, window, wnd_w, wnd_h, wnd_depth);
				if (!draw)
				{
					printf("Bad window.\n");
					Exit();
					return;
				}
			}

			render(wnd_w, wnd_h);

			XCopyArea(display, draw, window, prGC, 0, 0, wnd_w, wnd_h, 0, 0);
			XFlush(display);
			break;
		case KeyPress:
			XLookupString(&event.xkey, NULL, 0, &key, NULL);
			key = translateKey(key);
			if (!key) break;

			keyboard(key);
			break;
		}

		while (XCheckMaskEvent(display, KeyPressMask | KeyReleaseMask, &event)) ;
	}

	Exit();
}

void xlibSetTitle(const char *title)
{
	strcpy(wnd_title, title);
	XStoreName(display, window, wnd_title);
	XSetIconName(display, window, wnd_title);
	XFlush(display);
}

void xlibSetColor(double red, double green, double blue)
{
	XColor clr;
	Colormap map = DefaultColormap(display, ScrNum);

	clr.red = red * 65535;
	clr.green = green * 65535;
	clr.blue = blue * 65535;
	clr.flags = DoRed | DoGreen | DoBlue;

	XAllocColor(display, map, &clr);
	XSetForeground(display, prGC, clr.pixel);
}

void xlibDrawString(int x, int y, const char *string, ...)
{
	char buf[4096];

	va_list argptr;
	va_start(argptr, string);
	vsnprintf(buf, sizeof(buf), string, argptr);
	va_end(argptr);
	XDrawString(display, draw, prGC, x, y, buf, strlen(buf));
}

void xlibDrawPoint(int x, int y)
{
	XDrawPoint(display, draw, prGC, x, y);
}

void xlibDrawLine(int x_start, int y_start, int x_end, int y_end)
{
	XDrawLine(display, draw, prGC, x_start, y_start, x_end, y_end);
}

void xlibDrawRectangle(int x_top_left, int y_top_left, int x_bottom_right, int y_bottom_right)
{
	XDrawRectangle(display, draw, prGC, x_top_left, y_top_left,
		    x_bottom_right - x_top_left + 1,
		    y_bottom_right - y_top_left + 1);
}

void xlibFillRectangle(int x_top_left, int y_top_left, int x_bottom_right, int y_bottom_right)
{
	XFillRectangle(display, draw, prGC, x_top_left, y_top_left,
		    x_bottom_right - x_top_left + 2,
		    y_bottom_right - y_top_left + 2);
}

void xlibDrawTriangle(int x_1, int y_1, int x_2, int y_2, int x_3, int y_3)
{
	XPoint points[3];

	points[0].x = x_1; points[0].y = y_1;
	points[1].x = x_2; points[2].y = y_2;
	points[2].x = x_3; points[2].y = y_3;
	XDrawPoints(display, draw, prGC, points, 3, CoordModeOrigin);
}

void xlibFillTriangle(int x_1, int y_1, int x_2, int y_2, int x_3, int y_3)
{
	XPoint points[3];

	points[0].x = x_1; points[0].y = y_1;
	points[1].x = x_2; points[2].y = y_2;
	points[2].x = x_3; points[2].y = y_3;
	XFillPolygon(display, draw, prGC, points, 3, Convex, CoordModeOrigin);
}

#else

#include <windows.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include "xlib.h"

/* windows
*/

static int wnd_x = 0;
static int wnd_y = 0;
static int wnd_w = 800;
static int wnd_h = 600;
static int done = 0;
static char wnd_title[1024];

static HWND window;
static HDC hdc;
static HPEN hpen;
static HBRUSH hbrush;
static double r = 0.0, g = 0.0, b = 0.0;

typedef void (* PFRENDER) (int width, int height);
typedef void (* PFKEYBOARD) (int key);

static PFRENDER render = NULL;
static PFKEYBOARD keyboard = NULL;

void xlibInitPosition(int x, int y)
{
	wnd_x = x;
	wnd_y = y;
}

void xlibInitWindowSize(int width, int height)
{
	wnd_w = width;
	wnd_h = height;
}

void xlibRenderFunc(void (*func)(int width, int height))
{
	render = func;
}

void xlibKeyboardFunc(void (*func)(int key))
{
	keyboard = func;
}

void xlibPostRedisplay(void)
{
	PostMessage(window, WM_PAINT, 0, 0);
}

void xlibPostExit(void)
{
	done = 1;
}

static void Exit(void)
{
	if (window)
	{
		ReleaseDC(window, hdc);
		DestroyWindow(window);
	}
}

static int translateKey(int key)
{
	int ret;

	switch (key)
	{
	case VK_ESCAPE:		ret = KEY_ESC;			break;
	case VK_TAB:		ret = KEY_TAB;			break;
	case VK_RETURN:		ret = KEY_RETURN;		break;
	case VK_BACK:		ret = KEY_BACKSPACE;	break;
	case VK_DELETE:		ret = KEY_DELETE;		break;
	case VK_HOME:		ret = KEY_HOME;			break;
	case VK_END:		ret = KEY_END;			break;
	case VK_PRIOR:		ret = KEY_PGUP;			break;
	case VK_NEXT:		ret = KEY_PGDOWN;		break;
	case VK_LEFT:		ret = KEY_LEFT;			break;
	case VK_RIGHT:		ret = KEY_RIGHT;		break;
	case VK_UP:			ret = KEY_UP;			break;
	case VK_DOWN:		ret = KEY_DOWN;			break;
	case VK_SHIFT:		ret = KEY_SHIFT;		break;
	case VK_MENU:		ret = KEY_ALT;			break;
	case VK_CONTROL:	ret = KEY_CTRL;			break;
	case VK_F1:			ret = KEY_F1;			break;
	case VK_F2:			ret = KEY_F2;			break;
	case VK_F3:			ret = KEY_F3;			break;
	case VK_F4:			ret = KEY_F4;			break;
	case VK_F5:			ret = KEY_F5;			break;
	case VK_F6:			ret = KEY_F6;			break;
	case VK_F7:			ret = KEY_F7;			break;
	case VK_F8:			ret = KEY_F8;			break;
	case VK_F9:			ret = KEY_F9;			break;
	case VK_F10:		ret = KEY_F10;			break;
	case VK_F11:		ret = KEY_F11;			break;
	case VK_F12:		ret = KEY_F12;			break;
	default:
		ret = MapVirtualKey(key, 2);
		if (strchr("1234567890-=", ret))
		{
/*			if (KEY_SHIFT)
			{
				if (ret == '1') ret = '!';
				else if (ret == '2') ret = '@';
				else if (ret == '3') ret = '#';
				else if (ret == '4') ret = '$';
				else if (ret == '5') ret = '%';
				else if (ret == '6') ret = '^';
				else if (ret == '7') ret = '&';
				else if (ret == '8') ret = '*';
				else if (ret == '9') ret = '(';
				else if (ret == '0') ret = ')';
				else if (ret == '-') ret = '_';
				else if (ret == '=') ret = '+';
			}
*/		}
		else if (isascii(ret))
		{
			ret = tolower(ret);
/*			if (KEY_SHIFT) ret = toupper(ret);
			else ret = tolower(ret);
*/		}
		else ret = 0;
	}

	return ret;
}

LRESULT CALLBACK windowProc(HWND window, UINT message, WPARAM wparam, LPARAM lparam)
{
	PAINTSTRUCT ps;
	HPEN hpenOld;
	HBRUSH hbrushOld;

	switch (message)
	{
	case WM_SIZE:
		wnd_w = LOWORD(lparam);
		wnd_h = HIWORD(lparam);
		return 0;
	case WM_SYSKEYDOWN:
	case WM_KEYDOWN:
		keyboard(translateKey((int)wparam));
		return 0;
	case WM_DESTROY:
		PostQuitMessage(0);
		return 0;
	case WM_PAINT:
		BeginPaint(window, &ps);

		hpen = CreatePen(PS_SOLID, 1, RGB(r, g, b));
		hpenOld = SelectObject(hdc, hpen);

		hbrush = CreateSolidBrush(RGB(r, g, b));
		hbrushOld = SelectObject(hdc, hbrush);

		render(wnd_w, wnd_h);

		SelectObject(hdc, hpenOld);
		DeleteObject(hpen);
		SelectObject(hdc, hbrushOld);
		DeleteObject(hbrush);

		EndPaint(window, &ps);
		return 0;
	}
	return DefWindowProc(window, message, wparam, lparam);
}

void xlibMainLoop(const char *title)
{
	HINSTANCE hInstance;
	WNDCLASS wc;
	RECT window_rect;
	MSG msg;

	sprintf(wnd_title, title);

	hInstance = GetModuleHandle(NULL);

	wc.style = CS_HREDRAW | CS_VREDRAW;
	wc.lpfnWndProc = (WNDPROC)windowProc;
	wc.cbClsExtra = 0;
	wc.cbWndExtra = 0;
	wc.hInstance = hInstance;
	wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);
	wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
	wc.lpszMenuName = NULL;
	wc.lpszClassName = "Xlib";

	if (!RegisterClass(&wc))
	{
		printf("RegisterClass() failed.\n");
		Exit();
		return;
	}

	window_rect.left = 0;
	window_rect.top = 0;
	window_rect.right = wnd_w;
	window_rect.bottom = wnd_h;

	AdjustWindowRectEx(&window_rect, WS_OVERLAPPEDWINDOW, 0, 0);

	wnd_w = window_rect.right - window_rect.left;
	wnd_h = window_rect.bottom - window_rect.top;

	window = CreateWindowEx(0, "Xlib", title, WS_OVERLAPPEDWINDOW, wnd_x, wnd_y,
							wnd_w, wnd_h, NULL, NULL, hInstance, NULL);

	if (!window)
	{
		printf("Cannot create window.\n");
		Exit();
		return;
	}

	hdc = GetDC(window);

	ShowWindow(window, SW_NORMAL);
	UpdateWindow(window);
	SetForegroundWindow(window);
	SetFocus(window);

	while (!done && GetMessage(&msg, window, 0, 0))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}

	Exit();
}

void xlibSetTitle(const char *title)
{
	strcpy(wnd_title, title);
	SetWindowText(window, title);
}

void xlibSetColor(double red, double green, double blue)
{
	r = red;
	g = green;
	b = blue;
	DeleteObject(hpen);
	hpen = CreatePen(PS_SOLID, 1, RGB(r * 255, g * 255, b * 255));
	SelectObject(hdc, hpen);

	DeleteObject(hbrush);
	hbrush = CreateSolidBrush(RGB(r * 255, g * 255, b * 255));
	SelectObject(hdc, hbrush);

	SetBkMode(hdc, TRANSPARENT);
	SetTextColor(hdc, RGB(r * 255, g * 255, b * 255));
}

void xlibDrawString(int x, int y, const char *string, ...)
{
	char buf[4096];

	va_list argptr;
	va_start(argptr, string);
	_vsnprintf(buf, sizeof(buf), string, argptr);
	va_end(argptr);
	TextOut(hdc, x, y, buf, (int)strlen(buf));
}

void xlibDrawPoint(int x, int y)
{
	SetPixel(hdc, x, y, RGB(r * 255, g * 255, b * 255));
}

void xlibDrawLine(int x_start, int y_start, int x_end, int y_end)
{
printf("111111111111111111111111111111111111111111111\n");
	MoveToEx(hdc, x_start, y_start, NULL);
	LineTo(hdc, x_end, y_end);
}

void xlibDrawRectangle(int x_top_left, int y_top_left, int x_bottom_right, int y_bottom_right)
{
	MoveToEx(hdc, x_top_left, y_top_left, NULL);
	LineTo(hdc, x_bottom_right, y_top_left);
	LineTo(hdc, x_bottom_right, y_bottom_right);
	LineTo(hdc, x_top_left, y_bottom_right);
	LineTo(hdc, x_top_left, y_top_left);
}

void xlibFillRectangle(int x_top_left, int y_top_left, int x_bottom_right, int y_bottom_right)
{
	Rectangle(hdc, x_top_left, y_top_left, x_bottom_right, y_bottom_right);
}

void xlibDrawTriangle(int x_1, int y_1, int x_2, int y_2, int x_3, int y_3)
{
	MoveToEx(hdc, x_1, y_1, NULL);
	LineTo(hdc, x_2, y_2);
	LineTo(hdc, x_3, y_3);
	LineTo(hdc, x_1, y_1);
}

void xlibFillTriangle(int x_1, int y_1, int x_2, int y_2, int x_3, int y_3)
{
	POINT points[3];

	points[0].x = x_1; points[0].y = y_1;
	points[1].x = x_2; points[2].y = y_2;
	points[2].x = x_3; points[2].y = y_3;
	Polyline(hdc, points, 3);
}

#endif
