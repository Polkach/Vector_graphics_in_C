#ifndef __X_LIB_H__
#define __X_LIB_H__

void xlibInitPosition(int x, int y);

void xlibInitWindowSize(int width, int height);

void xlibRenderFunc(void (*func)(int width, int height));

void xlibKeyboardFunc(void (*func)(int key));

void xlibPostRedisplay(void);

void xlibPostExit(void);

void xlibMainLoop(const char *title);

void xlibSetTitle(const char *title);

/* Rendering functions */

void xlibSetColor(double red, double green, double blue);

void xlibDrawString(int x, int y, const char *string, ...);

void xlibDrawPoint(int x, int y);

void xlibDrawLine(int x_start, int y_start, int x_end, int y_end);

void xlibDrawRectangle(int x_top_left, int y_top_left, int x_bottom_right, int y_bottom_right);

void xlibFillRectangle(int x_top_left, int y_top_left, int x_bottom_right, int y_bottom_right);

void xlibDrawTriangle(int x_1, int y_1, int x_2, int y_2, int x_3, int y_3);

void xlibFillTriangle(int x_1, int y_1, int x_2, int y_2, int x_3, int y_3);

/* Keys definition */

#define KEY_ESC			256
#define KEY_TAB			257
#define KEY_RETURN		258
#define KEY_BACKSPACE		259
#define KEY_DELETE		260
#define KEY_HOME		261
#define KEY_END			262
#define KEY_PGUP		263
#define KEY_PGDOWN		264
#define KEY_LEFT		265
#define KEY_RIGHT		266
#define KEY_UP			267
#define KEY_DOWN		268
#define KEY_SHIFT		269
#define KEY_CTRL		270
#define KEY_ALT			271
#define KEY_F1			272
#define KEY_F2			273
#define KEY_F3			274
#define KEY_F4			275
#define KEY_F5			276
#define KEY_F6			277
#define KEY_F7			278
#define KEY_F8			279
#define KEY_F9			280
#define KEY_F10			281
#define KEY_F11			282
#define KEY_F12			283

#endif /* __X_LIB_H__ */
