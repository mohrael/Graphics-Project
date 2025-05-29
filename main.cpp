#include <windows.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>
#include <commdlg.h>
#include <locale>
#include <codecvt>
#include <climits>
#include <stack>

using namespace std;

// Menu command IDs
#define CLEAR 101
#define BLACK 102
#define RED 103
#define GREEN 104
#define BLUE 105
#define YELLOW 106
#define PINK 107
#define DDA_LINE 201
#define MIDPOINT_LINE 202
#define PARAMETRIC_LINE 203
#define CARTESIAN_CIRCLE 301
#define POLAR_CIRCLE 302
#define ITERPOLAR_CIRCLE 303
#define BRES_CIRCLE 304
#define MODBRES_CIRCLE 305
#define SAVE 401
#define LOAD 402
#define FILL_RECURSIVE  501
#define FILL_NONRECURSIVE  502
#define FILL_CANCEL   503


enum ShapeType { LINE ,CIRCLE , FILL};
enum LineAlgorithm { LINE_DDA, LINE_MIDPOINT, LINE_PARAMETRIC };
enum CircleAlgorithm { CIRCLE_CARTESIAN, CIRCLE_POLAR ,CIRCLE_ITERPOLAR,CIRCLE_BRES , CIRCLE_MODBRES ,FILL_CIRCLE_QUARTER};
enum FillAlgorithm { FLOOD_FILL_RECURSIVE, FLOOD_FILL_NONRECURSIVE };


struct Shape {
    ShapeType type;
    int x1, y1, x2, y2;
    COLORREF color;
    LineAlgorithm linealgo;
    CircleAlgorithm circleAlgo;
    COLORREF fillColor;
    COLORREF currColor;
    FillAlgorithm fillAlgo;
    int quarter = 1;
};
vector<Shape> shapes;
bool waitingForSecondClick = false;

//line algorithms
//drawing line DDA

void DDAline(HDC hdc , int x1 , int y1 ,int x2, int y2 , COLORREF c){
    int dx = x2 -x1;
    int dy = y2 -y1;
    SetPixel(hdc , x1 , y1 , c);
    if (abs(dx)>=abs(dy)){
        int x = x1;
        int xinc = dx>0?1:-1;
        double y=y1;
        double yinc=(double)dy/dx*xinc;
        while(x!= x2){
            x+=xinc;
            y+=yinc;
            SetPixel(hdc , x ,round(y) ,c);
        }
    }
    else{
        int y = y1;
        int yinc = dy>0?1:-1;
        double x=x1;
        double xinc=(double)dx/dy *yinc;
        while(y!= y2){
            x+=xinc;
            y+=yinc;
            SetPixel(hdc , round(x) ,y ,c);

        }
    }
}
//line berseham
void LineBersenham(HDC hdc , int x1 , int y1 , int x2 , int y2 , COLORREF c) {
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int x = x1;
    int y = y1;
    int xinc = (x2 > x1) ? 1 : -1;
    int yinc = (y2 > y1) ? 1 : -1;
    SetPixel(hdc, x, y, c);
    if (dx >= dy) {
        int d = 2 * dy - dx;
        int d1 = 2 * dy;
        int d2 = 2 * (dy - dx);
        while (x != x2) {
            if (d < 0) {
                d += d1;
            } else {
                y += yinc;
                d += d2;
            }
            x += xinc;
            SetPixel(hdc, x, y, c);
        }
    }
    else {
        int d = 2 * dx - dy;
        int d1 = 2 * dx;
        int d2 = 2 * (dx - dy);
        while (y != y2) {
            if (d < 0) {
                d += d1;
            } else {
                x += xinc;
                d += d2;
            }
            y += yinc;
            SetPixel(hdc, x, y, c);
        }
    }
}

//parametric polynomial line
void ParametricLine(HDC hdc , int x1 , int y1 , int x2 , int y2 ,int points, COLORREF c){
    int alpha1 = x2-x1;
    int alpha2 = y2-y1;
    double step = 1.0 /(points-1);
    SetPixel(hdc , x1 ,y1 ,c);
    for(double t=0.0 ; t<=1.0 ; t+=step){
        int x = round((alpha1*t)+x1);
        int y = round((alpha2*t)+y1);
        SetPixel(hdc , x, y ,c);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//circle algorithms:
//draw 8 points
void Draw8(HDC hdc , int xc , int yc , int x ,int y , COLORREF c){
    SetPixel(hdc ,xc+x ,yc+y ,c);
    SetPixel(hdc ,xc-x ,yc+y ,c);
    SetPixel(hdc ,xc+x ,yc-y ,c);
    SetPixel(hdc ,xc-x ,yc-y ,c);
    SetPixel(hdc ,xc+y ,yc+x ,c);
    SetPixel(hdc ,xc-y ,yc+x ,c);
    SetPixel(hdc ,xc+y ,yc-x ,c);
    SetPixel(hdc ,xc-y ,yc-x ,c);
}

//normal (cartesian equation) circle
void CircleDrawing(HDC hdc , int xc ,int yc , int r , COLORREF c){
//using second octal
    int x=0;
    int y=r;
    Draw8(hdc , xc, yc,x,y ,c);
    while(x<=y){
        x++;
        y=round(sqrt((double)(r*r)- x*x));
        Draw8(hdc , xc, yc,x,y ,c);
    }
}

//normal polar
void PolarCircle(HDC hdc ,int xc , int yc , int r , COLORREF c){
    //first octal
    double theta = 0.0;
    double dtheta=1.0/r;
    int x=r;
    int y=0;
    Draw8(hdc , xc ,yc ,x,y ,c);
    while(x>y){
        theta+=dtheta;
        x=r*cos(theta);
        y=r*sin(theta);
        Draw8(hdc , xc,yc,x,y,c);
    }
}


//iterative polar

void IterPolar(HDC hdc , int xc , int yc , int r , COLORREF c){
    //first octal
    double dtheta=1.0/r;
    double x=r;
    double y=0;
    Draw8(hdc , xc ,yc ,x ,y ,c);
    double st =sin(dtheta);
    double ct = cos(dtheta);
    while(x >y){
        double x1=x*ct -y*st;
        y=x*st +y*ct;
        x=x1;
        Draw8(hdc , xc ,yc ,round(x) ,round(y) ,c);
    }
}
//bresenham circle
void BresCircle(HDC hdc,int xc,int yc, int r,COLORREF color){
    //second octal
    int x=0,y=r;
    int d=1-r;
    Draw8(hdc,xc,yc,x,y,color);
    while(x<y)
    {
        if(d<0) d += 2*x+2,x++;
        else
        {
            d+=2*(x-y)+5;
            x++;
            y--;
        }
        Draw8(hdc,xc,yc,x,y,color);
    }
}
//modified bresenham
void ModBresCircle(HDC hdc , int xc , int yc , int r , COLORREF c){
    int d=1-r;
    int d1=3;
    int d2 = 5-2*r;
    //second octal
    int x=0 , y=r;
    Draw8(hdc , xc, yc, x,y,c);
    while(x<y){
        if(d<0){
            x++;
            d+=d1;
            d1+=2;
            d2+=2;
            Draw8(hdc , xc, yc, x,y,c);
        }
        else{
            x++;
            y--;
            d+=d2;
            d1+=2;
            d2+=4;
            Draw8(hdc , xc, yc, x,y,c);
        }
    }
}
//Filling Circle with lines after taking filling quarter from user

typedef struct {int xleft,xright;}table[1000];
struct Point{
    int x,y;
    Point(int x=0,int y=0):x(x),y(y){};
};

void DrawHorizontalLine(HDC hdc, int xLeft, int xRight, int y, COLORREF color) {
    for (int x = xLeft; x <= xRight; ++x) {
        SetPixel(hdc, x, y, color);
    }
}


void initTable(table t){
    for (int i = 0; i <1000 ; ++i) {
        t[i].xleft=INT_MAX;
        t[i].xright=-INT_MAX;
    }
}


void edgeToTable(Point p1 , Point p2, table t) {
    if (p2.y == p1.y)return; // ignore horizontal edges
    if (p1.y > p2.y)swap(p1, p2);
    double x = p1.x;
    int y = p1.y;
    double m = (float) (p2.x - x) / (p2.y - p1.y);
    while (y < p2.y) {
        if (x < t[y].xleft)t[y].xleft = (int)ceil(x);
        if (x > t[y].xright)t[y].xright = (int)floor(x);
        y++;
        x += m;
    }
}
void polygonToTable(table t , int n , Point p[]){
    Point v1 = p[n-1];
    for (int i = 0; i < n; ++i) {
        Point v2=p[i];
        edgeToTable(v1,v2,t);
        v1=p[i];
    }
}
void tableToScreen(HDC hdc,table t,COLORREF c){
    for (int i = 0; i < 1000; ++i) {
        if(t[i].xleft<t[i].xright)
            DrawHorizontalLine(hdc,t[i].xleft,t[i].xright ,i ,c);
    }
}

void confex(HDC hdc , Point p[], int n,COLORREF c){
    table t;
    initTable(t);
    polygonToTable(t,n,p);
    tableToScreen(hdc,t,c);
}

void generateQuarterPolygon(HDC hdc,int xc, int yc, int R, int quarter, vector<Point> &polygon) {
    int x = 0, y = R;
    int d = 1 - R;
    while (x <= y) {
        switch (quarter) {
            case 1:
                polygon.push_back({xc + x, yc - y});
                polygon.push_back({xc + y, yc - x});
                break;
            case 2:
                polygon.push_back({xc - x, yc - y});
                polygon.push_back({xc - y, yc - x});
                break;
            case 3:
                polygon.push_back({xc - x, yc + y});
                polygon.push_back({xc - y, yc + x});
                break;
            case 4:
                polygon.push_back({xc + x, yc + y});
                polygon.push_back({xc + y, yc + x});
                break;
        }
        if (d < 0) d += 2 * x + 3;
        else {
            d += 2 * (x - y) + 5;
            y--;
        }
        x++;
    }

    // Connect arc to center to form triangle-like shape
    polygon.push_back({xc, yc});
}



//////////////////////////////////////////////////////////////////////////////////////////
//Filling Algorithms:
//1.Recursive floof fill
void FloodFill1(HDC hdc , int x ,int y , COLORREF currCol , COLORREF FillCol ) {
    COLORREF c = GetPixel(hdc, x, y);
    if (c != currCol || c == FillCol) return;
    else if (c == currCol) {
        SetPixel(hdc, x, y, FillCol);
        FloodFill1(hdc, x + 1, y, currCol, FillCol);  //right
        FloodFill1(hdc, x, y + 1, currCol, FillCol);  //down
        FloodFill1(hdc, x - 1, y, currCol, FillCol);   //left
        FloodFill1(hdc, x, y - 1, currCol, FillCol);    //up
    }
}

//2.non-Recursive flood fill
struct point{
    int x ; int y;
    point(int x=0 , int y=0):x(x),y(y){}
};
void FloodFill2(HDC hdc , int x ,int y , COLORREF currCol , COLORREF FillCol ) {
    stack<point> points;
    points.push(point(x, y));
    while (!points.empty()) {
        point p = points.top();
        points.pop();

        COLORREF c = GetPixel(hdc, p.x, p.y);
        if (c == FillCol || c != currCol) continue;
        SetPixel(hdc, p.x, p.y, FillCol);
        points.emplace(p.x + 1, p.y);
        points.emplace(p.x - 1, p.y);
        points.emplace(p.x, p.y + 1);
        points.emplace(p.x, p.y - 1);

    }
}

///////////////////////////////////////////////////////////////////////////////////////////

int startX, startY;
LineAlgorithm currentAlgo = LINE_DDA;
CircleAlgorithm currentCircleAlgo = CIRCLE_CARTESIAN;
FillAlgorithm currentFillAlgo = FLOOD_FILL_NONRECURSIVE;

COLORREF currentColor = RGB(0, 0, 0);
void DrawShape(HDC hdc, const Shape& shape) {
    if (shape.type == LINE) {
        switch (shape.linealgo) {
            case LINE_DDA:
                DDAline(hdc, shape.x1, shape.y1, shape.x2, shape.y2, shape.color);
                break;
            case LINE_MIDPOINT:
                LineBersenham(hdc, shape.x1, shape.y1, shape.x2, shape.y2, shape.color);
                break;
            case LINE_PARAMETRIC:
                ParametricLine(hdc, shape.x1, shape.y1, shape.x2, shape.y2, 1000, shape.color);
                break;
        }
    } else if (shape.type == CIRCLE) {
        int r = shape.x2; // We stored radius in x2
        switch (shape.circleAlgo) {
            case CIRCLE_CARTESIAN:
                CircleDrawing(hdc, shape.x1, shape.y1, r, shape.color);
                break;
            case CIRCLE_POLAR:
                PolarCircle(hdc, shape.x1, shape.y1, r, shape.color);
                break;
            case CIRCLE_ITERPOLAR:
                IterPolar(hdc, shape.x1, shape.y1, r, shape.color);
                break;
            case CIRCLE_BRES:
                BresCircle(hdc, shape.x1, shape.y1, r, shape.color);
                break;
            case CIRCLE_MODBRES:
                ModBresCircle(hdc, shape.x1, shape.y1, r, shape.color);
                break;
            case FILL_CIRCLE_QUARTER:
                vector<Point> polygon;
                BresCircle(hdc, shape.x1, shape.y1, r, shape.color);  // optional: visualize the boundary
                generateQuarterPolygon(hdc, shape.x1, shape.y1, r, shape.quarter, polygon);
                confex(hdc, polygon.data(), polygon.size(), shape.color);
                break;

        }
    }
    else if (shape.type == FILL) {
        if (shape.fillAlgo == FLOOD_FILL_RECURSIVE) {
            FloodFill1(hdc, shape.x1, shape.y1, shape.currColor, shape.fillColor);
        } else if (shape.fillAlgo == FLOOD_FILL_NONRECURSIVE) {
            FloodFill2(hdc, shape.x1, shape.y1, shape.currColor, shape.fillColor);
        }
    }
}


    bool GetFileNameFromDialog(HWND hwnd, LPWSTR filename, DWORD flags, bool saveDialog = false) {
        OPENFILENAMEW ofn = {};
        ofn.lStructSize = sizeof(ofn);
        ofn.hwndOwner = hwnd;
        ofn.lpstrFilter = L"Binary Files\0*.dat\0All Files\0*.*\0";
        ofn.lpstrFile = filename;
        ofn.nMaxFile = MAX_PATH;
        ofn.Flags = flags;
        ofn.lpstrDefExt = L"dat";

        if (saveDialog)
            return GetSaveFileNameW(&ofn);
        else
            return GetOpenFileNameW(&ofn);
    }


    void SaveShapesToFile(HWND hwnd) {
        wchar_t filename[MAX_PATH] = L"";
        if (GetFileNameFromDialog(hwnd, filename, OFN_OVERWRITEPROMPT, true)) {
            // Convert wide string to narrow string (UTF-8)
            wstring_convert<codecvt_utf8<wchar_t>> converter;
            string narrowFilename = converter.to_bytes(filename);

            ofstream out(narrowFilename, ios::binary);
            if (!out.is_open()) return;

            size_t size = shapes.size();
            out.write(reinterpret_cast<const char *>(&size), sizeof(size));
            for (const auto &shape: shapes) {
                out.write(reinterpret_cast<const char *>(&shape), sizeof(Shape));
            }
            out.close();
        }
    }

    void LoadShapesFromFile(HWND hwnd) {
        wchar_t filename[MAX_PATH] = L"";
        if (GetFileNameFromDialog(hwnd, filename, OFN_FILEMUSTEXIST, false)) {
            // Convert wide string to narrow string (UTF-8)
            wstring_convert<codecvt_utf8<wchar_t>> converter;
            string narrowFilename = converter.to_bytes(filename);

            ifstream in(narrowFilename, ios::binary);
            if (!in.is_open()) return;

            size_t size;
            in.read(reinterpret_cast<char *>(&size), sizeof(size));
            shapes.resize(size);
            for (size_t i = 0; i < size; ++i) {
                in.read(reinterpret_cast<char *>(&shapes[i]), sizeof(Shape));
            }
            in.close();
            InvalidateRect(hwnd, NULL, TRUE);  // Redraw window
        }
    }


    bool fillingMode = false;
    bool drawingCircle = false;
    COLORREF currentFillColor = RGB(255, 0, 0);
    COLORREF currentBorderColor = RGB(0, 0, 0);
    LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {
        HDC hdc;
        switch (msg) {
            case WM_CREATE: {                                      //menu
                HMENU MainMenu = CreateMenu();
                //colors
                HMENU Colors = CreateMenu();
                AppendMenu(Colors, MF_STRING, BLACK, "Black");
                AppendMenu(Colors, MF_STRING, RED, "Red");
                AppendMenu(Colors, MF_STRING, GREEN, "Green");
                AppendMenu(Colors, MF_STRING, BLUE, "Blue");
                AppendMenu(Colors, MF_STRING, YELLOW, "Yellow");
                AppendMenu(Colors, MF_STRING, PINK, "pink");

                //line
                HMENU LineAlgorithms = CreateMenu();
                AppendMenu(LineAlgorithms, MF_STRING, DDA_LINE, "DDA");
                AppendMenu(LineAlgorithms, MF_STRING, MIDPOINT_LINE, "Midpoint");
                AppendMenu(LineAlgorithms, MF_STRING, PARAMETRIC_LINE, "Parametric");
                //circle
                HMENU CircleAlgorithms = CreateMenu();
                AppendMenu(CircleAlgorithms, MF_STRING, CARTESIAN_CIRCLE, "Cartesian Eq Circle");
                AppendMenu(CircleAlgorithms, MF_STRING, POLAR_CIRCLE, "Normal Polar Circle");
                AppendMenu(CircleAlgorithms, MF_STRING, ITERPOLAR_CIRCLE, "Iterative Polar Circle");
                AppendMenu(CircleAlgorithms, MF_STRING, BRES_CIRCLE, "Bresenham Circle");
                AppendMenu(CircleAlgorithms, MF_STRING, MODBRES_CIRCLE, "Modified Bresenham Circle");
                AppendMenu(CircleAlgorithms, MF_STRING, FILL_CIRCLE_QUARTER, "Fill circle");

                // Flood Fill
                HMENU FillAlgorithms = CreateMenu();
                AppendMenu(FillAlgorithms, MF_STRING, FILL_RECURSIVE, "Recursive Flood Fill");
                AppendMenu(FillAlgorithms, MF_STRING, FILL_NONRECURSIVE, "Non-Recursive Flood Fill");


                //main menu
                AppendMenu(MainMenu, MF_POPUP, (UINT_PTR) Colors, "Color");
                AppendMenu(MainMenu, MF_POPUP, (UINT_PTR) LineAlgorithms, "Line Algorithm");
                AppendMenu(MainMenu, MF_POPUP, (UINT_PTR) CircleAlgorithms, "Circle Algorithm");
                AppendMenu(MainMenu, MF_POPUP, (UINT_PTR)FillAlgorithms, "Flood Fill");
                AppendMenu(FillAlgorithms, MF_STRING, FILL_CANCEL, "Cancel Fill Mode");



                AppendMenu(MainMenu, MF_STRING, 401, "Save");
                AppendMenu(MainMenu, MF_STRING, 402, "Load");
                AppendMenu(MainMenu, MF_STRING, CLEAR, "Clear");

                SetMenu(hwnd, MainMenu);
            }
                break;


            case WM_COMMAND: {
                switch (LOWORD(wp)) {
                    case CLEAR:
                        shapes.clear();
                        InvalidateRect(hwnd, NULL, TRUE);
                        break;
                    case SAVE:
                        SaveShapesToFile(hwnd);
                        break;
                    case LOAD:
                        LoadShapesFromFile(hwnd);
                        InvalidateRect(hwnd, NULL, TRUE); // Redraw
                        break;
                    case BLACK:
                        currentColor = RGB(0, 0, 0);
                        break;
                    case RED:
                        currentColor = RGB(255, 0, 0);
                        break;
                    case GREEN:
                        currentColor = RGB(0, 255, 0);
                        break;
                    case BLUE:
                        currentColor = RGB(0, 0, 255);
                        break;
                    case YELLOW:
                        currentColor = RGB(255, 255, 0);
                        break;
                    case PINK:
                        currentColor = RGB(255, 192, 203);
                        break;
                    case DDA_LINE:
                        currentAlgo = LINE_DDA;
                        break;
                    case MIDPOINT_LINE:
                        currentAlgo = LINE_MIDPOINT;
                        break;
                    case PARAMETRIC_LINE:
                        currentAlgo = LINE_PARAMETRIC;
                        break;
                    case CARTESIAN_CIRCLE:
                        drawingCircle = true;
                        currentCircleAlgo = CIRCLE_CARTESIAN;
                        break;
                    case POLAR_CIRCLE:
                        drawingCircle = true;
                        currentCircleAlgo = CIRCLE_POLAR;
                        break;
                    case ITERPOLAR_CIRCLE:
                        drawingCircle = true;
                        currentCircleAlgo = CIRCLE_ITERPOLAR;
                        break;
                    case BRES_CIRCLE:
                        drawingCircle = true;
                        currentCircleAlgo = CIRCLE_BRES;
                        break;
                    case MODBRES_CIRCLE:
                        drawingCircle = true;
                        currentCircleAlgo = CIRCLE_MODBRES;
                        break;
                    case FILL_CIRCLE_QUARTER:
                        drawingCircle = true;
                        currentCircleAlgo = FILL_CIRCLE_QUARTER;
                        break;
                    case FILL_RECURSIVE:
                        currentFillAlgo = FLOOD_FILL_RECURSIVE;
                        fillingMode = true;
                        break;

                    case FILL_NONRECURSIVE:
                        currentFillAlgo = FLOOD_FILL_NONRECURSIVE;
                        fillingMode = true;
                        break;
                    case FILL_CANCEL:
                        fillingMode = false;
                        break;


                }
            }
                break;
            static int circleClickStage = 0;
            case WM_LBUTTONDOWN: {
                int x = LOWORD(lp);
                int y = HIWORD(lp);
                if (drawingCircle && currentCircleAlgo == FILL_CIRCLE_QUARTER) {
                    if (circleClickStage == 0) {
                        startX = x;
                        startY = y;
                        circleClickStage = 1;
                    } else if (circleClickStage == 1) {
                        int dx = x - startX;
                        int dy = y - startY;
                        int radius = (int) round(sqrt(dx * dx + dy * dy));
                        shapes.push_back(
                                {CIRCLE, startX, startY, radius, 0, currentColor, currentAlgo, currentCircleAlgo});
                        circleClickStage = 2;
                    } else if (circleClickStage == 2) {
                        // Determine quarter from mouse click relative to center
                        int quarter = 1;
                        if (x >= startX && y <= startY) quarter = 1;
                        else if (x <= startX && y <= startY) quarter = 2;
                        else if (x <= startX && y >= startY) quarter = 3;
                        else if (x >= startX && y >= startY) quarter = 4;
                        // Update last shape
                        shapes.back().quarter = quarter;
                        circleClickStage = 0;
                        drawingCircle = false;
                        InvalidateRect(hwnd, NULL, TRUE);
                    }
                }
                else if (fillingMode) {
                    hdc = GetDC(hwnd);
                    COLORREF borderColor = GetPixel(hdc, x, y);
                    shapes.push_back({FILL,x, y,0, 0,RGB(0, 0, 0),LINE_DDA,CIRCLE_CARTESIAN,currentFillColor,borderColor,currentFillAlgo,0});                    InvalidateRect(hwnd, NULL, TRUE);
                }
                else {
                    if (!waitingForSecondClick) {
                        startX = x;
                        startY = y;
                        waitingForSecondClick = true;
                    } else {
                        if (drawingCircle) {
                            int dx = x - startX;
                            int dy = y - startY;
                            int radius = (int) round(sqrt(dx * dx + dy * dy));
                            shapes.push_back(
                                    {CIRCLE, startX, startY, radius, 0, currentColor, currentAlgo, currentCircleAlgo});
                            drawingCircle = false;
                        } else {
                            shapes.push_back({LINE, startX, startY, x, y, currentColor, currentAlgo});
                        }
                        waitingForSecondClick = false;
                        InvalidateRect(hwnd, NULL, TRUE);
                    }
                }
            }
            break;

            case WM_PAINT: {
                PAINTSTRUCT ps;
                hdc = BeginPaint(hwnd, &ps);
                for (auto &shape: shapes) {
                    DrawShape(hdc, shape);
                }
                EndPaint(hwnd, &ps);
            }
                break;
            case WM_CLOSE:
                DestroyWindow(hwnd);
                break;
            case WM_DESTROY:
                PostQuitMessage(0);
                break;
            default:
                return DefWindowProc(hwnd, msg, wp, lp);
        }
        return 0;
    }

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow) {
    WNDCLASS wc = {};                   //window class (window type)
    wc.cbClsExtra = 0;
    wc.cbWndExtra = 0;
    wc.hbrBackground = (HBRUSH) GetStockObject(WHITE_BRUSH);   //background color(WHITE_BRUSH or LTGRAY_BRUSH ...)
    wc.hCursor = LoadCursor(NULL, IDC_CROSS);                //cursor shape  IDC_...
    wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
    wc.lpszClassName = "MyClass";                     //L means this is unicode string
    wc.lpszMenuName = NULL;
    wc.lpfnWndProc = WndProc;                 //msg handler(by default)
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.hInstance = hInstance;
    if (!RegisterClass(&wc)) {
        MessageBox(NULL, "Window Registration Failed!", "Error", MB_ICONERROR);
        return 0;
    }
    HWND hwnd = CreateWindow(
            "MyClass", "Your Drawing Helper", WS_OVERLAPPEDWINDOW,
            CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, NULL, NULL, hInstance, NULL);
    if (!hwnd) {
        MessageBox(NULL, "Window Creation Failed!", "Error", MB_ICONERROR);
        return 0;
    }
    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);
    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0)) {        //msg loop
        TranslateMessage(&msg);
        DispatchMessage(&msg);         //call window procedure
    }
    return (int) msg.wParam;
}