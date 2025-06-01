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
#define BEZIER_CURVE 306
#define HERMITE_CURVE 307
#define CARDINAL_SPLINE 308
//#define FILL_CIRCLE_QUARTER 309
#define SAVE 401
#define LOAD 402
#define FILL_RECURSIVE  501
#define FILL_NONRECURSIVE  502
#define FILL_CANCEL   503
#define CONVEX 504
#define NON_CONVEX 505
#define ELLIPSE_POLAR 506
#define ELLIPSE_DIRECT 508
#define ELLIPSE_MID 507



enum ShapeType { LINE ,CIRCLE , FILL, ConvexType, Curves, Ellipsee};
enum LineAlgorithm { LINE_DDA, LINE_MIDPOINT, LINE_PARAMETRIC };
enum EllipseAlgorithm{Ellipse_Polar,Ellipse_Direct,Ellipse_MID};
enum CircleAlgorithm { CIRCLE_CARTESIAN, CIRCLE_POLAR ,CIRCLE_ITERPOLAR,CIRCLE_BRES , CIRCLE_MODBRES ,FILL_CIRCLE_QUARTER};
enum FillAlgorithm { FLOOD_FILL_RECURSIVE, FLOOD_FILL_NONRECURSIVE };
enum ConvexAlgorithm { CONVEX_FILL, NONCONVEX_FILL };
enum CurveAlgorithm{BezierCurve, HermiteCurve,SplineCardinal };


struct Point{
    int x,y;
    Point(int x=0,int y=0):x(x),y(y){};
};

struct Shape {
    ShapeType type;
    int x1, y1, x2, y2;
    COLORREF color;
    LineAlgorithm lineAlgorithm;
    CircleAlgorithm circleAlgo;
    COLORREF fillColor;
    COLORREF currColor;
    int quarter = 1;
    FillAlgorithm fillAlgo;
    ConvexAlgorithm convexAlgo;
    CurveAlgorithm currCurveAlgo;
    int rx,ry;
    vector<POINT> P;
    vector<Point>point;
    vector<Point> vertices;
    int n;
    double c;
    EllipseAlgorithm currEllipseAlgo;

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
void generateQuarterPolygon(HDC hdc,int xc, int yc, int R, int quarter, vector<Point> &polygon) {
    int x = 0, y = R;
    int d = 1 - R;
    while (x <= y) {
        switch (quarter) {
            case 1:
                polygon.emplace_back(xc + x, yc - y);
                polygon.emplace_back(xc + y, yc - x);
                break;
            case 2:
                polygon.emplace_back(xc - x, yc - y);
                polygon.emplace_back(xc - y, yc - x);
                break;
            case 3:
                polygon.emplace_back(xc - x, yc + y);
                polygon.emplace_back(xc - y, yc + x);
                break;
            case 4:
                polygon.emplace_back(xc + x, yc + y);
                polygon.emplace_back(xc + y, yc + x);
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
    polygon.emplace_back(xc, yc);
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

void FloodFill2(HDC hdc , int x ,int y , COLORREF currCol , COLORREF FillCol ) {
    stack<Point> points;
    points.push(Point(x, y));
    while (!points.empty()) {
        Point p = points.top();
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
//Convex And non-Convex Algorithms
// 1.Convex

void polygonToTable(table t , int n ,const Point *p){
    Point v1 = p[n-1];
    for (int i = 0; i < n; ++i) {
        Point v2=p[i];
        edgeToTable(v1,v2,t);
        v1=p[i];
    }
}

void drawConvex(HDC hdc , const Point *p, int n, COLORREF c){
    table t;
    initTable(t);
    polygonToTable(t,n,p);
    tableToScreen(hdc,t,c);
}

//2. Non_Convex
struct Edge{
    double x;
    int yMax;
    double invSlope; //  ->1/m
    Edge(double x, int yMax, double invSlope) : x(x), yMax(yMax), invSlope(invSlope) {}
};
vector<Edge>edgeTable[1000];
vector<Edge>activeEdgeTable;

void buildEdgeTable(Point polygon[],int n){
    for (int i = 0; i < 1000; ++i) edgeTable[i].clear();
    for (int i = 0; i < n; ++i) {
        Point p1=polygon[i];
        Point p2=polygon[(i+1)%n];
        if (p1.y == p2.y) continue; // ignore horizontal

        if (p1.y > p2.y) swap(p1, p2); // p1 is the lower one

        double invSlope = (double)(p2.x - p1.x) / (p2.y - p1.y);
        edgeTable[p1.y].emplace_back(p1.x, p2.y, invSlope);
    }
}

void nonConvex(HDC hdc, Point polygon[], int n, COLORREF c){
    buildEdgeTable(polygon ,n);
    int y=0;
    while (y < 1000 && edgeTable[y].empty()) y++; // find first non-empty scanline

    activeEdgeTable.clear();
    while(!activeEdgeTable.empty() || !edgeTable[y].empty()){
        //step3: Append new edges to active
        for(const Edge& e:edgeTable[y])
            activeEdgeTable.push_back(e);

        //step4: Remove edges where y=ymax
        activeEdgeTable.erase(remove_if(activeEdgeTable.begin(),activeEdgeTable.end(),[y](Edge& e){
            return e.yMax==y;}),activeEdgeTable.end());

        //step5: sort active
        sort(activeEdgeTable.begin(),activeEdgeTable.end(),[](Edge a,Edge b){
            return a.x<b.x;
        });

        // Step 6: Draw spans between pairs of edges
        for (size_t i = 0; i + 1 < activeEdgeTable.size(); i += 2) {
            int xStart = (int)ceil(activeEdgeTable[i].x);
            int xEnd = (int)floor(activeEdgeTable[i + 1].x);
            DrawHorizontalLine(hdc, xStart, xEnd, y, c);
        }

        y++;

        for(auto& edge:activeEdgeTable)
            edge.x+=edge.invSlope;


    }
}

////////////////////////////////////////////////////////////////////////////////////////////
//Hermite & Bezier Curve
int Round(double x) {
    return (int)(x + 0.5);
}
class Vector2
{
public:
    double x, y;
    Vector2(double _x=0, double _y=0): x(_x), y(_y) {}
    Vector2(const Vector2& a){x=a.x;y=a.y;}

};
Vector2 operator-(Vector2& a,Vector2 &b)
{
    double m=a.x-b.x;
    double n=a.y-b.y;
    Vector2 c(m,n);
    return c;
}
Vector2 operator*(double z,Vector2 a)
{
    double m=a.x*z;
    double n=a.y*z;
    Vector2 c(m,n);
    return c;
}
class Vector4
{
public:
    double x1;
    double u1;
    double x2;
    double u2;
    double v[4];
    Vector4(double a=0,double b=0,double c=0,double d=0)
    {
        v[0]=x1=a;
        v[1]=u1=b;
        v[2]=x2=c;
        v[3]=u2=d;
    }
    Vector4(const Vector4& v2)
    {
        memcpy(v,v2.v,4*sizeof(double));
    }
    double&operator[](int index)
    {
        return v[index];
    }
};
class Matrix4
{
public:
    double x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15;
    double m[4][4];
    Matrix4(double a0=0,double a1=0,double a2=0,double a3=0,double a4=0,double a5=0,double a6=0,double a7=0,double a8=0,double a9=0,double a10=0,double a11=0,double a12=0,double a13=0,double a14=0,double a15=0)
    {
        m[0][0]=x0=a0;
        m[0][1]=x1=a1;
        m[0][2]=x2=a2;
        m[0][3]=x3=a3;
        m[1][0]=x4=a4;
        m[1][1]=x5=a5;
        m[1][2]=x6=a6;
        m[1][3]=x7=a7;
        m[2][0]=x8=a8;
        m[2][1]=x9=a9;
        m[2][2]=x10=a10;
        m[2][3]=x11=a11;
        m[3][0]=x12=a12;
        m[3][1]=x13=a13;
        m[3][2]=x14=a14;
        m[3][3]=x15=a15;
    }
    Matrix4(const Matrix4& m2)
    {
        memcpy(m,m2.m,16*sizeof(double));
    }
    double&operator()(int index1,int index2)
    {
        return m[index1][index2];
    }
};
Vector4 operator*(Matrix4 m1,Vector4& v1)
{
    Vector4 res;
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            res[i]+=(v1[j]*m1(i,j));
    return res;
}
double operator*(Vector4 &v1,Vector4& v2)
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3];
}

void DrawHermiteCurve(HDC hdc,double x1,double y1,double u1,double u2,double x2,double y2,double v1,double v2,int n,COLORREF c)
{
    double dt=1.0/n;
    Vector4 Gx(x1,u1,x2,u2);
    Vector4 Gy(y1,v1,y2,v2);
    Matrix4 m(1,0,0,0,0,1,0,0,-3,-2,3,-1,2,1,-2,1);
    Vector4 Vx(m*Gx);
    Vector4 Vy(m*Gy);
    double x,y;
    for(double t=0; t<1; t+=dt)
    {
        Vector4 Vt(1,t,t*t,t*t*t);
        x=Vt*Vx;
        y=Vt*Vy;
        SetPixel(hdc,Round(x),Round(y),c);

    }
}
void DrawHermiteCurve2(HDC hdc,Vector2& P1,Vector2& T1,Vector2& P2,Vector2& T2,int n,COLORREF c)
{
//    double dt=1.0/n;
    Vector4 Gx(P1.x,T1.x,P2.x,T2.x);
    Vector4 Gy(P1.y,T1.y,P2.y,T2.y);
    Matrix4 m(1,0,0,0,0,1,0,0,-3,-2,3,-1,2,1,-2,1);
    Vector4 Vx(m*Gx);
    Vector4 Vy(m*Gy);
    double x,y;
    for(double t=0; t<=1; t+=0.001)
    {
        Vector4 Vt(1,t,t*t,t*t*t);
        x=Vt*Vx;
        y=Vt*Vy;
        SetPixel(hdc,Round(x),Round(y),c);

    }
}


void DrawBezierCubicCurve(HDC hdc,double x1,double y1,double x2,double y2,double x3,double y3,double x4,double y4,int n,COLORREF c)
{
    double u1=3*(x2-x1);
    double v1=3*(y2-y1);
    double u2=3*(x4-x3);
    double v2=3*(y4-y3);
    DrawHermiteCurve(hdc,x1,y1,u1,u2,x4,y4,v1,v2,n,c);
}
///////////////////////////////////////////////////////
//17.Cardinal Spline Curve
void cardinalSplineCurve(HDC hdc, vector<POINT> points, int n , double c, COLORREF color) {
    std::vector<Vector2> p(n);
    for (int i = 0; i < n; ++i) {
        p[i] = Vector2(points[i].x, points[i].y);
    }
    Vector2 *T = new Vector2[n];
    for (int i = 1; i < n - 2; ++i)
        T[i] = c * (p[i + 1] - p[i - 1]);
    T[0] = c * (p[1] - p[0]);
    T[n - 1] = c * (p[n - 1] - p[n - 2]);
    for (int i = 0; i < n - 2; ++i)
        DrawHermiteCurve2(hdc, p[i], T[i], p[i + 1], T[i + 1], n, color);
    delete[] T;
}
////////////////////////////////////////////////////////////////////////////////
//Ellipse Algorithms [Direct, polar and midpoint]
//1.Ellipse Direct
int Max(int x,int y)
{
    if(x>y)
        return x;
    return y;
}
void Draw8EllipsePoints(HDC hdc, int xc, int yc, int x, int y, COLORREF color) {
    SetPixel(hdc, xc + x, yc + y, color);
    SetPixel(hdc, xc - x, yc + y, color);
    SetPixel(hdc, xc + x, yc - y, color);
    SetPixel(hdc, xc - x, yc - y, color);
}
//direct
void DrawEllipse_Direct(HDC hdc, int xc, int yc, int rx, int ry, COLORREF c) {
    int x = 0;
    int steps = (int)(rx / sqrt(1 + (double)(ry * ry) / (rx * rx)));
    while (x <= steps) {
        double y = ry * sqrt(1.0 - (double)(x * x) / (rx * rx));
        Draw8EllipsePoints(hdc, xc, yc, x, round(y), c);
        x++;
    }

    int y = 0;
    steps = (int)(ry / sqrt(1 + (double)(rx * rx) / (ry * ry)));
    while (y <= steps) {
        double xVal = rx * sqrt(1.0 - (double)(y * y) / (ry * ry));
        Draw8EllipsePoints(hdc, xc, yc, round(xVal), y, c);
        y++;
    }
}
//Polar
void DrawEllipse_Polar(HDC hdc,int xc,int yc,int A,int B,COLORREF c)
{
    double dtheta=1.0/Max(A,B);
    for(double theta=0; theta<6.28; theta+=dtheta)
    {
        double x=xc+A*cos(theta);
        double y=yc+B*sin(theta);
        SetPixel(hdc,Round(x),Round(y),c);
    }
}

//Midpoint
void DrawEllipse_Midpoint(HDC hdc, int xc, int yc, int rx, int ry, COLORREF color) {
    int x = 0, y = ry;
    int rx2 = rx * rx;
    int ry2 = ry * ry;
    int two_rx2 = 2 * rx2;
    int two_ry2 = 2 * ry2;

    int px = 0;
    int py = two_rx2 * y;

    // Region 1
    int p = round(ry2 - (rx2 * ry) + (0.25 * rx2));
    while (px < py) {
        Draw8EllipsePoints(hdc, xc, yc, x, y, color);
        x++;
        px += two_ry2;
        if (p < 0) {
            p += ry2 + px;
        } else {
            y--;
            py -= two_rx2;
            p += ry2 + px - py;
        }
    }

    // Region 2
    p = round(ry2 * (x + 0.5) * (x + 0.5) +
              rx2 * (y - 1) * (y - 1) -
              rx2 * ry2);

    while (y >= 0) {
        Draw8EllipsePoints(hdc, xc, yc, x, y, color);
        y--;
        py -= two_rx2;
        if (p > 0) {
            p += rx2 - py;
        } else {
            x++;
            px += two_ry2;
            p += rx2 - py + px;
        }
    }
}



///////////////////////////////////////////////////////////////////////////////////////////
int startX, startY;
LineAlgorithm currentAlgo = LINE_DDA;
CircleAlgorithm currentCircleAlgo = CIRCLE_CARTESIAN;
FillAlgorithm currentFillAlgo = FLOOD_FILL_NONRECURSIVE;
ConvexAlgorithm currConvex= CONVEX_FILL;
CurveAlgorithm currCurveAlgo = BezierCurve;
static int circleClickStage = 0;
EllipseAlgorithm currEllipse=Ellipse_MID;
COLORREF currentColor = RGB(0, 0, 0);
Point ellipseCenter;
bool waitingForndClick = false;

void DrawShape(HDC hdc, const Shape& shape) {
    if (shape.type == LINE) {
        switch (shape.lineAlgorithm) {
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
    } else if(shape.type == Ellipsee){
        switch (shape.currEllipseAlgo) {
            case Ellipse_Direct:
                DrawEllipse_Direct(hdc, shape.x1, shape.y1, shape.rx, shape.ry, shape.color);
                break;
            case Ellipse_Polar:
                DrawEllipse_Polar(hdc, shape.x1, shape.y1, shape.rx, shape.ry, shape.color);
                break;
            case Ellipse_MID:
                DrawEllipse_Midpoint(hdc, shape.x1, shape.y1, shape.rx, shape.ry, shape.color);
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
                drawConvex(hdc, polygon.data(), polygon.size(), shape.color);
                break;

        }
    }
    else if (shape.type == FILL) {
        if (shape.fillAlgo == FLOOD_FILL_RECURSIVE) {
            FloodFill1(hdc, shape.x1, shape.y1, shape.currColor, shape.fillColor);
        } else if (shape.fillAlgo == FLOOD_FILL_NONRECURSIVE) {
            FloodFill2(hdc, shape.x1, shape.y1, shape.currColor, shape.fillColor);
        }
//        else if(shape.fillAlgo ==  Fill_Quarter) {
//            circleClickStage=2;
//            int dx = shape.x1 - startX;
//            int dy = shape.y1 - startY;
//            int r = (int) round(sqrt(dx * dx + dy * dy));
//
////                BresCircle(hdc, shape.x1, shape.y1, r, shape.color);
////            DrawQuarterLine(hdc,shape.x1,shape.y1,shape.x2,shape.y2,shape.quarter,shape.color);
////            BresCircleWithFillingLines(hdc, shape.x1, shape.y1,r , shape.color, shape.quarter);
//        }
    }

    else if(shape.type == ConvexType){
        Point* pts = new Point[shape.n];
        for (int i = 0; i < shape.n; ++i) {
            pts[i] = Point(shape.P[i].x, shape.P[i].y);
        }

        if(shape.convexAlgo == CONVEX_FILL){
            drawConvex(hdc, shape.point.data(), shape.n, shape.fillColor);
        }
        else if(shape.convexAlgo == NONCONVEX_FILL){
            nonConvex(hdc, pts, shape.n, shape.fillColor);
        }
        delete[] pts;

    }

    else if(shape.type == Curves){
        if(shape.currCurveAlgo == SplineCardinal){
            cardinalSplineCurve(hdc, shape.P, shape.n, shape.c, shape.color);
//            DrawBezierCubicCurve(hdc,shape.x1,shape.y1,shape.x2,shape.y2,shape.x3,shape.y3,shape.x4,shape.y4,4,shape.color);

        }
        else if(shape.currCurveAlgo == HermiteCurve){
//            DrawHermiteCurve(hdc,)
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

bool drawingEllipse=false;
bool fillingMode = false;
bool drawingCircle = false;
bool ConvexMode = false;
bool drawingCurve=false;
int idx=0,numPoints=6;


Vector2 *P=new Vector2[numPoints];
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
//                AppendMenu(FillAlgorithms, MF_STRING, FILL_CIRCLE_QUARTER, "Fill circle");

                //Convex
                HMENU ConvexAlgorithms = CreateMenu();
                AppendMenu(ConvexAlgorithms, MF_STRING, CONVEX, "Convex");
                AppendMenu(ConvexAlgorithms, MF_STRING, NON_CONVEX, "Non-Convex");

                //ELLIPSE
                HMENU EllipseAlgorithms = CreateMenu();
                AppendMenu(EllipseAlgorithms, MF_STRING, ELLIPSE_DIRECT , "Direct Ellipse");
                AppendMenu(EllipseAlgorithms, MF_STRING, ELLIPSE_POLAR, "Polar Ellipse");
                AppendMenu(EllipseAlgorithms, MF_STRING, ELLIPSE_MID, "Midpoint Ellipse");

                HMENU CurveAlgorithms = CreateMenu();
                AppendMenu(CurveAlgorithms, MF_STRING, BEZIER_CURVE, "Bezier Curve");
                AppendMenu(CurveAlgorithms, MF_STRING, HERMITE_CURVE, "Hermite Curve");
                AppendMenu(CurveAlgorithms, MF_STRING, CARDINAL_SPLINE, "Cardinal Spline");

                //main menu
                AppendMenu(MainMenu, MF_POPUP, (UINT_PTR) Colors, "Color");
                AppendMenu(MainMenu, MF_POPUP, (UINT_PTR) LineAlgorithms, "Line Algorithm");
                AppendMenu(MainMenu, MF_POPUP, (UINT_PTR) CircleAlgorithms, "Circle Algorithm");
                AppendMenu(MainMenu, MF_POPUP, (UINT_PTR)FillAlgorithms, "Flood Fill");
                AppendMenu(MainMenu, MF_POPUP, (UINT_PTR)ConvexAlgorithms, "Convex Algorithm");
                AppendMenu(MainMenu, MF_POPUP, (UINT_PTR)EllipseAlgorithms, "Ellipse Algorithm");
                AppendMenu(MainMenu, MF_POPUP, (UINT_PTR)CurveAlgorithms, "Curve Algorithm");

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
                    case BEZIER_CURVE:
                        drawingCurve=true;
                        currCurveAlgo = BezierCurve;
                        break;
                    case HERMITE_CURVE:
                        drawingCurve=true;
                        currCurveAlgo = HermiteCurve;
                        break;
                    case CARDINAL_SPLINE:
                        drawingCurve=true;
                        currCurveAlgo = SplineCardinal;
                        break;
                    case CONVEX:
                        ConvexMode=true;
                        currConvex=CONVEX_FILL;
                        break;
                    case NON_CONVEX:
                        ConvexMode=true;
                        currConvex=NONCONVEX_FILL;
                        break;
                    case ELLIPSE_MID:
                        drawingEllipse= true;
                        currEllipse=Ellipse_MID;
                        break;
                    case ELLIPSE_POLAR:
                        drawingEllipse= true;
                        currEllipse=Ellipse_Polar;
                        break;
                    case ELLIPSE_DIRECT:
                        drawingEllipse= true;
                        currEllipse=Ellipse_Direct;
                        break;
                    case FILL_CANCEL:
                        fillingMode = false;
                        ConvexMode=false;
                        break;
                }
            }
                break;
                static Point polygon[5];
                static int count = 0;
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
                    shapes.push_back(
                            {FILL, x, y, 0, 0, RGB(0, 0, 0), LINE_DDA, CIRCLE_CARTESIAN, currentFillColor, borderColor,
                             1, currentFillAlgo});
                }
                else if (ConvexMode){
                    polygon[count].x = LOWORD(lp);
                    polygon[count].y = HIWORD(lp);
                    if (count ==4)
                    {
                        hdc = GetDC(hwnd);
                        drawConvex(hdc, polygon, 5, currentColor);
                        ReleaseDC(hwnd, hdc);
                        count = 0;
                    }
                    else {
                        count++;
                    }
                }
//                else if(currEllipse){
//                    if (!waitingForndClick) {
//                        startX = LOWORD(lp);
//                        startY = HIWORD(lp);
//                        waitingForndClick = true;
//                    } else {
//                        int rx = abs(startX - x);
//                        int ry = abs(startY-y);
//
//                        Shape ellipseShape;
//                        ellipseShape.type = Ellipsee;
//                        ellipseShape.x1 = startX;
//                        ellipseShape.y1 = startY;
//                        ellipseShape.rx = rx;
//                        ellipseShape.ry = ry;
//                        ellipseShape.color = currentColor;
//                        ellipseShape.currEllipseAlgo = currEllipse;
//                        shapes.push_back(ellipseShape);
//                        drawingEllipse = false;
//                    }
//                    waitingForSecondClick = false;
//                    InvalidateRect(hwnd, NULL, TRUE);
//                }
                else if(drawingCurve){
                    if(idx<numPoints)
                    {
                        P[idx].x=LOWORD(lp);
                        P[idx].y=HIWORD(lp);
                        idx++;
                    }
                    else {
                        idx=0;
                        hdc = GetDC(hwnd);
                        std::vector<POINT> pointVec;
                        for (int i = 0; i < numPoints; ++i) {
                            POINT pt;
                            pt.x = static_cast<LONG>(P[i].x);
                            pt.y = static_cast<LONG>(P[i].y);
                            pointVec.push_back(pt);
                        }
                        Shape s;
                        s.type = Curves;
                        s.currCurveAlgo = SplineCardinal;
                        s.P = pointVec;
                        s.n = numPoints;

                        s.c = 1.0; // or any tension value you want
                        s.color = currentColor;
                        shapes.push_back(s);
                        ReleaseDC(hwnd, hdc);
                        InvalidateRect(hwnd, NULL, TRUE);
                        drawingCurve= false;
                    }

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
                        }
                        else if(drawingEllipse){
                            int rx = abs(x - startX);
                            int ry = abs(y - startY);

                            Shape ellipseShape;
                            ellipseShape.type = Ellipsee;
                            ellipseShape.x1 = startX;
                            ellipseShape.y1 = startY;
                            ellipseShape.rx = rx;
                            ellipseShape.ry = ry;
                            ellipseShape.color = currentColor;  // or ellipseShape.currColor
                            ellipseShape.currEllipseAlgo = currEllipse;
                            shapes.push_back(ellipseShape);
                            drawingEllipse = false;

                        }
                        else {
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