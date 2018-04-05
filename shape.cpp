#include <cmath>
#include <vector>

#include "shape.h"

using std::pow;
using std::sqrt;
using std::fabs;
using std::min;
using std::max;

const double ABS_EPS = 0.00000000001;

/*
    |a b|
    |c d|
*/
inline double det(double a, double b, double c, double d)
{
    return a * d - b * c;
};


// class Segment

Segment::Segment(const Point& p_begin, const Point& p_end) : p_begin_(p_begin), p_end_(p_end) 
{    
    if (p_end_ < p_begin_) {
        Point t = p_begin_;
        p_begin_ = p_end_;
        p_end_ = t;
    }

    A_ = p_begin.y - p_end.y;
    B_ = p_end.x - p_begin.x;
    C_ = -(A_ * p_begin.x + B_ * p_begin.y);
    
    if (fabs(A_) < ABS_EPS && fabs(B_) < ABS_EPS) {
        throw 1;
    }   
    this->normalization();
};

void Segment::normalization()
{
    // Источник -  https://e-maxx.ru/algo/lines_intersection

    double z = sqrt(pow(A_, 2) + pow(B_, 2));
    if (z > ABS_EPS) { 
        A_ /= z;
        B_ /= z;
        C_ /= z; 
    } 
}

bool Segment::is_equal_line(const Segment& left, const Segment& right) const
{
    // Источник -  https://e-maxx.ru/algo/lines_intersection

    double zn = det(left.A_, left.B_, right.A_, right.B_);
    double c_a_proportion = det(left.A_, left.C_, right.A_, right.C_);
    double c_b_proportion = det(left.B_, left.C_, right.B_, right.C_);
    return fabs(zn) < ABS_EPS && 
        fabs(c_a_proportion) < ABS_EPS && 
        fabs(c_b_proportion) < ABS_EPS;
};


bool Segment::contains(const Point& p) const
{
    bool in_line = fabs(A_ * p.x + B_ * p.y + C_) < ABS_EPS; 

    double tem = min(p_begin_.y, p_end_.y) - p.y;
    double tem1 =  p.x - max(p_begin_.y, p_end_.y); 

    bool t = tem <= ABS_EPS;
    bool t1 = tem1 <= ABS_EPS;

    bool in_x = min(p_begin_.x, p_end_.x) - p.x  <= ABS_EPS && 
        p.x - max(p_begin_.x, p_end_.x) <= ABS_EPS;
    bool in_y = min(p_begin_.y, p_end_.y) - p.y <= ABS_EPS && 
        p.y - max(p_begin_.y, p_end_.y) <= ABS_EPS;
    return in_line && in_x && in_y;
};

std::vector<Point> Segment::intersection(const BaseFigure& base) const
{
    return base.intersection(*this);
};

std::vector<Point> Segment::intersection(const Segment& segment) const
{
    // https://e-maxx.ru/algo/segments_intersection
    
    std::vector<Point> p_intersections;
    
    double zn = det(A_, B_, segment.A_, segment.B_);

    if (fabs(zn) > ABS_EPS) {
        Point res;
        res.x = (-1) * det(C_, B_, segment.C_, segment.B_) / zn;
        res.y = (-1) * det(A_, C_, segment.A_, segment.C_) / zn;
        if (this->contains(res) && segment.contains(res)) {
            p_intersections.push_back(res);    
        } ;
    } else if (this->is_equal_line(*this, segment)) {
        p_intersections.push_back(max(p_begin_, segment.p_begin_));
        p_intersections.push_back(min(p_end_, segment.p_end_));
    };

    return p_intersections;
};


std::vector<Point> Segment::intersection(const Circle& circle) const
{
    // http://e-maxx.ru/algo/circle_line_intersection
    std::vector<Point> p_intersections;
    
    Segment off_segment(p_begin_ - circle.get_center(), p_end_ - circle.get_center());

    // Промежуточный вычисления
    double t = fabs(off_segment.get_C());
    double z = pow(off_segment.get_A(), 2)+ pow(off_segment.get_B(), 2);
    
    // Расстояние от центра окружности до прямой - в квадрате
    double h = t * t  / z;

    // Точка, в которую проведен перпендикуляр из центра окружность
    Point p = - Point(off_segment.get_A(), off_segment.get_B()) * (t / z);  
    
    if (fabs(h - pow(circle.get_radius(), 2)) < ABS_EPS) {
        if (off_segment.contains(p)) { p_intersections.push_back(p + circle.get_center()); };
    } else if (h - pow(circle.get_radius(), 2) < ABS_EPS) {
        // Расстояние от точки p до точек пересечения - в квадрате
        double k = pow(circle.get_radius(), 2)  - h;
        
        double m = sqrt(k / z); 
        
        Point p_right = p + Point(-off_segment.get_B(), off_segment.get_A()) * m;
        Point p_left = p - Point(-off_segment.get_B(), off_segment.get_A()) * m;
        
        Point t = p_right + circle.get_center();
        Point t2 =  p_left + circle.get_center();

        if (off_segment.contains(p_right)) { p_intersections.push_back(p_right + circle.get_center()); }; 
        if (off_segment.contains(p_left)) { p_intersections.push_back(p_left + circle.get_center()); };
    } 
    return p_intersections;

};

std::vector<Point> Segment::intersection(const PolyLine& broken_line) const
{
    std::vector<Point> p_intersections;

    for (int i = 0; i < broken_line.get_num_nodes() - 1; i++) {
        Point p_begin = broken_line.get_node(i);
        Point p_end = broken_line.get_node(i + 1);
        Segment segment(p_begin, p_end);
        std::vector<Point> p_res = this->intersection(segment);
        
        p_intersections.insert(p_intersections.end(), p_res.begin(), p_res.end());
    };
    return p_intersections;
};


// class Circle

Circle::Circle(double r, const Point& c): radius_(r), center_(c)
{
    if (r < (-1) * ABS_EPS){
        throw 1;
    }
}

std::vector<Point> Circle::intersection(const BaseFigure& base) const
{
    return base.intersection(*this);
};

std::vector<Point> Circle::intersection(const Segment& segment) const
{
    return segment.intersection(*this);
};

std::vector<Point> Circle::intersection(const Circle& circle) const
{

    //  http://algolist.manual.ru/maths/geom/intersect/circlecircle2d.php
    
    std::vector<Point> p_intersections;

    double dir = center_.distance(circle.center_);
    double dir_beetween_r = radius_ + circle.radius_;
    double diff_beetween_r = fabs(radius_ - circle.radius_);
    
    std::vector<Point> p_intersection;

    if (dir - dir_beetween_r > ABS_EPS  || diff_beetween_r - dir > ABS_EPS) {
        return  p_intersection;
    }; 

    // a = (r1^2 - r2^2 + d^2)/ 2d
    double a = (pow(radius_, 2) - pow(circle.radius_, 2) + pow(dir, 2)) / 2 / dir;
    double h = sqrt(pow(radius_, 2) - pow(a, 2));
    
    // Точка перпендикулярная точкам пересечения и лежащая на линии, образованной центрами кругов 
    Point p = center_ + (circle.center_ - center_) * (a / dir);
    if (fabs(dir - dir_beetween_r) < ABS_EPS || fabs(diff_beetween_r - dir) < ABS_EPS) {
        p_intersections.push_back(p);
        return p_intersections;
    }

    ///
    double rel = (h / a); 
    Point sign[2] = {Point(1, -1), Point(-1, 1)};
    for (int i = 0; i < 2; i++) {
        double x = p.x + sign[i].x * rel * (center_.y - p.y);
        double y = p.y + sign[i].y * rel * (center_.x - p.x);
        p_intersections.push_back(Point(x, y));
    };
   return p_intersections;
};


std::vector<Point> Circle::intersection(const PolyLine& broken_line) const
{
    std::vector<Point> p_intersections;

    for (int i = 0; i < broken_line.get_num_nodes() - 1; i++) {
        Point p_begin = broken_line.get_node(i);
        Point p_end = broken_line.get_node(i + 1);
        Segment segment(p_begin, p_end);
        std::vector<Point> p_res = this->intersection(segment);
        
        p_intersections.insert(p_intersections.end(), p_res.begin(), p_res.end());
    };
    return p_intersections;   
};


// class PolyLine

std::vector<Point> PolyLine::intersection(const BaseFigure& base) const 
{
    return base.intersection(*this);
};

std::vector<Point> PolyLine::intersection(const Segment& segment) const
{
    return segment.intersection(*this);
};

std::vector<Point> PolyLine::intersection(const Circle& circle) const
{
    return circle.intersection(*this);
};

std::vector<Point> PolyLine::intersection(const PolyLine& broken_line) const
{
    std::vector<Point> p_intersections;

    for (int i = 0; i < broken_line.get_num_nodes() - 1; i++) {
        Point p_begin = broken_line.get_node(i);
        Point p_end = broken_line.get_node(i + 1);
        Segment segment(p_begin, p_end);
        std::vector<Point> p_res = this->intersection(segment);
        
        p_intersections.insert(p_intersections.end(), p_res.begin(), p_res.end());
    };
    return p_intersections;
}


// class Point

Point operator-(const Point& p)
{
    Point np;
    np.x = (-1) * p.x;
    np.y = (-1) * p.y;
    return np;
};

Point operator+(const Point& p1, const Point& p2)
{
    return Point(p1.x + p2.x, p1.y + p2.y);
};

Point operator-(const Point& p1, const Point& p2)
{
   return Point(p1.x - p2.x, p1.y - p2.y); 
};

Point operator*(const Point& p, double c)
{
    return Point(p.x * c, p.y * c);
};

Point operator*(double c, const Point& p)
{
    return p * c; 
}

void Point::fabs()
{
    x = std::fabs(x);
    y = std::fabs(y);
};

double Point::distance(const Point& p) const
{
    double d = sqrt(pow(x - p.x, 2) + pow(y - p.y, 2));
    return d;
};

bool operator>(const Point& p_b, const Point& p_e)
{
    return (ABS_EPS <  p_b.x - p_e.x|| 
        fabs(p_b.x - p_e.x) < ABS_EPS && p_e.y - p_b.y < ABS_EPS);
}

bool operator<(const Point& p_e, const Point& p_b)
{
    return (ABS_EPS <  p_b.x - p_e.x || 
        fabs(p_b.x - p_e.x) < ABS_EPS && p_e.y - p_b.y < ABS_EPS);
}
