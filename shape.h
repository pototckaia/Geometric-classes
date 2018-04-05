#pragma once

#include <vector>


class Point
{
public:
    double x, y;
    
    Point(): x(0), y(0) {};
    Point(double x, double y) : x(x), y(y) {};
    Point(const Point& a) : x(a.x), y(a.y) {};
        
    void fabs();
    double distance(const Point&) const;
    void set(double x, double y) { this->x = x; this->y = y; };
    
    friend Point operator-(const Point&);
    friend Point operator+(const Point&, const Point&); 
    friend Point operator-(const Point&, const Point&); 
    friend Point operator*(const Point&, double);
    friend Point operator*(double, const Point&);
    
    friend bool operator>(const Point&, const Point&);
    friend bool operator<(const Point&, const Point&);
    
    Point& operator=(const Point& right)
    {
        //проверка на самоприсваивание
        if (this == &right) {
            return *this;
        }
        x = right.x;
        y = right.y;
        return *this;
    };
};


class Segment;
class Circle;
class PolyLine;

class BaseFigure
{
public:
    BaseFigure() {};
    virtual ~BaseFigure() {};

    virtual std::vector<Point> intersection (const BaseFigure&) const =0;
    virtual std::vector<Point> intersection (const Segment&) const =0;
    virtual std::vector<Point> intersection (const Circle&) const =0;
    virtual std::vector<Point> intersection (const PolyLine&) const =0;
};


class Segment: public BaseFigure
{
public:
    Segment(const Point&, const Point&);
    virtual ~Segment() {};

    std::vector<Point> intersection(const BaseFigure&) const override;
    std::vector<Point> intersection(const Segment&) const override;
    std::vector<Point> intersection(const Circle&) const override;
    std::vector<Point> intersection(const PolyLine&) const override;
        
    bool is_equal_line(const Segment&, const Segment&) const;   
    bool contains(const Point&) const;

    double get_A() { return A_; };
    double get_B() { return B_; };
    double get_C() { return C_; };
    
    const Point& get_point_begin() const { return p_begin_; };
    const Point& get_point_end() const { return p_end_; };

private:
    double A_, B_, C_; // Ax + By + C = 0
    Point p_begin_, p_end_;
    void normalization();
};


class Circle: public BaseFigure
{
public:
    Circle(double, const Point&);
    virtual ~Circle() {};

    std::vector<Point> intersection(const BaseFigure&) const override;
    std::vector<Point> intersection(const Segment&) const override;
    std::vector<Point> intersection(const Circle&) const override;
    std::vector<Point> intersection(const PolyLine&) const override;
        
    double get_radius() const { return radius_; };
    const Point& get_center() const { return center_; };

private:
    double radius_;
    Point center_;
};



class PolyLine: public BaseFigure
{
public:
    PolyLine(const std::vector<Point> nodes) : nodes_(nodes), num_nodes_(nodes.size()) {}
    virtual ~PolyLine() {};

    std::vector<Point> intersection(const BaseFigure&) const override;
    std::vector<Point> intersection(const Segment&) const override;
    std::vector<Point> intersection(const Circle&) const override;
    std::vector<Point> intersection(const PolyLine&) const override;

    double get_num_nodes() const { return num_nodes_; };
    const Point& get_node(int i) const { return nodes_[i]; };
private:
    std::vector<Point> nodes_;
    int num_nodes_;
};
