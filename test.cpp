#define CATCH_CONFIG_MAIN

#include <iostream>

#include "catch.hpp.h"
#include "shape.h"


typedef std::vector<Point> points;

void alert_point(Point p, Point expected)
{
    CHECK(p.x == Approx(expected.x));
    REQUIRE(p.y == Approx(expected.y));
}

void alert_intersection(BaseFigure& b_right, BaseFigure& b_left, std::vector<Point> expected)
{    
    std::vector<Point> result = b_right.intersection(b_left);
    REQUIRE(result.size() == expected.size());
    for (int i = 0; i < result.size(); i++) {
        alert_point(result[i], expected[i]);
    };
}

void alert_intersection(Point p1_begin, Point p1_end, Point p2_begin, Point p2_end, points expected)
{
    Segment s1(p1_begin, p1_end), s2(p2_begin, p2_end);
    alert_intersection(s1, s2, expected);
};

void alert_intersection(double r1, Point c1, double r2, Point c2, points expected)
{
    Circle c(r1, c1), c_(r2, c2);
    alert_intersection(c, c_, expected);
};


TEST_CASE("d") {
    std::vector<std::shared_ptr<BaseFigure>> p;

    p.emplace_back(new Circle(34, Point(12, 12)));
    p.emplace_back(new Segment(Point(12, 12), Point(23, 23)));
    p.emplace_back(new Circle(1134, Point(12, 12)));
    p.emplace_back(new Segment(Point(122, 122), Point(-122, 22)));

    for (int i = 0; i < p.size(); i++) {
        for (int j = i + 1; j < p.size(); ++j) {
            (*p[i]).intersection(*p[j]);
        }
    }


}

TEST_CASE("create object and initializing variables", "[create], [initializing]") 
{    
    SECTION("Line")
    {
        Point p1(1.3455, 6755.000);
        Point p2(-344.5, 0.99);
        Segment s1(p1, p2);

        alert_point(s1.get_point_begin(), p2);
        alert_point(s1.get_point_end(), p1);

        REQUIRE(s1.get_A() == Approx(0.998691548));
        REQUIRE(s1.get_B() == Approx(-0.05113895));
        REQUIRE(s1.get_C() == Approx(344.099865808));
    };
    
    SECTION("Circle")
    {
        Point p3(3473.34334, -994.33334);
        double r = 358745.49;
        Circle c1(r, p3);

        alert_point(c1.get_center(), p3);
        REQUIRE(c1.get_radius() == Approx(r));
    };
    
    SECTION("Broken Line")
    {
        std::vector<Point> nodes;
        for (int i = 0; i < 20; i++) {
            nodes.push_back(Point(std::rand(), std::rand()));
        };
        PolyLine br(nodes);

        for (int i = 0; i < nodes.size(); i++) {
            alert_point(br.get_node(i), nodes[i]);
        };
        REQUIRE(br.get_num_nodes() == nodes.size());
    }
};

void alert_intersection(Point p_begin, Point p_end, double r, Point center, points expected){
    Segment segment(p_begin, p_end);
    Circle circle(r, center);
    alert_intersection(segment, circle, expected);
}

TEST_CASE("intersect segment-segment", "[intersection], [segment]")
{
    points expected;        

    SECTION("intersect point not contain in segments")
    {
        alert_intersection(
            Point(21.0, 226),
            Point(43, 112),
            Point(58, 115),
            Point(98, 24),
            expected
        );
    };

    SECTION("lines parallel")
    {
        Point p1(78.9883, 1.8775);
        Point p2(-54.31, -343.99);
        alert_intersection(
            p1, p2, 
            (-758.26) * p1,
            (-785.26) * p2,
            expected
        );
        
        p1.set(833.9883, 109.8775);
        p2.set(444.31, 323.99);
        alert_intersection(
            p1, p2,
            Point(-100, -100) + p1,
            Point(-100, -100) + p2,
            expected
        );
    };
    SECTION("line equivalents")
    {
        alert_intersection(
            Point(89.33333, 89.3333), 
            Point(903.444444, 903.444444),
            Point(-5666.666, -5666.666),
            Point(-7.333393, -7.333393),
            expected
        );
    };
    
    SECTION("simple intersect")
    {      
        expected = {Point(7.74343, 8.04034)};
        alert_intersection(
            Point(3.874203 + 3 * (8.4534534 - 3.874203),  38.38574 + 3 * (9.84752 - 38.38574)),
            Point(3.874203, 38.38574),
            Point(45.45734, 45.756345),
            Point(-45.344, -47.5234),
            expected
        );

        expected = {Point(4.16395, 0.46455)};
        alert_intersection(
            Point(0.3 * 12.284, 0.3 * 1.5485),
            Point(0.3 * 15.5485, 0.3 *  1.5485),
            Point(4.16395, 0.46455),
            Point(0.3 * 13.879848, 0.3 *  20.5154),
            expected
        );
    };

    SECTION("segments on one line")
    {
        expected = {Point(-0.55, 6), Point(23.024556, 6)};
        alert_intersection(
            Point(-12.5555, 6),
            Point(23.024556, 6),
            Point(43.5464, 6),
            Point(-0.55, 6),
            expected
        );

        expected = {Point(-5, -5.354643), Point(-5, 65.34345)};
        alert_intersection(
            Point(-5, 89.345342),
            Point(-5, -90.435234),
            Point(-5, 65.34345),
            Point(-5, -5.354643),
            expected
        );
    }
    

};



TEST_CASE("intersect circle circle")
{
    points expected;
    
    SECTION("simple intersect")
    {
        expected = {Point(88.7591, 102.447), Point(94.0211, 91.478)};
        alert_intersection(
            10,
            Point(98.54646, 100.39545),
            7,
            Point(88.26684, 95.464),
            expected
        );

        expected = {Point(-421.554, -27.242), Point(-418.292, 34.1061)};
        alert_intersection(
            44.554,
            Point(-452.15, 5.1455),
            32.54,
            Point(-430.6456, 4.00212),
            expected
        );
    };

    SECTION("no intersect")
    {
        expected = {};
        alert_intersection(
            4.554,
            Point(-452.15, 5.1455),
            25.54,
            Point(-9.6456, 9.00212),
            expected
        );
    };

    SECTION("nested circle")
    {   
        expected = {};
        alert_intersection(
            44.554,
            Point(-452.15, 5.1455),
            1.54,
            Point(-430.6456, 4.00212),
            expected
        );
    };

    SECTION("external contact")
    {
        expected = {Point(-27.0/50.0, 21111.0/2000)};
        alert_intersection(
            5.54,
            Point(5, 10.5555),
            10.56,
            Point(-11.1, 10.5555),
            expected
        );
    };
    SECTION("internal contact")
    {
        expected = {Point(-50.54, 10.55555)};
        alert_intersection(
            55.54,
            Point(5, 10.5555),
            45.55555,
            Point(-4.98445, 10.5555),
            expected
        );  
    };
}



TEST_CASE("intersect circle segment")
{
    points expected;
    SECTION("simple intersect")
    {
        expected = {Point(-6.85125, 12.3884), Point(19.8913, -7.66844)};
        alert_intersection(
            Point(23.0, -10),
            Point(-25.0, 26),
            23.0,
            Point(16.0, 15.0),
            expected
        );
    };

    SECTION("intersect with one point and segment inside circle")
    {
        expected =  {Point(-30.8029, 38.7436)};
        alert_intersection(
            Point(23.0054, -4.51211), 
            Point(-89.5879, 86), 
            45.55555, 
            Point(4.98445, 10.5555), 
            expected
        );
    };  

    SECTION("no intersect")
    {
        expected = {};
        alert_intersection(
            Point(15, -10),
            Point(-25, 26),
            12, 
            Point(16, 15),
            expected
        );
    };

    SECTION("contact")
    {
        expected = {Point(4.0, 2.0)};
        alert_intersection(
            Point(4.0, 5.0),
            Point(4.0, 1.0),
            2.0,
            Point(2.0, 2.0),
            expected
        );
    };
};

TEST_CASE("intersect segment with broken line") 
{
    points nodes = {Point(0.0, 0.0), Point(4.0, 4.0), Point(2.0, 0.0), Point(2.0, 4.0)};
    PolyLine broken_line(nodes);

    Point a(0.0, 3.0), b(4.0, 3.0);
    Segment seg(a, b);

    points expected = {Point(3.0, 3.0), Point(3.5, 3.0), Point(2.0, 3.0)};
    alert_intersection(broken_line, seg, expected);
};

TEST_CASE("intersect broken line with broken line") 
{
    points nodes1 = {Point(6.0, 4.0), Point(0.0, 1.0), Point(5.0, 1.0)};
    points nodes2 = {Point(1.0, 3.0), Point(4.0, 0.0), Point(4.0, 4.0)};
    PolyLine bline1(nodes1), bline2(nodes2);

    points expected = {Point(2.0, 2.0), Point(4.0, 3.0), Point(3.0, 1.0), Point(4.0, 1.0)};
    alert_intersection(bline1, bline2, expected);
}



TEST_CASE("intersect broken line with circle") 
{
    points nodes = {Point(2.0, 5.0), Point(2.0, -1.0), Point(4.0, -1.0), Point(4.0, 5.0)};
    PolyLine broken_line(nodes);
    Point a(2.0, 2.0);
    Circle circle(2.0, a);

    points expected = {Point(2.0, 0.0), Point(2.0, 4.0), Point(4.0, 2.0)};
    alert_intersection(broken_line, circle, expected);
}; 
