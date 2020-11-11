SetFactory("OpenCASCADE");

r = 1;
a = 1;
b = 10;
c = 2;
d = 4;

H = 0.5;

Point(1) = {0, 0, 0};
Point(2) = {r, 0, 0};
Point(3) = {0, r, 0};

Point(4) = {r/Sqrt(2), r/Sqrt(2), 0};

Point(5) = {0, a+r, 0};
Point(6) = {a+r, 0, 0};
Point(7) = {a+r, a+r, 0};

Circle(1) = {3, 1, 4};
Circle(2) = {4, 1, 2};
Line(3) = {2, 6};
Line(4) = {6, 7};
Line(5) = {7, 5};
Line(6) = {5, 3};

Line(7) = {4, 7};

Point(8) = {a+r+b, 0, 0};
Point(9) = {a+r+b, c, 0};
Point(10) = {a+r+b, a+r+d, 0};
Point(11) = {c, a+r+d, 0};
Point(12) = {0, a+r+d, 0};

Line(8) = {6, 8};
Line(9) = {8, 9};
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 5};
Line(14) = {7, 9};
Line(15) = {7, 11};

Curve Loop(1) = {2, 3, 4, -7};
Plane Surface(1) = {1};
Curve Loop(2) = {7, 5, 6, 1};
Plane Surface(2) = {2};

Curve Loop(3) = {8, 9, -14, -4};
Plane Surface(3) = {3};

Curve Loop(4) = {15, 12, 13, -5};
Plane Surface(4) = {4};

Curve Loop(5) = {14, 10, 11, -15};
Plane Surface(5) = {5};
