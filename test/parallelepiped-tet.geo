SetFactory("OpenCASCADE");

h = 10;     // base length
l = 20;     // height

Point(1) = {0, h/2, h/2};
Point(2) = {l, h/2, h/2};
Point(3) = {l, h/2, -h/2};
Point(4) = {0, h/2, -h/2};
Point(5) = {0, -h/2, h/2};
Point(6) = {l, -h/2, h/2};
Point(7) = {l, -h/2, -h/2};
Point(8) = {0, -h/2, -h/2};

Point(12) = {0, 0, 0};
Point(13) = {l, 0, 0};
Point(14) = {0, h/2, 0};
Point(15) = {0, 0, h/2};

Line(1) = {5, 15};
Line(2) = {15, 1};
Line(3) = {1, 14};
Line(4) = {14, 4};
Line(5) = {4, 8};
Line(6) = {8, 5};
Line(9) = {5, 6};
Line(10) = {6, 2};
Line(11) = {2, 3};
Line(12) = {3, 7};
Line(13) = {7, 6};
Line(14) = {8, 7};
Line(15) = {4, 3};
Line(16) = {1, 2};

Line Loop(1) = {2, 3, 4, 5, 6, 1};
Plane Surface(1) = {1};

Line Loop(2) = {9, 10, -16, -2, -1};
Plane Surface(2) = {2};

Line Loop(3) = {-13, -10, -11, -12};
Plane Surface(3) = {3};

Line Loop(4) = {15, 12, -14, -5};
Plane Surface(4) = {4};

Line Loop(5) = {11, -15, -4, -3, 16};
Plane Surface(5) = {5};

Line Loop(6) = {14, 13, 9, 6};
Plane Surface(6) = {6};

Point {12} In Surface {1};
Point {13} In Surface {3};

Surface Loop(1) = {6, 4, 5, 3, 2, 1};
Volume(1) = {1};

// pretty names
Physical Point("O") = {12};
Physical Point("A") = {13};
Physical Point("B") = {14};
Physical Point("C") = {15};
Physical Point("D") = {2};

Physical Surface("x=0") = {1};
//Physical Surface("z=h/2") = {2};
Physical Surface("x=l") = {3};
//Physical Surface("z=-h/2") = {4};
//Physical Surface("y=h/2") = {5};
//Physical Surface("y=-h/2") = {6};

Physical Volume("bulk") = {1};

// meshing options
Mesh.Algorithm = 8;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.HighOrderOptimize = 2;

Mesh.CharacteristicLengthMin = Mesh.ElementOrder * l/16;
Mesh.CharacteristicLengthMax = Mesh.ElementOrder * l/16;
