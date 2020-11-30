SetFactory("OpenCASCADE");

l = 60;
D = 28;
r = 1.0;

a = 6;
b = 1;

Point(1) = {0, 0, 0};
Point(2) = {D/2-b-r, 0, 0};
Point(4) = {(D/2-b-r)*Sqrt(2)/2, (D/2-b-r)*Sqrt(2)/2, 0};
Point(5) = {0, (D/2-b-r), 0};
Point(6) = {1.2*a, 0, 0};
Point(7) = {a, a, 0};
Point(8) = {0, 1.2*a, 0};

Line(2) = {6, 2};
Line(3) = {6, 7};
Line(4) = {4, 7};

Circle(6) = {2, 1, 4};


Point(14) = {D/2 - r, 0, 0};
Point(15) = {D/2 - r*Sqrt(2)/2, 0, r*Sqrt(2)/2};
Point(16) = {D/2, 0, r};
Point(17) = {D/2, 0, 0};

Point(19) = {D/2, 0, b+r};

Line(18) = {2, 14};
Circle(19) = {14, 17, 15};
Circle(20) = {15, 17, 16};
Line(21) = {16, 19};

Curve Loop(1) = {6, 4, -3, 2};
Plane Surface(1) = {1};

Line(22) = {7, 8};
Line(23) = {8, 5};
Circle(24) = {4, 1, 5};

Curve Loop(2) = {24, -23, -22, -4};
Plane Surface(2) = {2};

Line(25) = {1, 6};
Line(26) = {1, 8};

Curve Loop(4) = {22, -26, 25, 3};
Plane Surface(3) = {4};

Extrude {0, 0, b+r} {
  Surface{1}; Surface{2}; Surface{3};
}


Line(43) = {19, 20};
Line(44) = {15, 20};

Curve Loop(17) = {44, -27, 18, 19};
Plane Surface(16) = {17};

Curve Loop(18) = {44, -43, -21, -20};
Plane Surface(17) = {18};


Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Surface{16}; Surface{17}; 
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Surface{26}; Surface{22}; 
}

Coherence;

Extrude {0, 0, l-r-b} {
  Surface{15,20,24,30,36}; Layers{ {8,6,4}, {0.1,0.4,1} }; Recombine;
}

Coherence;


Transfinite Curve {40, 37, 51, 81, 32, 30, 44, 43, 62, 63} = 6+1;
Transfinite Curve {41, 33, 42, 39, 25, 26, 3, 53, 46, 52, 45, 59, 79, 72, 56, 76, 67, 74, 64, 22} = 8+1;
Transfinite Curve {54, 55, 48, 47, 50, 49} = 6+1;
Transfinite Curve {80, 73, 77, 61, 58, 69, 60, 57, 68} = 6+1;
Transfinite Curve {75, 78, 66, 71, 65, 70} = 6+1;

Transfinite Surface {:};
Transfinite Volume {:};

Mesh.RecombineAll = 1;
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
Mesh.SecondOrderIncomplete = 0;


Physical Surface("u") = {42, 39, 60, 53, 22, 13, 44};
Physical Surface("v") = {50, 14, 45, 18, 28, 56, 33};
Physical Surface("w") = {40, 26, 19, 23, 3};
Physical Surface("load") = {47, 54, 51, 61, 58};
Physical Volume("bulk") = {1:12};
