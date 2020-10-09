// NAFEMS LE10 benchmark geometry & mesh
SetFactory("OpenCASCADE");

// geometry parameters
a = 1000;
b = 2750;
c = 3250;
d = 2000;
h = 600;

// define the four points
Point(1) = {0, a, 0};
Point(2) = {0, b, 0};
Point(3) = {c, 0, 0};
Point(4) = {d, 0, 0};

// join them with the ellipses and straight edges
Line(1) = {1, 2};
Ellipse (2) = {0,0,0, c, b, 0, Pi/2};
Line(3) = {3, 4};
Ellipse (4) = {0,0,0, d, a, 0, Pi/2};

// merge the points
Coherence;

// create the surface
Curve Loop(1) = {1, -2, 3, 4};
Plane Surface(1) = {1};

// extrude (twice to get the mid-plane edge
Extrude {0, 0, 0.5*h} {
  Surface{1};
}

Extrude {0, 0, 0.5*h} {
  Surface{6};
}

// define physical groups
Physical Surface("DCD'C'") = {9, 4};
Physical Surface("ABA'B'") = {7, 2};
Physical Surface("BCB'C'") = {8, 3};
Physical Surface("upper") = {11};
Physical Curve("midplane") = {9};
Physical Volume("bulk") = {1,2};

Mesh.ElementOrder = 2;
Mesh.SecondOrderIncomplete = 1;
Mesh.RecombineAll = 1;

// the original 6x4x2 NAFEMS "fine" mesh is obtained with CharacteristicLengthFactor = 1
Transfinite Curve{2,9,17,4,12,20} = 6/Mesh.CharacteristicLengthFactor + 1;
Transfinite Curve{1,7,15,3,11,19} = 4/Mesh.CharacteristicLengthFactor + 1;
Transfinite Curve{5,13,6,14,8,16,10,18} = 1/Mesh.CharacteristicLengthFactor + 1;

Transfinite Surface {1,2,3,4,5,6,7,8,9,10,11};
Transfinite Volume {1,2};

