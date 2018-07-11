//
SetFactory("OpenCASCADE");

Mesh.Algorithm = 6;
Mesh.ElementOrder = 2;

t = 1;
R = 20;
c = 1;

Point(1) = {0,0,0, 0.25*c*t};
Point(2) = {0,0,t, 0.25*c*t};
Cylinder(3) = {0,0,0, 0,0,t, R};

Point { 1 } In Surface { 3 };
Point { 2 } In Surface { 2 };

Mesh.CharacteristicLengthMax = c*t;

Physical Point("00") = {1};
Physical Point("01") = {2};
Physical Point("10") = {4};
Physical Point("11") = {3};
//+
Physical Volume("bulk") = {3};
