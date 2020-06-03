// NAFEMS LE1 plane-stress benchmark geometry & mesh
SetFactory("OpenCASCADE");

// create the geometry
a = 1.00;
b = 2.75;
c = 3.25;
d = 2.00;

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

// define physical groups
Physical Point("A",1) = {1};
Physical Point("B",2) = {2};
Physical Point("C",3) = {3};
Physical Point("D",4) = {4};
Physical Curve("AB",5) = {1};
Physical Curve("BC",6) = {2};
Physical Curve("CD",7) = {3};
Physical Curve("DA",8) = {4};
Physical Surface("bulk",9   ) = {1};

// ask for second-order curved complete elements (i.e. quad9)
Mesh.ElementOrder = 2;
Mesh.RecombineAll = 1;
Mesh.SecondOrderIncomplete = 0;
Mesh.SecondOrderLinear = 0;
Mesh.CharacteristicLengthMax = 0.04;

// use transfinite algoritmh to obtain a structured grid
Transfinite Curve{1,3} = 2/(Mesh.CharacteristicLengthMax * Mesh.CharacteristicLengthFactor);
Transfinite Curve{2,4} = 3/(Mesh.CharacteristicLengthMax * Mesh.CharacteristicLengthFactor);
Transfinite Surface{1};
