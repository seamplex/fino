// based on gmsh/demo/boolean/compsolid.geo
SetFactory("OpenCASCADE");

x = 0;
y = 0;
z = 0;
dx = 2;
dy = 2;
dz = 2;
x2 = x+dx;
y2 = 0;
z2 = 0;
dx2 = 1;
dy2 = 1;
dz2 = 3;

Block(1) = {x,y,z, x+dx,y+dy,z+dz};
Block(2) = {x2,y2,z2, x2+dx2,y2+dy2,z2+dz2};

Coherence;

Physical Volume("one") = {1};
Physical Volume("two") = {2};
//+
Physical Surface("fixed") = {1};
//+
Physical Surface("load") = {13};

lc = 0.15;
Mesh.CharacteristicLengthMax = lc;

Mesh.Algorithm = 6;
Mesh.ElementOrder = 2;

Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.HighOrderOptimize = 2;


// refinements
Field[1] = Distance;
Field[1].EdgesList = {7};
Field[1].NNodesByEdge = 50;
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 0.3*lc;
Field[2].LcMax = lc;
Field[2].DistMin = 1 * lc;
Field[2].DistMax = 3 * lc;

Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = {3};
