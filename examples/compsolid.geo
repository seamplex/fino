// this works only with Gmsh >= 2.16.1
// based on gmsh/demo/boolean/compsolid.geo

SetFactory("OpenCASCADE");

Mesh.Algorithm = 6;
Mesh.CharacteristicLengthMax = 0.15;
Mesh.ElementOrder = 2;
Mesh.Optimize = 1;

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

f() = BooleanFragments { Volume{1}; Delete; }{ Volume{2}; Delete; };

Physical Volume("one") = {1};
Physical Volume("two") = {2};
//+
Physical Surface("fixed") = {1};
//+
Physical Surface("load") = {13};
