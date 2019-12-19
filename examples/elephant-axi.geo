// an structured rectangle which creates an axi-symmetric
// solid cylinder of radius r and length (height) z
r = 0.5;
l = 2;
lc = 1e-2;

SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, r, l};

Physical Curve("bottom", 1) = {1};
Physical Curve("top", 2) = {3};
Physical Surface("bulk", 3)  = {1};

Transfinite Line {1,3} = 1+r/lc;
Transfinite Line {2,4} = 1+l/lc;
Transfinite Surface "*";

Mesh.RecombineAll = 1;
Mesh.ElementOrder = 2;
