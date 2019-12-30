//+
SetFactory("OpenCASCADE");

a = 1;
l = 5;
lc = a/10;
Box(1) = {0, -0.5*a, -0.5*a, l,  a, a};

Mesh.Algorithm = 6;
Mesh.ElementOrder = 2;
Mesh.SecondOrderIncomplete = 0;

Physical Surface("fixed") = {1};
Physical Surface("load") = {2};
Physical Volume("bulk") = {1};


Transfinite Line {1,2,3,4,5,6,7,8} = a/lc;
Transfinite Line {9,10,11,12} = l/lc;
Transfinite Surface "*";
Transfinite Volume "*";

Mesh.RecombineAll = 1;
Mesh.Recombine3DAll = 1;
Mesh.Recombine3DLevel = 2;
