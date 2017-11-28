//+
SetFactory("OpenCASCADE");
//+
a = 1;
l = 10;
lc = a/6;
Box(1) = {0, -0.5*a, -0.5*a, l,  a, a};

// Characteristic Length {1:8} = a/4;
// Mesh.CharacteristicLengthMin = a/4;
// Mesh.CharacteristicLengthMax = a/4;

Mesh.Algorithm = 6;
Mesh.ElementOrder = 1;

//+
Physical Surface("fixed") = {1};
//+
Physical Line("uno") = {6};
//+
Physical Line("dos") = {8};
//+
Physical Surface("load") = {2};
//+
Physical Volume("bulk") = {1};


Transfinite Line {1,2,3,4,5,6,7,8} = a/lc;
Transfinite Line {9,10,11,12} = l/lc;
Transfinite Surface "*";
Transfinite Volume "*";

Mesh.RecombineAll = 1;
Mesh.Recombine3DAll = 1;
Mesh.Recombine3DLevel = 2;

Mesh.Optimize = 1;

