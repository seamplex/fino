a = 1;
n = 4;

SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, a, a, a};

Physical Surface("left") = {1};
Physical Surface("right") = {2};
Physical Surface("front") = {3};
Physical Surface("back") = {4};
Physical Surface("bottom") = {5};
Physical Surface("top") = {6};
Physical Volume("bulk") = {1};

Transfinite Line "*" = n;
Transfinite Surface "*";
Transfinite Volume "*";

//Mesh.Algorithm = 6;
Mesh.ElementOrder = 2;
Mesh.SecondOrderIncomplete = 1;

Mesh.RecombineAll = 1;
Mesh.Recombine3DAll = 1;
Mesh.Recombine3DLevel = 2;

Mesh.Optimize = 1;

