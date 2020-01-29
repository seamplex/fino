SetFactory("OpenCASCADE");
Box(1) = {0,0,0, 1,1,1};


Transfinite Line "*" = 1;
Transfinite Surface "*";
Transfinite Volume "*";

Mesh.RecombineAll = 1;
Mesh.Recombine3DAll = 1;
Mesh.MshFileVersion = 2.2;

Physical Surface("left") = {1};
Physical Surface("right") = {2};
Physical Volume("bulk")  = {1};
