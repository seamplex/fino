SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 50, 5};

Physical Curve("left") = {4};
Physical Curve("bottom") = {1};
Physical Curve("right") = {2};
Physical Curve("top") = {3};
Physical Surface("bulk") = {1};

n = 2;
Transfinite Curve{1,3} = 8*n+1;
Transfinite Curve{2,4} = n+1;
Transfinite Surface{1};

Mesh.ElementOrder = 1;
Mesh.RecombineAll = 1;
