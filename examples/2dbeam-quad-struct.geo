SetFactory("OpenCASCADE");
l = 50;
h = 5;
n = 3;
Rectangle(1) = {0, 0, 0, l, h};

Physical Curve("left") = {4};
Physical Curve("bottom") = {1};
Physical Curve("right") = {2};
Physical Curve("top") = {3};
Physical Surface("bulk") = {1};

Transfinite Curve{1,3} = l/h*n+1;
Transfinite Curve{2,4} = n+1;
Transfinite Surface{1};

Mesh.ElementOrder = 1;
Mesh.RecombineAll = 1;
