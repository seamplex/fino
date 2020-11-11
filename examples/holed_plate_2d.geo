Merge "holed_plate.geo";

n1 = 24;
n2 = 32;
n3 = 12;

Transfinite Curve{1:7,9,12} = n1 + 1;
Transfinite Curve{8,14,11}  = n2 + 1;
Transfinite Curve{13,15,10} = n3 + 1;
Transfinite Surface{:};

Mesh.RecombineAll = 1;
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
Mesh.SecondOrderIncomplete = 0;

Physical Curve("u") = {6, 13};
Physical Curve("v") = {3, 8};
Physical Curve("load") = {10, 9};
Physical Surface("bulk") = {1, 2, 3, 4, 5};
