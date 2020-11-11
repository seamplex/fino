Merge "holed_plate.geo";

n1 = 12;
n2 = 16;
n3 = 6;
n4 = 4;

Transfinite Curve{1:7,9,12} = n1 + 1;
Transfinite Curve{8,14,11}  = n2 + 1;
Transfinite Curve{13,15,10} = n3 + 1;
Transfinite Surface{:};

Extrude {0, 0, H} {
  Surface{:}; Layers{n4}; Recombine;
}

Mesh.RecombineAll = 1;
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
Mesh.SecondOrderIncomplete = 0;

Physical Surface("u") = {21, 12};
Physical Surface("v") = {7, 15};
Physical Surface("w") = {1, 2, 3, 5, 4};
Physical Surface("load") = {16, 23};
Physical Volume("bulk") = {1,2,3,4,5};
