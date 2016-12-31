//
Merge "two-material-beam.brep";

Transfinite Line {1,2,3,4,6,9,11,12,14,17,19,20} = 10;
Transfinite Line {5,7,8,10,13,15,16,18} = 50;

Transfinite Surface "*";
Transfinite Volume "*";

Mesh.RecombineAll = 1;

Physical Volume ("left") = 1;
Physical Volume ("right") = 2;

Physical Surface("fixed_left") = {1};
Physical Surface("fixed_right") = {11};
Physical Surface("load") = {3, 9};
