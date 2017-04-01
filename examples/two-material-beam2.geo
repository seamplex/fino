//
Merge "two-material-beam.brep";

n = 3;
Transfinite Line {1,2,3,4,6,9,11,12,14,17,19,20} = n;
Transfinite Line {5,7,8,10,13,15,16,18} = 4*n;

Transfinite Surface "*";
Transfinite Volume "*";

Mesh.ElementOrder = 2;

Physical Volume ("left") = 1;
Physical Volume ("right") = 2;

Physical Surface("fixed_left") = {1};
Physical Surface("fixed_right") = {11};
Physical Surface("load") = {3, 9};
