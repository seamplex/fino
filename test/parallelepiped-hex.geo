SetFactory("OpenCASCADE");

h = 10;     // base length
l = 20;     // height
Box(1) = {0, 0, 0, l, -h/2, -h/2};
Box(2) = {0, 0, 0, l, -h/2, +h/2};
Box(3) = {0, 0, 0, l, +h/2, +h/2};
Box(4) = {0, 0, 0, l, +h/2, -h/2};
Coherence;


// pretty names
Physical Point("O") = {3};
Physical Point("A") = {7};
Physical Point("B") = {14};
Physical Point("C") = {10};
Physical Point("D") = {15};

Physical Surface("x=0") = {1,7,12,17};
//Physical Surface("z=h/2") = {11,16};
Physical Surface("x=l") = {2,8,13,18};
//Physical Surface("z=-h/2") = {5,20};
//Physical Surface("y=h/2") = {14,19};
//Physical Surface("y=-h/2") = {3,9};

Physical Volume("bulk") = {1,2,3,4};

Transfinite Curve{1,2,3,4,23,29,30,13,14,15,21,22,5,6,7,8,16,17,18,24,25,26,31,32} = 4/Mesh.ElementOrder + 1;
Transfinite Curve{9,10,11,12,19,20,28,27,33} = 16/Mesh.ElementOrder + 1;
Transfinite Surface{:};
Transfinite Volume{:};

Mesh.RecombineAll = 1;
