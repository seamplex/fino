SetFactory("OpenCASCADE");
Rectangle(1) = {-1/2, -1/2, 0, 1, 1};

Physical Point("one") = {1};
Physical Point("two") = {2};
Physical Point("three") = {3};
Physical Point("four") = {4};

Physical Curve("bottom") = {1};
Physical Curve("top") = {3};
Physical Curve("left") = {4};
Physical Curve("right") = {2};

Physical Surface("bulk", 3)  = {1};

Transfinite Line {1,2,3,4} = 1;
Transfinite Surface "*";

Mesh.RecombineAll = 1;
Mesh.MshFileVersion = 2.2;
