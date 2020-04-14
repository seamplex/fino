// create the continuous rectangle
SetFactory("OpenCASCADE");
Rectangle(1) = {-1/2, -1/2, 0, 1, 1};

// name the edges to set Neumann BCs
Physical Curve("bottom") = {1};
Physical Curve("top") = {3};
Physical Curve("left") = {4};
Physical Curve("right") = {2};

// name the corner points to set Dirichlet BCs
Physical Point("one") = {1};
Physical Point("two") = {2};
Physical Point("three") = {3};
Physical Point("four") = {4};

// one physical surface for the bulk of the material
Physical Surface("bulk", 3)  = {1};

Mesh.RecombineAll = 1;       // ask for quads instead of triangs
