// NAFEMS Benchmark #1: elliptical membrane
SetFactory("OpenCASCADE");
Merge "le1.brep";

Physical Point("A") = {1};
Physical Point("B") = {2};
Physical Point("C") = {3};
Physical Point("D") = {4};

Physical Curve("AB") = {1};
Physical Curve("BC") = {2};
Physical Curve("CD") = {3};
Physical Curve("DA") = {4};

Physical Surface("bulk") = {1};

Mesh.CharacteristicLengthMax = 0.05;

// unstructured grid
Mesh.Algorithm = 8;
Mesh.Optimize = 1;
Mesh.HighOrderOptimize = 2;

// structured grid
Transfinite Curve{1,3} = 2/Mesh.CharacteristicLengthMax;
Transfinite Curve{2,4} = 3/Mesh.CharacteristicLengthMax;
Transfinite Surface{1};

Mesh.ElementOrder = 2;
Mesh.RecombineAll = 1;
