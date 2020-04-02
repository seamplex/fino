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
Mesh.Algorithm = 6;
Mesh.Optimize = 1;
Mesh.ElementOrder = 2;
Mesh.RecombineAll = 1;
