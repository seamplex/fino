SetFactory("OpenCASCADE");
Geometry.OCCTargetUnit = "MM";
a() = ShapeFromFile("tensile-test-specimen.step");

Mesh.CharacteristicLengthMax = 3;
Mesh.Algorithm = 8;
Mesh.ElementOrder = 2;
Mesh.Optimize = 1;

Physical Surface("left") = {1};
Physical Surface("right") = {7};

Physical Volume("bulk") = {1};
