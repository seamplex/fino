SetFactory("OpenCASCADE");
Geometry.OCCTargetUnit = "MM";
a() = ShapeFromFile("tensile-test-specimen.step");

Mesh.CharacteristicLengthMax = 3;
Mesh.Algorithm = 6;
Mesh.ElementOrder = 2;

Physical Surface("left") = {1};
Physical Surface("right") = {7};

Physical Volume("bulk") = {1};
