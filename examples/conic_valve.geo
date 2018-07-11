//
SetFactory("OpenCASCADE");
Geometry.OCCTargetUnit = "MM";
a() = ShapeFromFile("conic_valve.brep");

Mesh.CharacteristicLengthMax = 0.25;
Mesh.Algorithm = 6;
Mesh.ElementOrder = 2;

Physical Volume("bulk") = {1};
Physical Surface("base") = {1};
Physical Surface("cone") = {3};
Physical Surface("top") = {6};
Physical Surface("other") = {2, 4, 5};
