//
SetFactory("OpenCASCADE");

// Geometry.OCCTargetUnit = "MM";
// a() = ShapeFromFile("square-beam.step");

Box(1) = {0, -h/2, -h/2, l, h, h};


If ( struct != 0 )
 Transfinite Line {1:8} = h/lc;
 Transfinite Line {9:12} = l/lc;
 Transfinite Surface "*";
 Transfinite Volume "*";
 Mesh.RecombineAll = 1;
 Mesh.SecondOrderIncomplete = 0;
Else
 Mesh.CharacteristicLengthMax = lc;
EndIf

Mesh.Optimize = 1;
Mesh.Algorithm = 6;
Mesh.ElementOrder = order;

Physical Volume ("bulk") = 1;

Physical Surface("left") = {1};
Physical Surface("right") = {2};
