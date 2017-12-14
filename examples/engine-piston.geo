Geometry.OCCTargetUnit = "MM";
Merge "engine-piston.stp";

lc = 1.75;

Characteristic Length { 1: 104 } = lc;
Mesh.CharacteristicLengthMin = 0.9*lc;
Mesh.CharacteristicLengthMax = 1.1*lc;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.ElementOrder = 1;
Mesh.Algorithm = 2;
Mesh.Algorithm3D = 1;

Physical Surface("top", 1) = {70};
Physical Surface("ring 1", 2) = {69};
Physical Surface("ring 1 groove", 3) = {68,67,66};
Physical Surface("ring 2", 4) = {65};
Physical Surface("ring 2 groove", 5) = {64,63,62};
Physical Surface("ring 3", 6) = {61};
Physical Surface("ring 3 groove", 7) = {60,55,52,57,54,49};

Physical Surface("interior and skirt", 8) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,50,51,53,56,58,59};

Physical Volume("bulk") = {1};

