//
SetFactory("OpenCASCADE");
Cylinder(1) = {0, 0, 0, 1, 0, 0, 0.5};
Physical Volume("bulk") = {1};
Physical Surface("hot") = {3};
Physical Surface("cool") = {1, 2};

lc = 0.1;
Mesh.CharacteristicLengthMax = lc;

// refinements
Field[1] = Distance;
Field[1].FacesList = {3};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 0.5*lc;
Field[2].LcMax = lc;
Field[2].DistMin = 1 * lc;
Field[2].DistMax = 3 * lc;

Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = {3};


