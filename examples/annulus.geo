SetFactory("OpenCASCADE");

Disk(1) = {0,0,0, 2};
Disk(2) = {0,0,0, 1};
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }

Physical Curve("inner") = {2};
Physical Curve("outer") = {3};
Physical Surface("bulk") = {1};

Mesh.CharacteristicLengthMax = 0.1;
