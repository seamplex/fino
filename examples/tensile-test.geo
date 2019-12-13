SetFactory("OpenCASCADE");
Merge "tensile-test-specimen.step";

lc = 1.5;
Mesh.CharacteristicLengthMax = lc;

Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;

Mesh.Algorithm = 6;
Mesh.Algorithm3D = 1;

Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.HighOrderOptimize = 2;

Mesh.ElementOrder = 2;

// bcs
Physical Surface ("left", 1) =  {1};
Physical Surface ("right", 21) =  {7};
Physical Volume ("bulk", 1) =  {1};
