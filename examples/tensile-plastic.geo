// Sample geo file for meshing the tensile test specimen example

SetFactory("OpenCASCADE");
Merge "tensile-test-specimen.step";

lc = 5;
Mesh.CharacteristicLengthMax = lc;

Mesh.ElementOrder = 2;

Mesh.Algorithm = 6;
Mesh.Algorithm3D = 1;

Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.HighOrderOptimize = 2;


Physical Surface ("left", 1) =  {1};
Physical Surface ("right", 21) =  {7};
Physical Volume ("bulk", 1) =  {1};
