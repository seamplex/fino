//
SetFactory("OpenCASCADE");
l = 0.5*303e-3;
d = 1.948e-3;

n = 12;
Rectangle (1) = {0, 0, 0, d/2, l};
Extrude{ {0,1,0}, {0,0,0}, Pi }{ Surface{1}; Layers{n/2}; Recombine; }
Extrude{ {0,-1,0}, {0,0,0}, Pi }{ Surface{1}; Layers{n/2}; Recombine; }
Coherence;

Mesh.Optimize = 1;
Mesh.RecombineAll = 1;
Mesh.Algorithm = 8;
Mesh.ElementOrder = 2;
Mesh.SecondOrderIncomplete = 1;
Mesh.CharacteristicLengthMax = d/3;

Physical Volume ("bulk") = {1,2};
Physical Surface("fixed") = {2,6};
