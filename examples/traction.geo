//
Merge "square-beam.step";

lc = 1;
struct = 1;
hex = 1;

If ( struct != 0 )
 Transfinite Line {1,2,3,4,5,6,7,8} = 10/lc;
 Transfinite Line {9,10,11,12} = 100/lc;
 Transfinite Surface "*";
 Transfinite Volume "*";
Else
 Characteristic Length {1:8} = lc;
EndIf

If ( hex != 0 ) 
 Mesh.RecombineAll = 1;
 If ( struct == 0 )
   Mesh.Recombine3DAll = 1;
   Mesh.Recombine3DLevel = 2;
 EndIf
EndIf

Mesh.Optimize = 1;


Physical Volume ("bulk") = 1;

Physical Surface("left") = {1};
Physical Surface("right") = {2};
