//
Merge "square-beam.step";

lc = 2.5;
struct = 1;
hex = 0;
order = 1;

If ( struct != 0 )
 Transfinite Line {1,2,3,4,5,6,7,8} = 10/lc+1;
 Transfinite Line {9,10,11,12} = 100/lc+1;
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
Mesh.ElementOrder = order;


Physical Volume ("bulk") = 1;

Physical Surface("left") = {1};
Physical Surface("right") = {2};

