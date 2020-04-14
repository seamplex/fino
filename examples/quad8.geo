Merge "quad.geo";
Transfinite Line {1:4} = 1;      // make sure we have only one element
Transfinite Surface "*";
Mesh.ElementOrder = 2;           // ask for second-order...
Mesh.SecondOrderIncomplete = 1;  // ...incomplete elements (i.e. quad8)
