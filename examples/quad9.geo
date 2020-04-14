Merge "quad.geo";
Transfinite Line {1:4} = 1;      // make sure we have only one element
Transfinite Surface "*";
Mesh.ElementOrder = 2;           // ask for second-order...
Mesh.SecondOrderIncomplete = 0;  // ...complete elements (i.e. quad9)
