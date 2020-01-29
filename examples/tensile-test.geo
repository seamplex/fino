Merge "tensile-test-specimen.step"; // read the step file
Mesh.CharacteristicLengthMax = 1.5; // set the max element size lc
Mesh.ElementOrder = 2;              // ask for second-order elements

// define physical entities for BCs and materials
// the name in the LHS has to appear in the Fino input
// the number in the RHS is the numerical id of the entity in the CAD file
Physical Surface ("left") =  {1};   // left face, to be fixed
Physical Surface ("right") =  {7};  // right face, will hold the load
Physical Volume ("bulk") =  {1};    // bulk material elements
