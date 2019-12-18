// an structured rectangle which creates an axi-symmetric
// solid cylinder of radius r and length (height) z
r = 0.5;
l = 2;
lc = 1e-2;

SetFactory("OpenCASCADE");
Point(1) = {0, 0, 0};
Extrude {r, 0, 0} {
  Point{1}; Layers{1+(r/lc)}; Recombine;
}
Extrude {0, l, 0} {
  Curve{1}; Layers{1+(l/lc)}; Recombine;
}

Physical Curve("bottom", 1) = {1};
Physical Curve("top", 2) = {4};
Physical Surface("bulk", 3)  = {1};

Mesh.ElementOrder = 2;
