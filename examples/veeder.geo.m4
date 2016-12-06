//
//lc = 1;
//r = 15;
//H = 15;
h=H/2;


Point(1) = {0, 0, 0, lc};
Point(2) = {r, 0, 0, lc};
Point(3) = {-r, 0, 0, lc};
Point(4) = {0, r, 0, lc};
Point(5) = {0, -r, 0, lc};

Circle(1) = {3, 1, 4};
Circle(2) = {4, 1, 2};
Circle(3) = {2, 1, 5};
Circle(4) = {5, 1, 3};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// embebemos el centro en la superficie
Point{1} In Surface{6};

Extrude {0, 0, h} {
  Surface{6};
}
// embebemos el centro en la superficie
Point{7} In Surface{28};

Physical Surface("inf") = {6};
Physical Surface("sup") = {28};
Physical Surface("lat") = {27, 23, 19, 15};

Physical Volume("bulk") = {1};
Physical Point("origin") = {1};
Physical Point("up") = {7};
