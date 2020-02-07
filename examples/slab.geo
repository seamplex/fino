lc = 0.5;
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};

Line(1) = {1, 2};

Physical Point("left") = {1};
Physical Point("right") = {2};
Physical Line("bulk") = {1};
