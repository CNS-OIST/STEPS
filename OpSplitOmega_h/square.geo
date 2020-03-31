// Default element size is 0.03
// To override it, use CLI options -setnumber elt_size VALUE
If (!Exists(elt_size))
    elt_size = 0.03;
EndIf
Point(1) = {0, 0, 0, elt_size};
Point(2) = {1, 0, 0, elt_size};
Point(3) = {1, 1, 0, elt_size};
Point(4) = {0, 1, 0, elt_size};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Physical Surface("comp1") = {6};
Mesh 2;
