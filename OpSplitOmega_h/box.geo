// Default element size is 0.1
// To override it, use CLI options -setnumber elt_size VALUE
If (!Exists(elt_size))
    elt_size = .1;
EndIf
Point(1) = {0, 0, 0, elt_size};
Point(2) = {1, 0, 0, elt_size};
Point(3) = {1, 1, 0, elt_size};
Point(4) = {0, 1, 0, elt_size};
Point(5) = {0, 0, 1, elt_size};
Point(6) = {1, 0, 1, elt_size};
Point(7) = {1, 1, 1, elt_size};
Point(8) = {0, 1, 1, elt_size};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line(5) = {1, 5};
Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {4, 8};
Line(9) = {5, 6};
Line(10) = {6, 7};
Line(11) = {7, 8};
Line(12) = {8, 5};
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(13) = {13};
Line Loop(14) = {-1, 5, -12, -8};
Plane Surface(14) = {14};
Line Loop(15) = {-2, 8, -11, -7};
Plane Surface(15) = {15};
Line Loop(16) = {-3, 7, -10, -6};
Plane Surface(16) = {16};
Line Loop(17) = {-4, 6, -9, -5};
Plane Surface(17) = {17};
Line Loop(18) = {9, 10, 11, 12};
Plane Surface(18) = {18};
Surface Loop(19) = {13, 14, 15, 16, 17, 18};
Volume(19) = {19};
Physical Volume("comp1") = {19};
