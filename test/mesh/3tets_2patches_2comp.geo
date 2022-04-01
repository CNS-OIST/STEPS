SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthFactor = 10;

Mesh.MshFileVersion = 4.1;
Mesh.PartitionOldStyleMsh2 = 1;
Mesh.PartitionCreateGhostCells = 1;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(4) = {0, 0, 1, 1.0};
Point(5) = {-1, 0, 0, 1.0};
Point(6) = {0, -1, 0, 1.0};

Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 1};
Line(4) = {1, 5};
Line(5) = {5, 3};
Line(6) = {1, 2};
Line(7) = {2, 3};
Line(8) = {5, 4};
Line(9) = {4, 2};
Line(10) = {4, 6};
Line(11) = {6, 1};
Line(12) = {6, 2};

Curve Loop(1) = {8, -1, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {1, 9, -6};
Plane Surface(2) = {2};
Curve Loop(3) = {6, 7, 3};
Plane Surface(3) = {3};
Curve Loop(4) = {3, 4, 5};
Plane Surface(4) = {4};
Curve Loop(5) = {5, -2, -8};
Plane Surface(5) = {5};
Curve Loop(6) = {2, -7, -9};
Curve Loop(7) = {1, 2, 3};
Plane Surface(6) = {7};
Curve Loop(8) = {2, -7, -9};
Plane Surface(7) = {8};
Curve Loop(9) = {10, 11, 1};
Plane Surface(8) = {9};
Curve Loop(10) = {10, 12, -9};
Plane Surface(9) = {10};
Curve Loop(11) = {12, -6, -11};
Plane Surface(10) = {11};

Surface Loop(1) = {5, 1, 4, 3, 7, 2};
Volume(1) = {1};
Surface Loop(2) = {9, 8, 10, 2};
Volume(2) = {2};

Physical Surface("patch1", 1) = {4};
Physical Surface("patchInBetween", 2) = {2};
Physical Volume("comp1", 3) = {1};
Physical Volume("comp2", 4) = {2};
