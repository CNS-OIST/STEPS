SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthFactor = 4;
Mesh.MshFileVersion = 2.2;

Sphere(1) = {0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};

Physical Surface("patch1", 1) = {1};
Physical Volume("comp1", 2) = {1};