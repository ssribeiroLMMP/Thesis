// Gmsh project created on Mon Sep 16 15:47:18 2019
Lx = 350e-3;
Ly = 150e-3;
Lz = 700e-6;

Nx = 0.38*400;
Ny = 0.38*150;
Nz = 6;

Point(1) = {0,0,0};
Point(2) = {Lx,0,0};
Point(3) = {Lx,Ly,0};
Point(4) = {0,Ly,0};

Line(1) = {1,2}; Transfinite Line {1} = Nx Using Progression 1;
Line(2) = {2,3}; Transfinite Line {2} = Ny Using Progression 1;
Line(3) = {3,4}; Transfinite Line {3} = Nx Using Progression 1;
Line(4) = {4,1}; Transfinite Line {4} = Ny Using Progression 1;

Line Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1}; Transfinite Surface {1};

Extrude {0, 0, Lz} {
  Surface{1}; Layers{Nz};
}

Physical Surface("Inlet") = {13};
Physical Surface("Outlet") = {21};
Physical Surface("LeftWall") = {25};
Physical Surface("RightWall") = {17};
Physical Surface("TopWall") = {26};
Physical Surface("BottomWall") = {1};
Physical Volume("Volume") = {1};
