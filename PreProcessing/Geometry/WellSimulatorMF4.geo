// Gmsh project created on Mon Jul 15 18:34:03 2019
Lz0 = 6;
Lzf = 8;
Lr0 = 0.10;
Lrf = 0.22;
FLzf = 7.1;
FLz0 = 7.0;
MeshFactor = 4.0;
p1 = newp;
Point(p1) = {Lz0, Lr0, 1};
p2 = newp;
Point(p2) = {FLzf, Lr0, 1};
p3 = newp;
Point(p3) = {Lzf, Lr0, 1};
p4 = newp;
Point(p4) = {Lzf, Lrf, 1};
p5 = newp;
Point(p5) = {FLzf, Lrf, 1};
p6 = newp;
Point(p6) = {FLz0, Lrf, 1};
p7 = newp;
Point(p7) = {Lz0, Lrf, 1};
l1 = newl;
Line(l1) = {p1, p2}; Transfinite Line{l1} = 10*MeshFactor;
l2 = newl;
Line(l2) = {p2, p3}; Transfinite Line{l2} = 8*MeshFactor;
l3 = newl;
Line(l3) = {p3, p4}; Transfinite Line{l3} = 2*MeshFactor;
l4 = newl;
Line(l4) = {p4, p5}; Transfinite Line{l4} = 8*MeshFactor;
l5 = newl;
Line(l5) = {p5, p6};  Transfinite Line{l5} = 4*MeshFactor;
l6 = newl;
Line(l6) = {p6, p7};  Transfinite Line{l6} = 9*MeshFactor; 
l7 = newl;
Line(l7) = {p7, p1};  Transfinite Line{l7} = 4*MeshFactor;
ll1 = newll;
Line Loop(ll1) = {l1, l2, l3, l4, l5, l6, l7};
s1 = news;
Plane Surface(s1) = {ll1};
Physical Line("Inlet") = l7;
Physical Line("Outlet") = l5;
Physical Line("InnerPipe") = {l1,l2};
Physical Line("OuterWall") = {l4,l6};
Physical Line("BottomWall") = l3;
Physical Surface("Fluid") = s1;
