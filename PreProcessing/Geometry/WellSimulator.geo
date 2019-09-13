// Gmsh project created on Mon Jul 15 18:34:03 2019
Lz0 = 6;
Lzf = 8;
Lr0 = 0.10;
Lrf = 0.22;
FLzf = 7.1;
FLz0 = 7.0;
p1 = newp;
Point(p1) = {Lz0, Lr0, 0.01};
p2 = newp;
Point(p2) = {Lzf, Lr0, 0.01};
p3 = newp;
Point(p3) = {Lzf, Lrf, 0.01};
p4 = newp;
Point(p4) = {FLzf, Lrf, 0.01};
p5 = newp;
Point(p5) = {FLz0, Lrf, 0.01};
p6 = newp;
Point(p6) = {Lz0, Lrf, 0.01};
p7 = newp;
Point(p7) = {0.2*(Lzf-Lz0)/4 +Lz0, (Lrf-Lr0)/2+Lr0, 0.01};
p8 = newp;
Point(p8) = {3.8*(Lzf-Lz0)/4 +Lz0, (Lrf-Lr0)/2+Lr0, 0.01};
l1 = newl;
Line(l1) = {p1, p2};
l2 = newl;
Line(l2) = {p2, p3};
l3 = newl;
Line(l3) = {p3, p4};
l4 = newl;
Line(l4) = {p4, p5};
l5 = newl;
Line(l5) = {p5, p6};
l6 = newl;
Line(l6) = {p6, p1};
l7 = newl;
Line(l7) = {p7, p8};
ll1 = newll;
Line Loop(ll1) = {l1, l2, l3, l4, l5, l6};
s1 = news;
Plane Surface(s1) = {ll1};
Physical Line("Inlet") = l6;
Physical Line("Outlet") = l4;
Physical Line("InnerPipe") = l1;
Physical Line("OuterWall") = {l3,l5};
Physical Line("BottomWall") = l2;
Physical Surface("Fluid") = s1;
Point{p7} In Surface {s1};
Point{p8} In Surface {s1};
Line{l7} In Surface {s1};
//+
Characteristic Length {4, 5} = 0.001;
//+
Characteristic Length {6, 1} = 0.001;
//+
Characteristic Length {1, 2} = 0.004;
//+
Characteristic Length {3, 4} = 0.004;
//+
Characteristic Length {7, 8} = 0.008;
