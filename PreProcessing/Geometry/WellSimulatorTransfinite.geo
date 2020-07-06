// Gmsh project created on Mon Jul 15 18:34:03 2019
Lz0 = 6;
Lzf = 8;
Lr0 = 0.10;
Lrf = 0.22;
FLzf = 7.1;
FLz0 = 7.0;
p1 = newp;
Point(p1) = {Lz0, Lr0, 1};
p2 = newp;
Point(p2) = {FLz0, Lr0, 1};
p3 = newp;
Point(p3) = {FLzf, Lr0, 1};
p4 = newp;
Point(p4) = {Lzf, Lr0, 1};
p5 = newp;
Point(p5) = {Lzf, Lrf, 1};
p6 = newp;
Point(p6) = {FLzf, Lrf, 1};
p7 = newp;
Point(p7) = {FLz0, Lrf, 1};
p8 = newp;
Point(p8) = {Lz0, Lrf, 1};
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
Line(l6) = {p6, p7};   
l7 = newl;
Line(l7) = {p7, p8};  
l8 = newl;
Line(l8) = {p8, p1};  
ll1 = newll;
Line Loop(ll1) = {l1, l2, l3, l4, l5, l6, l7, l8};
s1 = news;
Plane Surface(s1) = {ll1}; 
Physical Line("Inlet") = l8;
Physical Line("Outlet") = l6;
Physical Line("InnerPipe") = {l1,l2,l3};
Physical Line("OuterWall") = {l5,l7};
Physical Line("BottomWall") = l4;
Physical Surface("Fluid") = s1;
//+
Transfinite Surface {10} = {8, 5, 4, 1};
//+
Transfinite Line {8, 4} = 20 Using Progression 1;
//+
Transfinite Line {7, 1} = 50 Using Progression 1;
//+
Transfinite Line {5, 3} = 45 Using Progression 1;
//+
Transfinite Line {6, 2} = 20 Using Progression 1;
//+
Recombine Surface {10};
