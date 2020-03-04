// Gmsh project created on Mon Jul 15 18:34:03 2019
// Parâmetros
Lx0 = 0;
Lxf = 0.770;
Ly0 = 0.0;
Lyf = 0.15;

// Define pontos do contorno
p1 = newp;
Point(p1) = {Lx0, Ly0, 0.1};
p2 = newp;
Point(p2) = {Lxf, Ly0, 0.1};
p3 = newp;
Point(p3) = {Lxf, Lyf, 0.1};
p4 = newp;
Point(p4) = {Lx0, Lyf, 0.1};

// Define Linhas do Contorno
l1 = newl;
Line(l1) = {p1, p2};
l2 = newl;
Line(l2) = {p2, p3};
l3 = newl;
Line(l3) = {p3, p4};
l4 = newl;
Line(l4) = {p4, p1};
ll1 = newll;

// Define Loop de linhas
Line Loop(ll1) = {l1, l2, l3, l4};

// Define quantos pontos tem na linha
Transfinite Line{l4} = 180;

// Define superfície do fluido
s0 = news;
Plane Surface(s0) = {ll1};

// Define cada contorno nomeado
Physical Line("Inlet") = l4;
Physical Line("Outlet") = l2;
Physical Line("BottomWall") = l1;
Physical Line("TopWall") = l3;
Physical Surface("Fluid") = s0;

