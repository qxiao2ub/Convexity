% solve 4 Lagrangian eqn_112021
%{
function F = SolveLag(p)

F(1) = (8*p(3)*exp ((8*p(4))/725 - (8*p(3))/725 + 68/145)*exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245)*exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455)*exp (228/125 - (4*p(3))/625))/(1225*(exp ((8*p(4))/725 - (8*p(3))/725 + 68/145) + 1)*(exp ((8*p(1))/1225 - (8*p(3))/1225+ 428/245) + 1)*(exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455) + 1)*(exp (228/125 - (4*p(3))/625) + 1)) - (8*p(3) * exp ((8*p(4))/725 - (8*p(3))/725 + 68/145)* exp ((16*p(1))/1225 - (16*p(3))/1225 + 856/245)*exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455)* exp (228/125 - (4*p(3))/625))/(1225*(exp ((8*p(4))/725 - (8*p(3))/725 + 68/145)+ 1)*(exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245) + 1)^2*(exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455)+ 1)*(exp (228/125 - (4*p(3))/625) + 1));
F(2) = (16*p(3)*exp ((8*p(4))/725 - (8*p(3))/725 + 68/145)*exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245)*exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455)*exp (228/125 - (4*p(3))/625))/(2275*(exp ((8*p(4))/725 - (8*p(3))/725 + 68/145) + 1)*(exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245) + 1)*(exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455) + 1)*(exp (228/125 - (4*p(3))/625) + 1)) - (16*p(3)*exp ((8*p(4))/725 - (8*p(3))/725 + 68/145)*exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245)*exp ((32*p(2))/2275 - (32*p(3))/2275 + 1312/455)*exp (228/125 - (4*p(3))/625))/(2275*(exp ((8*p(4))/725 - (8*p(3))/725 + 68/145) +1)*(exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245) +1)*(exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455) + 1)^2*(exp (228/125 - (4*p(3))/625) + 1));
F(3) = (exp ((8*p(4))/725 - (8*p(3))/725 + 68/145)*exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245)*exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455)*exp (228/125 - (4*p(3))/625))/((exp ((8*p(4))/725 - (8*p(3))/725 + 68/145) + 1)*(exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245) +  1)*(exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455) + 1)*(exp (228/125 - (4*p(3))/625) + 1)) - (357892*p(3)*exp ((8*p(4))/725 - (8*p(3))/725 + 68/145)*exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245)*exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455)*exp (228/125 - (4*p(3))/625))/(11545625*(exp ((8*p(4))/725 - (8*p(3))/725 + 68/145) + 1)*(exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245) + 1)*(exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455) + 1)*(exp (228/125 - (4*p(3))/625) + 1)) + (8*p(3)*exp ((16*p(4))/725 - (16*p(3))/725 + 136/145)*exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245)*exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455)*exp (228/125 - (4*p(3))/625))/(725*(exp ((8*p(4))/725 - (8*p(3))/725 + 68/145) + 1)^2*(exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245) + 1)*(exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455) + 1)*(exp (228/125 - (4*p(3))/625) + 1)) + (4*p(3)*exp ((8*p(4))/725 - (8*p(3))/725 + 68/145)*exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245)*exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455)*exp (456/125 - (8*p(3))/625))/(625*(exp ((8*p(4))/725 - (8*p(3))/725 + 68/145) + 1)*(exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245) + 1)*(exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455) +1)*(exp (228/125 - (4*p(3))/625) + 1)^2) + (8*p(3)*exp ((8*p(4))/725 - (8*p(3))/725 + 68/145)*exp ((16*p(1))/1225 - (16*p(3))/1225 + 856/245)*exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455)*exp (228/125 - (4*p(3))/625))/(1225*(exp ((8*p(4))/725 - (8*p(3))/725 + 68/145) + 1)*(exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245) + 1)^2*(exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455) + 1)*(exp (228/125 - (4*p(3))/625) + 1)) + (16*p(3)*exp ((8*p(4))/725 - (8*p(3))/725 + 68/145)*exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245)*exp ((32*p(2))/2275 - (32*p(3))/2275 + 1312/455)*exp (228/125 - (4*p(3))/625))/(2275*(exp ((8*p(4))/725 - (8*p(3))/725 + 68/145) + 1)*(exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245) + 1)*(exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455) +1)^2*(exp (228/125 - (4*p(3))/625) + 1));
F(4) = (8*p(3)*exp ((8*p(4))/725 - (8*p(3))/725 + 68/145)* exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245)*exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455)*exp (228/125 - (4*p(3))/625))/(725*(exp ((8*p(4))/725 - (8*p(3))/725 + 68/145) + 1)*(exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245) + 1)*(exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455) + 1)*(exp (228/125 - (4*p(3))/625) + 1)) - (8*p(3)*exp ((16*p(4))/725 - (16*p(3))/725 + 136/145)*exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245)*exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455)*exp (228/125 - (4*p(3))/625))/(725*(exp ((8*p(4))/725 - (8*p(3))/725 + 68/145) +1)^2*(exp ((8*p(1))/1225 - (8*p(3))/1225 + 428/245) +1)*(exp ((16*p(2))/2275 - (16*p(3))/2275 + 656/455) +1)*(exp (228/125 - (4*p(3))/625) + 1));

% F(1) = diff(obj,p1);
% F(2) = diff(obj,p2);
% F(3) = diff(obj,p3);
% F(4) = diff(obj,p4);

% following type code:
%}

% take derivative for obj include j_112521
%{
syms p1 p2 p3 p4 j1 j2 j3 j4

k = .8;
alpha0 = .63;
alpha1 = .66;
alpha2 = .91;
alpha3 = 1.73;
alpha4 = 1.56;
L00 = 1;
L11 = .98;
L22 = .91;
L33 = .84;
L44 = .58;
beta = -.004;
p0 = 12;
gamma = 0.03;

p0 = 12;

% original p5/10/15/20
% obj = p15*rdivide(exp(2*k*rdivide((alpha15-alpha0+beta*(p15-p0)+gama*pc15*15),L00)),1+exp(2*k*rdivide((alpha15-alpha0+beta*(p15-p0)+gama*pc15*15),L00)))...
% *rdivide(exp(2*k*rdivide((alpha15-alpha5+beta*(p15-p5)+gama*(pc15*15-pc5*5)),L11)),1+exp(2*k*rdivide((alpha15-alpha5+beta*(p15-p5)+gama*(pc15*15-pc5*5)),L11)))...
% *rdivide(exp(2*k*rdivide((alpha15-alpha10+beta*(p15-p10)+gama*(pc15*15-pc10*10)),L22)),1+exp(2*k*rdivide((alpha15-alpha10+beta*(p15-p10)+gama*(pc15*15-pc10*10)),L22)))...
% *rdivide(exp(2*k*rdivide((alpha15-alpha20+beta*(p15-p20)+gama*(pc15*15-pc20*20)),L44)),1+exp(2*k*rdivide((alpha15-alpha20+beta*(p15-p20)+gama*(pc15*15-pc20*20)),L44)));

obj_1 = p3*rdivide(exp(2*k*rdivide((alpha3-alpha0+beta*(p3-p0)+gamma*(j3)),L00)),1+exp(2*k*rdivide((alpha3-alpha0+beta*(p3-p0)+gamma*(j3)),L00)))...
*rdivide(exp(2*k*rdivide((alpha3-alpha1+beta*(p3-p1)+gamma*(j3-j1)),L11)),1+exp(2*k*rdivide((alpha3-alpha1+beta*(p3-p1)+gamma*(j3-j1)),L11)))...
*rdivide(exp(2*k*rdivide((alpha3-alpha2+beta*(p3-p2)+gamma*(j3-j2)),L22)),1+exp(2*k*rdivide((alpha3-alpha2+beta*(p3-p2)+gamma*(j3-j2)),L22)))...
*rdivide(exp(2*k*rdivide((alpha3-alpha4+beta*(p3-p4)+gamma*(j3-j4)),L44)),1+exp(2*k*rdivide((alpha3-alpha4+beta*(p3-p4)+gamma*(j3-j4)),L44))); %P3

%{
obj = p3*rdivide(exp(2*k*rdivide((alpha15-alpha0+beta*(p3-p0)+gama*(p7*15-p9*0.0001)),L00)),1+exp(2*k*rdivide((alpha15-alpha0+beta*(p3-p0)+gama*(p7*15-p9*0.0001)),L00)))...
*rdivide(exp(2*k*rdivide((alpha15-alpha5+beta*(p3-p1)+gama*(p7*15-p5*5)),L11)),1+exp(2*k*rdivide((alpha15-alpha5+beta*(p3-p1)+gama*(p7*15-p5*5)),L11)))...
*rdivide(exp(2*k*rdivide((alpha15-alpha10+beta*(p3-p2)+gama*(p7*15-p6*10)),L22)),1+exp(2*k*rdivide((alpha15-alpha10+beta*(p3-p2)+gama*(p7*15-p6*10)),L22)))...
*rdivide(exp(2*k*rdivide((alpha15-alpha20+beta*(p7-p4)+gama*(p7*15-p8*20)),L44)),1+exp(2*k*rdivide((alpha15-alpha20+beta*(p3-p4)+gama*(p7*15-p8*20)),L44)));
%}
dif_p1=diff(obj_1,p1)
dif_p2=diff(obj_1,p2)
dif_p3=diff(obj_1,p3)
dif_p4=diff(obj_1,p4)
dif_j1=diff(obj_1,j1)
dif_j2=diff(obj_1,j2)
dif_j3=diff(obj_1,j3)
dif_j4=diff(obj_1,j4)
%}
       
% take derivative for obj_112021
%{
syms p1 p2 p3 p4


k = .8;
alpha0 = .63;
alpha1 = .66;
alpha2 = .91;
alpha3 = 1.73;
alpha4 = 1.56;
L00 = 1;
L11 = .98;
L22 = .91;
L33 = .84;
L44 = .58;
beta = -.004;


p0 = 10;

obj = p3*rdivide(exp(2*k*rdivide((alpha3-alpha0+beta*(p3-p0)),L00)),1+exp(2*k*rdivide((alpha3-alpha0+beta*(p3-p0)),L00)))...
*rdivide(exp(2*k*rdivide((alpha3-alpha1+beta*(p3-p1)),L11)),1+exp(2*k*rdivide((alpha3-alpha1+beta*(p3-p1)),L11)))...
*rdivide(exp(2*k*rdivide((alpha3-alpha2+beta*(p3-p2)),L22)),1+exp(2*k*rdivide((alpha3-alpha2+beta*(p3-p2)),L22)))...
*rdivide(exp(2*k*rdivide((alpha3-alpha4+beta*(p3-p4)),L44)),1+exp(2*k*rdivide((alpha3-alpha4+beta*(p3-p4)),L44)));

diff(obj,p4)
%}

% bundle paper whole obj concave proof-4 bun ops-112522
% %{
syms p1 p2 p3 p4 j1 j2 j3 j4 z1 z2 z3 z4 z5 z6 z7 z8

k = .8;
alpha0 = .63;
alpha1 = .66;
alpha2 = .91;
alpha3 = 1.73;
alpha4 = 1.56;
L00 = 1;
L11 = .98;
L22 = .91;
L33 = .84;
L44 = .58;
beta = -.004;
p0 = 12;
gamma = 0.03;

obj_1 = p3*rdivide(exp(2*k*rdivide((alpha3-alpha0+beta*(p3-p0)+gamma*(j3)),L00)),1+exp(2*k*rdivide((alpha3-alpha0+beta*(p3-p0)+gamma*(j3)),L00)))...
*rdivide(exp(2*k*rdivide((alpha3-alpha1+beta*(p3-p1)+gamma*(j3-j1)),L11)),1+exp(2*k*rdivide((alpha3-alpha1+beta*(p3-p1)+gamma*(j3-j1)),L11)))...
*rdivide(exp(2*k*rdivide((alpha3-alpha2+beta*(p3-p2)+gamma*(j3-j2)),L22)),1+exp(2*k*rdivide((alpha3-alpha2+beta*(p3-p2)+gamma*(j3-j2)),L22)))...
*rdivide(exp(2*k*rdivide((alpha3-alpha4+beta*(p3-p4)+gamma*(j3-j4)),L44)),1+exp(2*k*rdivide((alpha3-alpha4+beta*(p3-p4)+gamma*(j3-j4)),L44))); %P3
%{
+ p1*rdivide(exp(2*k*rdivide((alpha1-alpha0+beta*(p1-p0)),L00)),1+exp(2*k*rdivide((alpha1-alpha0+beta*(p1-p0)),L00)))...
*rdivide(exp(2*k*rdivide((alpha1-alpha2+beta*(p1-p2)),L22)),1+exp(2*k*rdivide((alpha1-alpha2+beta*(p1-p2)),L22)))...
*rdivide(exp(2*k*rdivide((alpha1-alpha3+beta*(p1-p3)),L33)),1+exp(2*k*rdivide((alpha1-alpha3+beta*(p1-p3)),L33)))...
*rdivide(exp(2*k*rdivide((alpha1-alpha4+beta*(p1-p4)),L44)),1+exp(2*k*rdivide((alpha1-alpha4+beta*(p1-p4)),L44)))... %P1
+ p2*rdivide(exp(2*k*rdivide((alpha2-alpha0+beta*(p2-p0)),L00)),1+exp(2*k*rdivide((alpha2-alpha0+beta*(p2-p0)),L00)))...
*rdivide(exp(2*k*rdivide((alpha2-alpha1+beta*(p2-p1)),L11)),1+exp(2*k*rdivide((alpha2-alpha1+beta*(p2-p1)),L11)))...
*rdivide(exp(2*k*rdivide((alpha2-alpha3+beta*(p2-p3)),L33)),1+exp(2*k*rdivide((alpha2-alpha3+beta*(p2-p3)),L33)))...
*rdivide(exp(2*k*rdivide((alpha2-alpha4+beta*(p2-p4)),L44)),1+exp(2*k*rdivide((alpha2-alpha4+beta*(p2-p4)),L44)))... %P2
+ p4*rdivide(exp(2*k*rdivide((alpha4-alpha0+beta*(p4-p0)),L00)),1+exp(2*k*rdivide((alpha4-alpha0+beta*(p4-p0)),L00)))...
*rdivide(exp(2*k*rdivide((alpha4-alpha1+beta*(p4-p1)),L11)),1+exp(2*k*rdivide((alpha4-alpha1+beta*(p4-p1)),L11)))...
*rdivide(exp(2*k*rdivide((alpha4-alpha2+beta*(p4-p2)),L22)),1+exp(2*k*rdivide((alpha4-alpha2+beta*(p4-p2)),L22)))...
*rdivide(exp(2*k*rdivide((alpha4-alpha3+beta*(p4-p3)),L33)),1+exp(2*k*rdivide((alpha4-alpha3+beta*(p4-p3)),L33))); %P4
%}
% + p0*rdivide(exp(2*k*rdivide((alpha0-alpha1+beta*(p0-p1)),L11)),1+exp(2*k*rdivide((alpha0-alpha1+beta*(p0-p1)),L11)))...
% *rdivide(exp(2*k*rdivide((alpha0-alpha2+beta*(p0-p2)),L22)),1+exp(2*k*rdivide((alpha0-alpha2+beta*(p0-p2)),L22)))...
% *rdivide(exp(2*k*rdivide((alpha0-alpha3+beta*(p0-p3)),L33)),1+exp(2*k*rdivide((alpha0-alpha3+beta*(p0-p3)),L33)))...
% *rdivide(exp(2*k*rdivide((alpha0-alpha4+beta*(p0-p4)),L44)),1+exp(2*k*rdivide((alpha0-alpha4+beta*(p0-p4)),L44)))... %P0

result = hessian(obj_1,[p1,p2,p3,p4,j1,j2,j3,j4]);

B=[z1
  z2
  z3
  z4
  z5
  z6
  z7
  z8
  ];

semidef_1 = transpose(B)*result*B;

semidef_sim_1 = simplify(semidef_1);

% save long semi-definite result to txt
% %{
file1= fopen('semidef_objccv.txt','wt');
fprintf(file1,'%s\n',semidef_1);
fclose(file1);
%}

% save long simplified semi-definite result to txt
% %{
file1= fopen('semidef_sim_objccv.txt','wt');
fprintf(file1,'%s\n',semidef_sim_1);
fclose(file1);

%}

% convex proof only for P3_111921
 %{
syms p0 p1 p2 p3 p4 z0 z1 z2 z3 z4

k = .8;
alpha0 = .63;
alpha1 = .66;
alpha2 = .91;
alpha3 = 1.73;
alpha4 = 1.56;
L00 = 1;
L11 = .98;
L22 = .91;
L33 = .84;
L44 = .58;
beta = -.004;

obj = p3*rdivide(exp(2*k*rdivide((alpha3-alpha0+beta*(p3-p0)),L00)),1+exp(2*k*rdivide((alpha3-alpha0+beta*(p3-p0)),L00)))...
*rdivide(exp(2*k*rdivide((alpha3-alpha1+beta*(p3-p1)),L11)),1+exp(2*k*rdivide((alpha3-alpha1+beta*(p3-p1)),L11)))...
*rdivide(exp(2*k*rdivide((alpha3-alpha2+beta*(p3-p2)),L22)),1+exp(2*k*rdivide((alpha3-alpha2+beta*(p3-p2)),L22)))...
*rdivide(exp(2*k*rdivide((alpha3-alpha4+beta*(p3-p4)),L44)),1+exp(2*k*rdivide((alpha3-alpha4+beta*(p3-p4)),L44)));

result = hessian(obj,[p0,p1,p2,p3,p4]);

A=[z0
  z1
  z2
  z3
  z4
  ];

semidef = transpose(A)*result*A;

semidef_sim = simplify(semidef);

%}


% coding for https://www.youtube.com/watch?v=F3YoC5A6Avg 
%{
syms x1 x2 x3 z1 z2 z3
f = (x1-x2)^2+2*x3^2;
result_ex = hessian(f,[x1,x2,x3]);

B=[z1
  z2
  z3];

semidef_ex = transpose(B)*result_ex*B;

%}
