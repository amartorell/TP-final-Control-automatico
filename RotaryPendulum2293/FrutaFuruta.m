% datos del sistema, joint 1
m1 =0.01835; %0.0198; 
g=9.81;
T=725.926e-3;


aw1= 5.71e4; %sacado del grafico con link 1 bien armado
j1= 1/(aw1);   %obtenemos el momento de inercia con respecto al punto de pivote
L1=0.1035;
l1 = j1/(m1*g*((T/(2*pi))^2));
%cosntantes generales
b1=1e-6;
b2=1*1e-6;

% datos del sistema, joint 2
m2 = 0.00575; %3.19e-3;
% experimento d, sacamos la distancia al centro de masa del joint 2
l2=(2e-3)/(m2*9.81*sin(0.6592)); 
% experimento c, sacamos el periodo del pendulo obtenemos inercia de 2
j2=((T/(2*pi))^2)*m2*g*l2;

%linearizacion del modelo

J1= j1;
J2 = j2 + m2*l2^2;
J0= J1 + m2*L1;

m2l2 = m2*l2;       %producto utilizado varias veces

denom = (J0*J2-(m2l2*L1)^2);

a31=0;
a32= (g*m2l2^2*L1)/denom;
a33= (-b1*J2)/denom;
a34= (-b2*m2l2*L1)/denom;

a41=0;
a42= (g*m2l2*J0)/denom;
a43= (-b1*m2l2*L1)/denom;
a44=(-b2*J0)/denom;

b31=J2/denom;
b41=m2l2*L1/denom;
b32=m2l2*L1/denom;
b42=J0/denom;

A= [ 0 0 1 0;
    0 0 0 1;
    a31 a32 a33 a34;
    a41 a42 a43 a44];


B=[0 0;
   0 0;
   b31 b32;
   b41 b42];

B1=[0; 0 ; b31 ; b41]/1000;
   
C= [ 1 0 0 0];


rank(ctrb(A,B)); %cheque controlabilidad del sistema
rank(obsv(A,C)); %chequeo observavilidad del sistema

sys=ss(A,B1,C,0); 
Q=diag([1 8 1 8 0.05]);
Q2=diag([1 10 1 10]);
R=1;
[K,S,E]=lqi(sys,Q,R);
[K1,S1,E2]=lqr(sys,Q2,R);
Kp=K(1:4);
ki=K(5);

Aclp = (A-B1*Kp);
nsys=ss(Aclp,B1,C,0);