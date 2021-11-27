% clear all; close all;
% controller_on_off=1;
% datos del sistema, joint 1
clear all;close all;
controller_on_off=1;
m1 =0.01835; %0.0198; 
g=9.81;
T=725.926e-3;

tau_1=1;
aw1= 5.71e4; %sacado del grafico con link 1 bien armado
j1= 2*pi*1/(aw1);   %obtenemos el momento de inercia con respecto al punto de pivote
L1=0.1035;
l1 = j1/(m1*g*((T/(2*pi))^2));
%cosntantes generales
b1=1e-6;
b2=0.001e-3; %%1*1e-6;

% datos del sistema, joint 2
m2 = 0.00575; %3.19e-3;
% experimento d, sacamos la distancia al centro de masa del joint 2
l2=(2e-3)/(m2*9.81*sin(0.6592)); 
% experimento c, sacamos el periodo del pendulo obtenemos inercia de 2
j2=((T/(2*pi))^2)*m2*g*l2;

%linearizacion del modelo

J1= j1;
J2 = j2;% + m2*l2^2;
J0= J1 + m2*L1^2;

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

A = [ 0 0 1 0;
    0 0 0 1;
    a31 a32 a33 a34;
    a41 a42 a43 a44];


% B=[0 0;
%    0 0;
%    b31 b32;
%    b41 b42];
B1=[0; 0 ; b31 ; b41]/1000;

C= [ 1 0 0 0];


rank(ctrb(A,B1)); %cheque controlabilidad del sistema
rank(obsv(A,C)); %chequeo observavilidad del sistema

%% Control LQR y LQI
sys=ss(A,B1,C,0); 
Q=diag([1 8 1 8 0.05]); %para lqi
%%Q2=diag([1 10 1 10]); esta es para lqr
R=1;
[K,S,E]=lqi(sys,Q,R);
%%K1=lqr(sys,Q2,R);
Kp=K(1:4);
ki=K(5);

%discretizacion
Ts=1e-3;
sysD =c2d(sys,Ts,'zoh');
Ad=sysD.A;
Bd=sysD.B;
Cd=sysD.C;
Qd=diag([1 8 1 8 0.05]);
Rd=2;
[Kd,Sd,Ed]= lqi(sysD,Qd,Rd);
Kpd=Kd(1:4);
Kid=Kd(5);


Aclp = (A-B1*Kp);
nsys=ss(Aclp,B1,C,0);

%% observador de orden reducido
Aw= A;
Bw = B1;
Aaa=Aw(1:2,1:2);
Aab=Aw(1:2,3:4);
Aba=Aw(3:4,1:2);
Abb=Aw(3:4,3:4);

Ba=Bw(1:2);
Bb=Bw(3:4);

%para la L del observador reducido
Ao=Abb;
Co=Aab;

Lobs = place(Ao',Co',10*[E(1) E(1)]);
Ke= Lobs';
E3=eig(Ao-Ke*Co);

% Matrices equivalentes para simulink 

A_h = Abb-Ke*Aab;
B_h = A_h*Ke + Aba - Ke*Aaa;
F_h = Bb - Ke*Ba;
C_h = [0 0; 0 0; 1 0; 0 1];
D_h = [1 0; 0 1; Ke];

%% Loop shaping

%primer controlador para el sistema cascada
s = tf('s');
[num,dem]= ss2tf(A,B1,[0 1 0 0],0);
P1 = minreal(zpk(tf(num,dem)),0.01);
Cbeta = ((s+9.45)/(1+(s/42)));%*ksys;
tuvi =((s+9.45)/(1+(s/62)))*P1;
ksys = 1/db2mag(-12.3);
Lsec = P1*Cbeta*ksys;
figure();
margin(P1*Cbeta*ksys) %cheque MP y PH del sistema a lazo abierto
figure()
nyqlog(P1*Cbeta*ksys)

%%  se realiza la interconexion del sistema realimentando solo beta
Csec = Cbeta*ksys;
Cbeta = - minreal(Cbeta*ksys,0.01);
Cbeta = ss(Cbeta);

Cbeta.u = 'beta';
Cbeta.y = 'ut';
Gss= ss(A,B1,[1 0 0 0; 0 1 0 0],[0;0]);
Gss.u = 'u';
Gss.y = 'y';
Sum = sumblk('u = ut + up');
Sys = ss([],[],[],[1 0]);
Sys.u = 'y';
Sys.y = 'alpha';
Sys2 = ss([],[],[],[0 1]);
Sys2.u = 'y';
Sys2.y = 'beta';

Gb = connect(Gss,Cbeta,Sum,Sys,Sys2,'up','alpha');
Gb = zpk(Gb);

Calpha1 = ((s^2 + 52.59*s + 1239)*(s+9.449)/(6.925*(s+62)*(s+8.667)*(s+8.644)))/s;
Calpha2 = -((s+0.2)^2/(1+s/42));

Calpha = minreal(Calpha1*Calpha2,0.01);
Cprim = Calpha;
Kprim=1/db2mag(0);
PK=minreal(Calpha*Gb*Kprim,0.01);
figure()
margin(PK)
figure()
nyqlog(PK)
eig(1/(1-PK))


% Discretizacion del Loop Shaping

Ts = 1e-3;
CsecD = c2d(Csec,Ts,'zoh');
CprimD = c2d(Cprim,Ts,'zoh');




%figure();
%margin(Gb*Calpha)
