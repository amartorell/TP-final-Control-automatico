% datos del sistema, joint 1
m1a = 0.0011644908507951884 + 0.0011644908507951884 + 0.0011644908507951884 + 0.00058767552757679601 + 0.0010321252629775849 + 0.012717581311371812;
m1b = 0.00052672018519800228;
m1 = m1a + m1b ; 

aw1= 5.71e+4; %sacado del grafico con link 1 bien armado
j1= 1/(aw1);   %obtenemos el momento de inercia con respecto al punto de pivote
L1=0.1035;

%cosntantes generales
g=9.81;
b1=1e-6;
b2=1e-6;

% datos del sistema, joint 2
m2a= 0.00013195457400179582 + 0.00052672018519800228 + 0.00043271308549453917 + 0.00086673289514495463;
m2b= 0.00043271308549453917 + 0.0010695037303870972 + 0.0012730475598330083 + 0.00043271308549453917 + 0.00058767552757679558;
m2= m2a+m2b;

% experimento d, sacamos la distancia al centro de masa del joint 2
 l2=(2e-3)/(m2*9.81*sin(0.6592));
 
 
 
% experimento c, sacamos el periodo del pendulo obtenemos inercia de 2
T=725.926e-3;
j2=(T/(2*pi))^2*m2*g*l2;



%linearizacion del modelo

J1= j1;
J2 = j2 + m2*l2^2;
J0= J1 + m2*L1;

m2l2 = m2*l2;       %producto utilizado varias veces

denom = (J0*J2-(m2l2*L1)^2);

a32= (g*m2l2^2*L1)/denom;
a33= (-b1*J2)/denom;
a34= (-b2*m2l2*L1)/denom;

a42= (g*m2l2*J0)/denom;
a43= (-b1*m2l2*L1)/denom;
a44=(-b2*J0)/denom;

b31=J2/denom;








