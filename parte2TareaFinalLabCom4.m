clc
clear
pkg load signal;
pkg load optim;
[x, fs] = audioread('x2_U017.wav');

lx= length(x);

[y, fs] = audioread('y2_U017.wav');

%Graficando vectores
n=0:lx-1;
figure(1);
subplot(2,1,1);
stem(n,x);
title('X[n]');
subplot(2,1,2);
stem(n,y);
title('Y[n]');

%[b,a]=calculocoef2(x,y)
clc
lx=5000;
q=1:lx;
n=1:length(x);

%eq simultaneas
for i=2001:lx-1;
  i;
  c1(i)=x(i-2000);
  c2(i)=-y(i-750);
  b1(i)=y(i)-x(i);
endfor

C1=transpose(c1);
C2=transpose(c2);
B=transpose(b1);
D=[C1 C2];
lb= [-Inf; -Inf];
ub=[Inf; Inf];
R=lsqlin(D,B, [], [], [], [], lb, ub);

fprintf('EL valor para c/u de los coeficientes es: ');
fprintf('\n A= ')
disp(R(1,1));
fprintf('   B= ')
disp(R(2,1));

%comparacion evaluación de coeficientes
fprintf('comparacion de valores para las posiciones entre [2001] y [20016]');
for j=2001:2016
  m(j-2000)=x(j)+R(1,1)*(j-2000)+R(2,1)*y(j-750);
  l(j-2000)=y(j);
  s(j-2000)=l(j-2000)-m(j-2000);
endfor

M=transpose(m);
L=transpose(l);
S=transpose(s);

comparacion=[M L S]

b(1)=1;
b(2)=R(1,1);
a(1)=1;
a(2)=R(2,1);

fprintf('\n coeficientesde fracciones parciales: \n');
[r,p,k]=residuez(b,a);
r=transpose(r);
p=transpose(p);

b=1:96000;
b=b*0;
a=b;
b(1)=1;
b(2000)=R(1,1);
a(1)=1;
a(750)=R(2,1);
for i=1:lx
  z=i;
  HH(i)=sum(r./(1-p*z^(-1)))+sum(k*[1 ; z^(-1)]);
endfor

HH;
%analisis espectro

f=1/10;
N1=30;
N2=120;

%analis espectro X
X1=abs(fft(x,N1));
X2=abs(fft(x,N2));

%rango normalizado transformadas
F1x=[(0:N1-1)/N1];
F2x=[(0:N2-1)/N2];

%analisis expectro Y
Y1=abs(fft(y,N1));
Y2=abs(fft(y,N2));

%rango normaliado transformadas
F1y=[(0:N1-1)/N1];
F2y=[(0:N2-1)/N2];

%grafica funciones
% X Y
figure(2);
subplot(2,1,1);
stem(n,x);
title('X[n]');

subplot(2,1,2);
stem(n,y)
title('Y[n]');

%graficas espectrales

figure(3)
subplot(2,2,1);
stem(F1x,X1, '.');
title('Espectro par X[n] con N=30');

subplot(2,2,2);
stem(F2x, X2, ' . ')
title('Espectro para X[n] con N=120');

subplot (2,2,3);
stem(F1y,Y1, '. ')
title('Espectro par Y[n] con N=30');

subplot(2,2,4);
stem(F2y,X2, ' . ')
title('Espectro par Y[n] con N=120');

%Partes real e imaginaria

figure(4);
subplot(2,1,1);
stem(q, real(HH));
title('Real(H[z])');

subplot(2,1,2);
stem(q, imag(HH));
title('Imaginaria H[z]');

%Respuesta amplitud y frecuencia
figure(5)
y=filter(b,a,q);
plot(q,y);
title('Respuesta al impulso y escalón unitario');

%polos y ceros
figure(6);
zplane(b,a);
title('Diagrama de polos y ceros');
