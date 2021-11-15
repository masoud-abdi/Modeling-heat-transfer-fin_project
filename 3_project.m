clc
clear
%% exact
syms m x c
k=94;
r=.1;
q=150;
x=5;  %length
h2=10;
m2=2*h2*x/(k*r)
m=sqrt(m2)
xx=[.01 1 2 3 4 5];
z=2*m*x^.5;
teta=c*x^-.5*besseli(1,z);
f=diff(teta);
p=-(besseli(1, 2*m*x^(1/2)))/(2*x^(3/2)) - ((besseli(1, 2*m*x^(1/2))/(2*x) - (m*besseli(0, 2*m*x^(1/2)))/x^(1/2)))/x^(1/2);
c=(q/k)/p
zz=2*m*xx.^.5;
tetexact=c.*xx.^-.5.*besseli(1,zz)
%% approximate
A2=q/2/k/x
A0=(4*A2*x^3-(m2*A2*x^4/2))/(m2*x^2)
tetapproximate=A0+A2.*xx.^2
ARD=abs(tetexact-tetapproximate)./tetexact
