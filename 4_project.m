clc
clear
%% exact
syms m x h
k=94;
r=.1;
x=.5;  %length
h2=10;
m2=2*h2*x/(k*r)
m=sqrt(m2)
q=200;
hinf=150;
Tinf1=20;
Tinf2=100;
xx=[.01 .2 .3 .35 .4 .5];
z=2*m*x^.5;
teta=x^-.5*besseli(1,z);
f=diff(teta);
p=-(besseli(1, 2*m*x^(1/2)))/(2*x^(3/2)) - ((besseli(1, 2*m*x^(1/2))/(2*x) - (m*besseli(0, 2*m*x^(1/2)))/x^(1/2)))/x^(1/2);
c=(q/hinf-Tinf1+Tinf2)/(k*p/hinf+teta)
zz=2*m*xx.^.5;
tetexact=c.*xx.^-.5.*besseli(1,zz)
%% approximate
A2=m2*x^2/2*(q/hinf+Tinf2-Tinf1)/(2*x^3+m2*k*x^3/hinf+m2*x^4/2-m2*x^4/4)
A0=q/hinf-2*k*A2*x/hinf-A2*x^2-Tinf1+Tinf2
tetapproximate=A0+A2.*xx.^2
ARD=abs(tetexact-tetapproximate)./tetexact
