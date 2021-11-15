clc
clear
%% exact
k=94;
r=.1;
L=5;  %length
h2=10;
teta0=30;
m2=2*h2*L/(k*r)
m=sqrt(m2)
x=[.01 1 2 3 4 5];
c=teta0/(L^-.5*besseli(1, 2*m*L^(1/2)))
tetexact=c.*x.^-.5.*besseli(1, 2*m*x.^(1/2))
%% approximate
A2=((m2*L^2*teta0)/2)/(2*L^3+(m2*L^4)/2-(m2*L^4)/4)
A0=teta0-A2*L^2
tetapproximate=A0+A2*x.^2
%% ARD
ARD=abs(tetexact-tetapproximate)./tetexact