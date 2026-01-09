function [pstar,ustar,rhoLstar,cLstar,rhoRstar,cRstar] = FUN_Root_Bisection(qL,qR,gama)

% Author    : Abhyudaya Singh
% UIN       : 655691216
% NetID     : singh124@illinois.edu

% this function finds the root using the Bisection method
% the qk passed in this function should have qk=[rhok, uk, Pk];

itmax=100;
a = 0.9 * (qL(3) + qR(3))/2;
b = 1.1 * (qL(3) + qR(3))/2;
it = 0;

% Make sure f(a) is negative
while qR(2) - qL(2) + FUN_fk(a,qL,gama) + FUN_fk(a,qR,gama) > 0 && it < itmax
    a = a/2;
    it = it + 1;
end

it = 0;

% Make sure f(b) is positive
while qR(2) - qL(2) + FUN_fk(b,qL,gama) + FUN_fk(b,qR,gama) < 0 && it < itmax
    b = b*2;
    it = it + 1;
end

% Bisection algorithm

it = 0;
e = 1e-6;
pstar = (a + b)/2;
fpstar = qR(2) - qL(2) + FUN_fk(pstar,qL,gama) + FUN_fk(pstar,qR,gama);

while abs(fpstar) > e && (b - a)/2 > e && it < itmax
    pstar = (a + b)/2;
    fpstar = qR(2) - qL(2) + FUN_fk(pstar,qL,gama) + FUN_fk(pstar,qR,gama);
    if fpstar > 0
        b = pstar;
    else
        a = pstar;
    end
    it = it + 1;
end

ustar = (qR(2) + qL(2) + FUN_fk(pstar,qR,gama) - FUN_fk(pstar,qL,gama))/2;
rhoLstar = qL(1)*(pstar/qL(3))^(1/gama);
cLstar = sqrt(gama*pstar/rhoLstar);
rhoRstar = qR(1)*(pstar/qR(3))^(1/gama);
cRstar = sqrt(gama*pstar/rhoRstar);
end