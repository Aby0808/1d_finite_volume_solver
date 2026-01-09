function [Q] = FUN_Godunov_flux_recon(q1,q2,gama,xt)

% Author    : Abhyudaya Singh
% UIN       : 655691216
% NetID     : singh124@illinois.edu

% This function constructs the fluxes on the cell faces by using the Riemann solver

%converting vectors from [rho, rhou, rhoE] to [rho, u, P]
qL=[q1(1), q1(2)/q1(1), (gama-1)*(q1(3)-(0.5*q1(2)^2 /q1(1)))];
qR=[q2(1), q2(2)/q2(1), (gama-1)*(q2(3)-(0.5*q2(2)^2 /q2(1)))];

cL=sqrt(gama*qL(3)/qL(1));  % sound speed on left
cR=sqrt(gama*qR(3)/qR(1));  % sound speed on right

%% solving for P*
[pstar,ustar,rhoLstar,cLstar,rhoRstar,cRstar] = FUN_Root_Bisection(qL,qR,gama);


%% Calculating the fluxes

sL = qL(2)-cL*sqrt(((gama+1)*pstar)/(2*gama*qL(3)) + ((gama-1)/(2*gama)));
sR = qR(2)+cR*sqrt(((gama+1)*pstar)/(2*gama*qR(3)) + ((gama-1)/(2*gama)));

if xt<ustar
    if qL(3)>=pstar
        if xt<=qL(2)-cL
            q=qL;
        elseif xt<=ustar-cLstar
            rhofL = qL(1)*((2/(gama+1)) + ((gama-1)/(gama+1))*((qL(2)-xt)/cL))^(2/(gama-1));
            ufL = ((gama-1)/(gama+1)*(qL(2) + 2*((xt+cL)/(gama-1))));
            PfL = qL(3)*((2/(gama+1)) + ((gama-1)/(gama+1))*((qL(2)-xt)/cL))^(2*gama/(gama-1));
            q=[rhofL,ufL,PfL];
        else%if xt>ustar-cLstar && xt<=ustar
            q=[rhoLstar,ustar,pstar];
        end
    else
        if xt<=sL
            q=qL;
        else%if xt>sL && xt<=ustar
            q=[rhoLstar,ustar,pstar];
        end
    end

else
    if qR(3)>=pstar
        if xt<=ustar+cRstar
            q=[rhoRstar,ustar,pstar];
        elseif xt<=qR(2)+cR
            rhofR = qR(1)*((2/(gama+1)) + ((gama-1)/(gama+1))*((qR(2)-xt)/cR))^(2/(gama-1));
            ufR = ((gama-1)/(gama+1)*(qR(2) + 2*((xt-cR)/(gama-1))));
            PfR = qR(3)*((2/(gama+1)) + ((gama-1)/(gama+1))*((qR(2)-xt)/cR))^(2*gama/(gama-1));
            q=[rhofR,ufR,PfR];
        else
            q=qR;
        end
    else
        if xt<=sR
            q=[rhoRstar,ustar,pstar];
        else%if xt>sR
            q=qR;
        end
    end
end
Q=[q(1), q(1)*q(2), q(3)/(gama-1) + 0.5*q(1)*q(2)^2];       %[rho,rhou,rhoE]
end