function []=FUN_Finite_Volume_Solver(L1,BC1,q1in,L4,BC4,q4in,Nx,T,delt,ch)

% Author    : Abhyudaya Singh
% UIN       : 655691216
% NetID     : singh124@illinois.edu

% This code solves the 1D Euler's equation using Finite Volume method using Godunov's method for flux reconstruction

%% initialization

gama=1.4;
R=287;

L=L1+L4;
delx=L/(Nx);

N_dia=round(L4/delx);       % node at diaphram

q0=zeros(10000,Nx);
q1=zeros(10000,Nx);
q2=zeros(10000,Nx);

FL=zeros(3,Nx);
FR=zeros(3,Nx);

%% initalizing flow field
rho1in=q1in(1)/(R*q1in(2));
rho4in=q4in(1)/(R*q4in(2));
q0(1,1:N_dia) = ones(1,N_dia)*rho1in;
q1(1,1:N_dia) = q0(1,1:N_dia)*q1in(3);
q2(1,1:N_dia) = ones(1,N_dia)*(q1in(1)/(gama-1) + 0.5*rho1in*q1in(3)^2);

q0(1,N_dia+1:end) = ones(1,Nx-N_dia)*rho4in;
q1(1,N_dia+1:end) = q0(1,N_dia+1:end)*q4in(3);
q2(1,N_dia+1:end) = ones(1,Nx-N_dia)*(q4in(1)/(gama-1) + 0.5*rho4in*q4in(3)^2);


%% solution
fprintf('\nCalculating the solution!...')
fprintf('\n\titer \t Flow time')

t=zeros(10000,1);
ct=0;
while t(ct+1)<=T
    ct=ct+1;
    for n=1:Nx
        
        if n==1
            if BC1 == "SW"
                if q1(n)>0
                    qL=FUN_Godunov_flux_recon([q0(ct,n), 0, q2(ct,n)],[q0(ct,n),q1(ct,n),q2(ct,n)],gama,0);
                else
                    qL=FUN_Godunov_flux_recon([2*q0(ct,n)-q0(ct,n+1), 0, 2*q2(ct,n)-q2(ct,n+1)],[q0(ct,n),q1(ct,n),q2(ct,n)],gama,0);
                end
            elseif BC1 == "FF"
                if q1(n)>0
                    qL=FUN_Godunov_flux_recon([q0(ct,n), q1(ct,n), q2(ct,n)],[q0(ct,n),q1(ct,n),q2(ct,n)],gama,0);
                else
                    qL=FUN_Godunov_flux_recon([2*q0(ct,n)-q0(ct,n+1), 2*q1(ct,n)-q1(ct,n+1), 2*q2(ct,n)-q2(ct,n+1)],[q0(ct,n),q1(ct,n),q2(ct,n)],gama,0);
                end
            end
            qR=FUN_Godunov_flux_recon([q0(ct,n),q1(ct,n),q2(ct,n)],[q0(ct,n+1),q1(ct,n+1),q2(ct,n+1)],gama,0);
            PL=(gama-1)*(qL(3)-(0.5*qL(2)^2 /qL(1)));
            FL(1:3,n)=[qL(2), qL(2)^2 /qL(1) + PL, (qL(3)+PL)*qL(2)/qL(1)]';
            PR=(gama-1)*(qR(3)-(0.5*qR(2)^2 /qR(1)));
            FR(1:3,n)=[qR(2), qR(2)^2 /qR(1) + PR, (qR(3)+PR)*qR(2)/qR(1)]';
        
        elseif n==Nx

            if BC4 == "SW"
                if q1(n)>0
                    qR=FUN_Godunov_flux_recon([q0(ct,n),q1(ct,n),q2(ct,n)],[2*q0(ct,n)-q0(ct,n-1), 0, 2*q2(ct,n)-q2(ct,n-1)],gama,0);
                else
                    qR=FUN_Godunov_flux_recon([q0(ct,n),q1(ct,n),q2(ct,n)],[q0(ct,n),0,q2(ct,n)],gama,0);
                end
            elseif BC4 == "FF"
                if q1(n)>0
                    qR=FUN_Godunov_flux_recon([q0(ct,n),q1(ct,n),q2(ct,n)],[2*q0(ct,n)-q0(ct,n-1), 2*q1(ct,n)-q1(ct,n-1), 2*q2(ct,n)-q2(ct,n-1)],gama,0);
                else
                    qR=FUN_Godunov_flux_recon([q0(ct,n),q1(ct,n),q2(ct,n)],[q0(ct,n),q1(ct,n),q2(ct,n)],gama,0);
                end
            end
            qL=FUN_Godunov_flux_recon([q0(ct,n-1),q1(ct,n-1),q2(ct,n-1)],[q0(ct,n),q1(ct,n),q2(ct,n)],gama,0);
            PL=(gama-1)*(qL(3)-(0.5*qL(2)^2 /qL(1)));
            FL(1:3,n)=[qL(2), qL(2)^2 /qL(1) + PL, (qL(3)+PL)*qL(2)/qL(1)]';
            PR=(gama-1)*(qR(3)-(0.5*qR(2)^2 /qR(1)));
            FR(1:3,n)=[qR(2), qR(2)^2 /qR(1) + PR, (qR(3)+PR)*qR(2)/qR(1)]';

        else

            qL=FUN_Godunov_flux_recon([q0(ct,n-1),q1(ct,n-1),q2(ct,n-1)],[q0(ct,n),q1(ct,n),q2(ct,n)],gama,0);
            qR=FUN_Godunov_flux_recon([q0(ct,n),q1(ct,n),q2(ct,n)],[q0(ct,n+1),q1(ct,n+1),q2(ct,n+1)],gama,0);
            PL=(gama-1)*(qL(3)-(0.5*qL(2)^2 /qL(1)));
            FL(1:3,n)=[qL(2), qL(2)^2 /qL(1) + PL, (qL(3)+PL)*qL(2)/qL(1)]';
            PR=(gama-1)*(qR(3)-(0.5*qR(2)^2 /qR(1)));
            FR(1:3,n)=[qR(2), qR(2)^2 /qR(1) + PR, (qR(3)+PR)*qR(2)/qR(1)]';
        end
    end
    
        q0(ct+1,:) = q0(ct,:) - (delt/delx)*(FR(1,:)-FL(1,:));
        q1(ct+1,:) = q1(ct,:) - (delt/delx)*(FR(2,:)-FL(2,:));
        q2(ct+1,:) = q2(ct,:) - (delt/delx)*(FR(3,:)-FL(3,:));
        
        if ct>1
            fprintf(repmat('\b',1,lineLength))
        end
        lineLength = fprintf('\n\t%d \t %f',ct,t(ct));
        
        t(ct+1)=t(ct)+delt;

end

fprintf('\nCalculation complete!')

%% plotting
if ch == 1 || ch == 2
    if ch==1
        toro= table2array(readtable('Test1_Toro1.txt'));
    else
        toro= table2array(readtable('Test2_Toro1.txt'));
    end
    xana=toro(:,1);
    rhoana=toro(:,4);
    uana=toro(:,5);
    Pana=toro(:,6);
    
    figure(1)
    plot(-L4+delx:delx:L1,q0(ct,:))
    grid on
    hold on
    plot(xana,rhoana)
    xlabel('x')
    ylabel('rho')
    legend('Numerical solution','Toro')
    
    figure(2)
    plot(-L4+delx:delx:L1,q1(ct,:)./q0(ct,:))
    grid on
    hold on
    plot(xana,uana)
    xlabel('x')
    ylabel('u')
    legend('Numerical solution','Toro')
    
    figure(3)
    plot(-L4+delx:delx:L1,(gama-1)*(q2(ct,:)-(0.5*(q1(ct,:).^2)./q0(ct,:))))
    grid on
    hold on
    plot(xana,Pana)
    xlabel('x')
    ylabel('P')
    legend('Numerical solution','Toro')
end

if ch == 3 || ch == 4
    figure(4)
    contourf(-L4+delx:delx:L1,t(1:10:ct),(gama-1)*(q2(1:10:ct,:)-(0.5*(q1(1:10:ct,:).^2)./q0(1:10:ct,:))),100,'edgecolor','none')
    a=colorbar;
    ylabel(a,'P (Pa)','Rotation',90);
    xlabel('x(m)')
    ylabel('t(sec)')
    
    figure(5)
    contourf(-L4+delx:delx:L1,t(1:10:ct),q0(1:10:ct,:),100,'edgecolor','none')
    a=colorbar;
    ylabel(a,'\rho (kg/m^{3})','Rotation',90);
    xlabel('x(m)')
    ylabel('t(sec)')
    
    figure(6)
    contourf(-L4+delx:delx:L1,t(1:10:ct),q1(1:10:ct,:)./q0(1:10:ct,:),100,'edgecolor','none')
    a=colorbar;
    ylabel(a,'U (m/s)','Rotation',90);
    xlabel('x(m)')
    ylabel('t(sec)')
    
    figure(7)
    contourf(-L4+delx:delx:L1,t(1:10:ct),((gama-1)*(q2(1:10:ct,:)-(0.5*(q1(1:10:ct,:).^2)./q0(1:10:ct,:))))./(R*q0(1:10:ct,:)),100,'edgecolor','none')
    a=colorbar;
    ylabel(a,'T (K)','Rotation',90);
    xlabel('x(m)')
    ylabel('t(sec)')
end
end