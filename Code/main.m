clc
clear
close all

% Author    : Abhyudaya Singh
% UIN       : 655691216
% NetID     : singh124@illinois.edu

%% Main

fprintf('\nEnter 1 for solving Toro case1 in Problem 2 \nEnter 2 for solving Toro case2 in Problem 2 \nEnter 3 for solving Problem 3 \nEnter 4 to enter custom conditions')
ch=input('\nEnter choice');

switch ch
    case 1
        L1=0.5;
        L4=0.5;
        BC1="SW";
        BC4="SW";
        q1in=[1,1/287,0];
        q4in=[0.1,0.1/(287*0.125),0];
        Nx=500;
        T=0.25;
        delt=0.0001;
        FUN_Finite_Volume_Solver(L1,BC1,q1in,L4,BC4,q4in,Nx,T,delt,ch);
    case 2
        L1=0.5;
        L4=0.5;
        BC1="SW";
        BC4="SW";
        q1in=[0.4,0.4/287,-2];
        q4in=[0.4,0.4/287,2];
        Nx=500;
        T=0.15;
        delt=0.0001;
        FUN_Finite_Volume_Solver(L1,BC1,q1in,L4,BC4,q4in,Nx,T,delt,ch);
    case 3
        L1=19;
        L4=1;
        BC1="SW";
        BC4="SW";
        q1in=[100*101325,293,0];
        q4in=[0.1*101325,293,0];
        Nx=500;
        T=0.08;
        delt=0.00001;
        FUN_Finite_Volume_Solver(L1,BC1,q1in,L4,BC4,q4in,Nx,T,delt,ch);
    case 4
        L1=input('\nEnter length of the driven section');   %driven section length length
        L4=input('\nEnter length of the driver section');   %driver section length length
                
        q4in=input('\nEnter the driver section initial conditions[p1,T1,u1]');
        q1in=input('\nEnter the driven section initial conditions[p4,T4,u4]');
        
        fprintf('\nEnter boundary conditoin type: \nSW-> solid wall condition \nFF-> far-field condition')
        BC1=input('\ndriven section boundary type','s');    %driven section boundary type
        BC4=input('\ndrive4 section boundary type','s');    %driven section boundary type
        
        Nx=input('\n Enter number of grid points');         %number of cells
        T=input('\nEnter final simulation time');           %total time
        delt=input('\nEnter the time step size');           %time step size
        FUN_Finite_Volume_Solver(L1,BC1,q1in,L4,BC4,q4in,Nx,T,delt,ch);
    otherwise
        fprintf('\n wrong choice')
end
