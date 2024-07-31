%%% Script for Michaelis Menten reaction.
%%% Note: function definitions are placed at the end of this file.

%%% Parameters
%%%     k1     k-1   k2  
prm = [ 10.0    19    11   ];

%%% Initial conditions (uM)
S0 = 16.0;
E0 = 1.0;
C0 = 0.0;
P0 = 0.0;

%%% set numerical accuracy options
options=odeset('RelTol',1e-6,'AbsTol',1e-9);

%%% solve with stiff-system solver ode15s
%%% [note: X(1)=S, X(2)=E, X(3)=C, X(4)=P]
t0 = 0;
tf = 2;
ics = [S0,E0,C0,P0];
[T,X]=ode15s(@lab1q1,[t0 tf],ics,options,prm);

S=X(:,1);
E=X(:,2);
C=X(:,3);
P=X(:,4);

%%% plots
figure;
plot(T,S,'LineWidth',1.5)
xlabel('Time (sec)');
ylabel('Concentration (uM)');
title('S')

figure;
plot(T,E,'LineWidth',1.5)
xlabel('Time (sec)');
ylabel('Concentration (uM)');
title('E')

%%%%%%%%%%%%%%%%%%%%%%%%
%%% function definitions 
%%%%%%%%%%%%%%%%%%%%%%%%

%%% define ode system for Michaelis Menten reaction
function dXdT = lab1q1(T,X,prm)
 
    %%% parameters
    k1  = prm(1);
    k_1 = prm(2);
    k2  = prm(3);
 
    %%% variables
    s = X(1); 
    e = X(2); 
    c = X(3); 
    p = X(4); 
 
    %%% the differential equations
    dXdT = zeros(4,1);    % a column vector
    dXdT(1) = -k1*e*s+k_1*c;
    dXdT(2) = -k1*e*s+(k_1+k2)*c;
    dXdT(3) = k1*e*s-(k_1+k2)*c;
    dXdT(4) = k2*c;
end