% Code écrit par Solal BITTOUN, Lilian DELORY et Maxime LELONG - MASTER SAR
% Dans le cadre du TP/DM du cours d'Estimation et Identification.

% Code qui permet d'appliquer fmincon a notre problème.


clear var; close all; clc;

TR.T = 10;     %Durée de la période
TR.N = 2;      % Nombre de période
TR.nh = 5;
TR.fe_q = 1e3;

% solve the problem
x0 = rand(2, 2*TR.nh+1);


Cij = [pi 0 0 0 0; 0 0 0 1 0];
Sij = [0 0 1 0 0; pi 0 0 0 0];

q0 = zeros(2,1);
x0 = [Cij, Sij, q0];
disp(x0);

[solution,fval,exitflag,output] = minimize_cost_function(TR,x0)
 
q0i = solution(:, 1);
Cij = solution(:, 2:TR.nh+1);
Sij = solution(:, TR.nh+2:2*TR.nh+1);