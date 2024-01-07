% Code écrit par Solal BITTOUN, Lilian DELORY et Maxime LELONG - MASTER SAR
% Dans le cadre du TP/DM du cours d'Estimation et Identification.

% Code qui permet de paramétrer fmincon pour notre porblème.

function [solution,fval,exitflag,output]=minimize_cost_function(TR, x0)
    % Paramètres
    nh = TR.nh;             % Nombre de termes dans la somme
    T = TR.T;             % Période en secondes
    fs = TR.fe_q;          % Fréquence d'échantillonnage en Hz
    duration = TR.T * TR.N;       % Durée totale du signal
    t = 0:1/fs:duration-1/fs;  % Vecteur temps

    % summary of the optimization problem-general form
    % fmincon attempts to solve problems of the form:
    %     min J(x)  subject to:  Aineq*x  <= Bineq, Aeq*x  = Beq (linear constraints)
    %      x                    C(x) <= 0, Ceq(x) = 0   (nonlinear constraints)
    %                               lb <= x <= ub        (bounds)
     
     
    Aineq=[];
    Bineq=[];
    Aeq=[];
    Beq=[];
    lb=[];
    ub=[];
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % keyboard
    % algorithm can be 'interior-point', 'SQP','active set', and 'trust region reflective'
    % set the options
    options = optimoptions(@fmincon,'Algorithm','interior-point','SpecifyObjectiveGradient',false,'Display','iter','MaxFunctionEvaluations',10000,'UseParallel',true,'OptimalityTolerance',1.0000e-10, 'StepTolerance',1.0000e-10)
    
    [solution,fval,exitflag,output] = fmincon(@(x)cost_function(x),x0,Aineq,Bineq,Aeq,Beq,lb,ub,[],options);
     
    function res = cost_function(x)
        % notice that this function can access the arguments of the function
        % minimize_cost_function()

        q0i = x(:, 1);
        Cij = x(:, 2:nh+1);
        Sij = x(:, nh+2:2*nh+1);
        
        qi = zeros(length(t), 2);
        vi = zeros(length(t), 2);
        ai = zeros(length(t), 2);
        % Calculs des vitesses et accélérations.
        phi_tot = [];
        for i = 1:length(t)
            for axe=1:2
                qi(i,axe) = q0i(axe);
                for j=1:nh
                    qi(i,axe) = qi(i,axe) + Cij(axe, j)*cos(2*pi*j*t(i)/T) + Sij(axe, j)*sin(2*pi*t(i)/T);
                    vi(i,axe) = vi(i,axe) + (((2*pi*j)/T)*(-Cij(axe,j)*sin((2*pi*j*t(i))/T) + Sij(axe,j)*cos((2*pi*j*t(i))/T)));
                    ai(i,axe) = ai(i,axe) + (((2*pi*j)/T)^(2)*(-Cij(axe,j)*cos((2*pi*j*t(i))/T) - Sij(axe,j)*sin((2*pi*j*t(i))/T)));
                end
            end
        end
        nb_samples = 8000;
        selected_time_idxs = randsample(length(t), nb_samples);
        for i = selected_time_idxs'
        
            
            q1 = qi(i,1);
            q2 = qi(i,2);
        
            dq1 = vi(i,1);
            dq2 = vi(i,2);
        
            ddq1 = ai(i,1);
            ddq2 = ai(i,2);
        
            % Appel à la fonction mod_dyn qui permet de calculer la matrice phi du
            % modèle linéaire.
            phi = mod_dyn(q1,q2,dq1,dq2,ddq1,ddq2);
            
            % On ajoute phi.
            phi_tot = [phi_tot; phi];
        
        end

        res = cond(phi_tot);
        if res == inf
            res = 1e50;
        end
         
    end
 
end