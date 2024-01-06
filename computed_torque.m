% Code écrit par Solal BITTOUN, Lilian DELORY et Maxime LELONG - MASTER SAR
% Dans le cadre du TP/DM du cours d'Estimation et Identification.

% Fonction qui permet d'obtenir le couple prédit avec les paramètres
% identifiés.


function torque = computed_torque(q1,q2,dq1,dq2,ddq1,ddq2,param)
    
     % Il y a neuf paramètres.
    phi = zeros(2,9);

    %  -----------  TERMES ASSOCiES A LA MATRICE MASSE  --------------
    % Seulement associés aux paramètres theta 1 jusqu'à theta 4

    % Associés à theta 1
    phi(1,1) = phi(1,1) + ddq1;
    phi(2,1) = phi(2,1) + 0;
    
    % Associés à theta 2
    phi(1,2) = phi(1,2) + ((sin(q2))^2)*ddq1;
    phi(2,2) = phi(2,2) + 0;

    % Associés à theta 3
    phi(1,3) = phi(1,3) + cos(q2)*ddq2*0.5;
    phi(2,3) = phi(2,3) + cos(q2)*ddq1*0.5;

    % Associés à theta 4
    phi(1,4) = phi(1,4) + 0;
    phi(2,4) = phi(2,4) + ddq2;



    %  -----------  TERMES ASSOCiES A LA MATRICE C  --------------
    % Seulement associés aux paramètres theta 2 et theta 3

    % Associés à theta 2
    phi(1,2) = phi(1,2) + sin(2*q2)*dq2*dq1;
    phi(2,2) = phi(2,2) -0.5*(dq1)^2 * sin(2*q2);

    % Associés à theta 3
    phi(1,3) = phi(1,3) - (dq2)^2 * sin(q2)*0.5 ;
    phi(2,3) = phi(2,3) + 0;


    %  -----------  TERMES ASSOCiES A LA GRAVITE G  --------------
    % Seulement associés aux paramètres theta 5 seulement
    
    % Modif ici
    % Associés à theta 5
    %phi(1,5) = phi(1,5) + ddq1*0.25;
    phi(2,3) = phi(2,3) - sin(q2)*9.81;

    %  -----------  TERMES ASSOCiES AUX FROTTEMENTS VISQUEUX   --------------
    % Seulement associés aux paramètres theta 6 et 7 seulement
    
    phi(1,6) = phi(1,6) + dq1;
    phi(2,7) = phi(2,7) + dq2;

    
    %  -----------  TERMES ASSOCiES AUX FORCES DE COULOUMB   --------------
    % Seulement associés aux paramètres theta 8 et 9 seulement

    % Associés à theta 8 et 9 - les frottements s'appliquent sur chacun des
    % segements

    r = 0.01;
    phi(1,8) = phi(1,8) + sign(dq1); %tanh(r*dq1);
    phi(2,9) = phi(2,9) + sign(dq2); %tanh(r*dq2);

    phi = phi(:,[1:4, 6:9]);

    torque = phi * param;

end