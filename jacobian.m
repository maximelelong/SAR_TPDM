
% Code écrit par Solal BITTOUN, Lilian DELORY et Maxime LELONG - MASTER SAR
% Dans le cadre du TP/DM du cours d'Estimation et Identification.
% Code "jacobian" permettant de réaliser le calcul de la jacobienne pour la
% méthode de descente de gradient.

function [torque,J] = jacobian(q1,q2,dq1,dq2,ddq1,ddq2,bg)

    torque1 = 0;
    torque2 = 0;
    J = zeros(2,10);

    m1 = 2;
    m2 = bg(2);

    l1 = 0.2;
    l2 = bg(4);

    I1 = 0.04;
    I2 = bg(6);

    fv1 = bg(7);
    fv2 = bg(8);

    fc1 = bg(9);
    fc2 = bg(10);

    % Calcul des couples
   
    torque1 = (m1*(l1^(2)) + I1 + m2*((0.5)^2) + m2*(l2^2)*(sin(q2)^2))*ddq1 + (0.5*l2*m2*cos(q2))*ddq2;

    torque2 = (0.5*l2*m2*cos(q2))*ddq1 + (m2*(l2^2)+I2)*ddq2;

    % Centrifuges et Coriolis 

    torque1 = torque1 + m2*l2*l2*sin(2*q2)*dq1*dq2 - 0.5*l2*m2*sin(q2)*(dq2^2);
    torque2 = torque2 - 0.5*m2*l2*sin(2*q2)*(dq1^2);

    % Frottements et gravité
    torque1 = torque1 + fv1*dq1 + fc1*sign(dq1);
    torque2 = torque2 + fv2*dq2 + fc2*sign(dq2) - l2*m2*9.81*sin(q2);

    torque = [torque1;torque2];

    % disp("Predicted Torque");
    % disp(torque);

    % Expression de la Jacobienne
    J(1,1) = (l1^2)*ddq1;
    J(1,2) = (0.5^2)*ddq1 + l2*l2*((sin(q2))^2)*ddq1 + (0.5*l2*cos(q2))*ddq2 + (sin(2*q2)*l2*l2*dq1*dq2) - 0.5*l2*(dq2^2)*sin(q2);
    J(1,3) = (2*l1*m1)*ddq1;
    J(1,4) = (2*l2*m2*((sin(q2))^2))*ddq1 + 0.5*m2*cos(q2)*ddq2 + 2*l2*m2*sin(2*q2)*dq1*dq2 -0.5*m2*(dq2^2)*sin(q2);
    J(1,5) = ddq1;
    J(1,6) = 0;
    J(1,7) = dq1;
    J(1,8) = 0;
    J(1,9) = sign(dq1);
    J(1,10) = 0;

    J(2,1) = 0;
    J(2,2) = 0.5*l2*cos(q2)*ddq1 + l2*l2*ddq2 -0.5*l2*l2*(dq1^2)*sin(2*q2) - l2*9.81*sin(q2);
    J(2,3) = 0;
    J(2,4) = 0.5*m2*cos(q2)*ddq1 + 2*l2*m2*ddq2 -0.5*2*l2*m2*(dq1^2)*sin(2*q2) - m2*9.81*sin(q2);
    J(2,5) = 0;
    J(2,6) = ddq2;
    J(2,7) = 0;
    J(2,8) = dq2;
    J(2,9) = 0;
    J(2,10) = sign(dq2);

end