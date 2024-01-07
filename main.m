
% Code écrit par Solal BITTOUN, Lilian DELORY et Maxime LELONG - MASTER SAR
% Dans le cadre du TP/DM du cours d'Estimation et Identification.
% Code "main" permettant de réaliser l'identification paramétrique en utilisant un estimateur MLE ou LS

%% Cleanup
clear variables; clc; close all;

%% Full robot 

% Modèle simplifié avec theta1 = m2. (1 sinon non = 0)
simplifie = 0;

TR.T = 10;     %Durée de la période
TR.N = 2;      % Nombre de période
TR.Q = [0 ; 0]; % Position de départ.

% Trajectoire 1
TR.C = [1.2 0 1 0.4 0.3; 0 0 1 0 0];
TR.S = [1.1 0 0.2 pi 0.3; 1.2 0.3 0.7 0.4 0.3];

% Trajectoire 2
%TR.C = [pi 0 0 0 0; 0 0 0 1 0];
%TR.S = [0 0 1 0 0; pi 0 0 0 0];



% % Trajectoire obtenue avec fmincon.
  %TR.C= [3.2904  , 0.7161  , 0.2157  , 1.0262  , 0.0879 ; -1.0219  , 0.2683 ,  0.3577 ,   0.4425 ,  1.1029];
% 
 %TR.S = [ 0.7021, 0.5119, -0.3817 ,  -0.9182 , -0.6793;  1.6586 ,  1.2675 , 0.4810  , 0.3084  , 0.0739];
% 
 %TR.Q = [-0.5287; 0.7020];


 

% Run script
option = 'full_robot'; % 'full_robot' or  'only_acquisition_hardware'
[q, tau] = myrobot(TR,option);

% Fréquence d'échantillonnage des positions angulaires
params.fe_q = 1e3;
% Fréquence d'échantillonnage des couples.
params.fe_tau = 2.5e3;


signals.q = q;
signals.tau = tau;


%% Plot signaux de positions et de couples.

a = figure();
subplot(211); plot(signals.q); xlabel('Time Sample'); ylabel('q'); axis tight 
legend('Axe 1', 'Axe 2');
xlabel("Echantillons")
ylabel("Position (en rad)")
title("Signaux de positions")
subplot(212); plot(signals.tau); xlabel('Time Sample'); ylabel('\tau'); axis tight
xlabel("Echantillons")
ylabel("Position (en rad)")
title("Signaux de couples")
legend('Axe 1', 'Axe 2');
sgtitle("Signaux de Base");

%% Filtrage 50 Hz sur les deux types de signaux

% Filtre pour les signaux de positions
filt_50hz_1k = designfilt('bandstopiir', ...       % Response type
       'PassbandFrequency1',48, ...    % Frequency constraints
       'StopbandFrequency1',49.5, ...
       'StopbandFrequency2',50.5, ...
       'PassbandFrequency2',52, ...
       'PassbandRipple1',0.1, ...         % Magnitude constraints
       'StopbandAttenuation',15, ...
       'PassbandRipple2',0.1, ...
       'DesignMethod','ellip', ...      % Design method
       'SampleRate',params.fe_q) ;              % Sample rate

% Filtre pour les signaux de couples.
filt_50hz_2_5K = designfilt('bandstopiir', ...       % Response type
       'PassbandFrequency1',48, ...    % Frequency constraints
       'StopbandFrequency1',49.5, ...
       'StopbandFrequency2',50.5, ...
       'PassbandFrequency2',52, ...
       'PassbandRipple1',0.1, ...         % Magnitude constraints
       'StopbandAttenuation',15, ...
       'PassbandRipple2',0.1, ...
       'DesignMethod','ellip', ...      % Design method
       'SampleRate',params.fe_tau) ;              % Sample rate

% On applique le filtre
signals.q_no_50 = filtfilt(filt_50hz_1k, q);
signals.tau_no_50 = filtfilt(filt_50hz_2_5K, tau);


%% Resample des signaux de couples.

signals.tau_1k = resample(signals.tau_no_50,params.fe_q,params.fe_tau);
signals.tau_1k = signals.tau_1k(1:end-5,:);

%% t - alignment - On décale de 4 échantillons les signaux de couples.
signals.tau_aligned = signals.tau_1k(4:end,:);
signals.q_aligned = signals.q_no_50(1:end-3,:);

%% FFT des signaux de positions
N = length(signals.q_aligned(:,1));
frequencies = (-N/2:N/2-1) * params.fe_q / N;
Y1 = fftshift(fft(signals.q_aligned(:,1)));
Y2 = fftshift(fft(signals.q_aligned(:,2)));

figure()
plot(frequencies, abs(Y1));
hold on 
plot(frequencies, abs(Y2));
xlabel("Fréquences (en Hz)")
ylabel("Amplitude (en dB)")
legend("Axe 1", "Axe 2")
title("FFT des signaux de postions");

%% Derivation numérique

% Paramètres du filtre passe-bas
fc = 5;  % Fréquence de coupure

% Création du filtre Butterworth d'ordre 8
[b, a] = butter(8, fc/(params.fe_q/2), 'low');

position = signals.q_aligned;

% Filtrage des signaux de position
position_filtree(:,1) = filtfilt(b, a, position(:,1));
position_filtree(:,2) = filtfilt(b, a, position(:,2));

% Calcul de la vitesse à partir du signal de position filtré
vitesse = diff(position_filtree) * params.fe_q; 

% Filtrage des signaux de vitesses 
vitesse_filtree(:,1) = filtfilt(b, a, vitesse(:,1));
vitesse_filtree(:,2) = filtfilt(b, a, vitesse(:,2));

% Calcul de l'accélération à partir du signal de vitesse
acceleration = diff(vitesse_filtree) * params.fe_q; 

% Acceleration filtrée 
acceleration_filtree(:,1) = filtfilt(b, a, acceleration(:,1));
acceleration_filtree(:,2) = filtfilt(b, a, acceleration(:,2));

% Filtrage des mesures de couples si nécessaire.

% Paramètres du filtre passe-bas
fc = 5;  % Fréquence de coupure
% Création du filtre Butterworth d'ordre 8
[b, a] = butter(8, fc/(params.fe_q/2), 'low');

% Application du filtre.
tau_filtree(:,1) =  filtfilt(b, a, signals.tau_aligned(:,1));
tau_filtree(:,2) =  filtfilt(b, a, signals.tau_aligned(:,2));

% Echelle de temps.
temps = 0:1/params.fe_q:((length(position(:,1))-1)/params.fe_q);

% Tracage des signaux de positions
figure;
subplot(3,1,1);
plot(temps, position(:,1), 'b');
hold on;
plot(temps, position(:,2), 'r');
hold on;
plot(temps, position_filtree(:,1), 'g',Linewidth=1.5);
hold on;
plot(temps, position_filtree(:,2), 'y',Linewidth=1.5);
hold on;

title('Signaux de position');
legend('Position 1', 'Position 2','Position 1 filtrée', 'Position 2 filtrée');
xlabel('Temps (s)');
ylabel('Position');
grid on;

% Tracage des signaux de vitesses
subplot(3,1,2);
plot(temps(1:end-1), vitesse(:,1), 'b');
hold on;
plot(temps(1:end-1), vitesse(:,2), 'r');
hold on;
plot(temps(1:end-1), vitesse_filtree(:,1), 'g',Linewidth=1.5);
hold on;
plot(temps(1:end-1), vitesse_filtree(:,2), 'y',Linewidth=1.5);


title('Signaux de vitesse');
legend('Vitesse 1', 'Vitesse 2','Vitesse 1 filtrée', 'Vitesse 2 filtrée' );
xlabel('Temps (s)');
ylabel('Vitesse');
grid on;

% Tracage des signaux d'accélérations.
subplot(3,1,3);
plot(temps(1:end-2), acceleration(:,1), 'b');
hold on;
plot(temps(1:end-2), acceleration(:,2), 'r');
title('Signaux d''accélération');
legend('Accélération 1', 'Accélération 2');
xlabel('Temps (s)');
ylabel('Accélération');
grid on;
sgtitle("Comparaison signaux filtrés et non filtrés")

%% Tracage des signaux de couples.
figure()

subplot(2,1,1);
plot(signals.tau_aligned(:,1), 'b');
hold on;
plot(tau_filtree(:,1), 'r');
xlabel('Echantillons')
ylabel('Couple (en Nm)')
title("Comparaison des couples mesurés et mesurés et filtrés de l'axe 1")
grid on;

subplot(2,1,2);
plot(signals.tau_aligned(:,2), 'b');
hold on;
plot(tau_filtree(:,2), 'r');
xlabel('Echantillons')
ylabel('Couple (en Nm)')
title("Comparaison des couples mesurés et mesurés et filtrés de l'axe 2")
grid on;

sgtitle("Comparaison signaux de couples filtrés et non filtrés");



%% Vitesses et accélerations issues de la formule de la trajectoire : 

% Paramètres
nh = 5;             % Nombre de termes dans la somme
T = TR.T;             % Période en secondes
fs = params.fe_q;          % Fréquence d'échantillonnage en Hz
duration = TR.T * TR.N;       % Durée totale du signal
t = 0:1/fs:duration-1/fs;  % Vecteur temps

% Coefficients
q0i = TR.Q;             % Valeurs initiales nulles
Cij = TR.C;             % Utilisation des coefficients de TR.C
Sij = TR.S;             % Utilisation des coefficients de TR.S

vi = zeros(length(t), 2);
ai = zeros(length(t), 2);

% Calculs des vitesses et accélérations.
for i = 1:length(t)
    for axe=1:2
        for j=1:nh
            vi(i,axe) = vi(i,axe) + (((2*pi*j)/T)*(-Cij(axe,j)*sin((2*pi*j*t(i))/T) + Sij(axe,j)*cos((2*pi*j*t(i))/T)));
            ai(i,axe) = ai(i,axe) + (((2*pi*j)/T)^(2)*(-Cij(axe,j)*cos((2*pi*j*t(i))/T) - Sij(axe,j)*sin((2*pi*j*t(i))/T)));
        end
    end
end


% Affichage des signaux
figure;

subplot(3,1,1);
plot(signals.q);
title('Signaux de Position');
xlabel('Echantillon')
ylabel('Position');
xlim([0 10000]);
ylim([-2 2]);


subplot(2,1,1);
plot(t,vi);
title('Vitesse');
xlabel('Temps (en s)')
ylabel('Vitesse');
legend('Axe1 1','Axe 2')

subplot(2,1,2);
plot(t,ai);
title('Accélération');
xlabel("Temps (en s)")
ylabel('Accélération');
legend('Axe 1', 'Axe 2')

sgtitle("Signaux de vitesses et d'accélérations issus de l'expression de la trajectoire");

%% Comparaison des signaux de vitesses et d'acceleration
figure()

subplot(2,1,1);
plot(vi,Linewidth=1.5);
hold on;
plot(vitesse_filtree)
title('Vitesse');
xlabel("Echantillon");
ylabel('Vitesse');
legend("Littérale 1","Littérale 2","Derivation numérique 1","Dérivation numérique 2");

subplot(2,1,2);
plot(ai,Linewidth=1.5);
hold on;
plot(acceleration_filtree)
title('Accélération');
xlabel("Echantillon");
ylabel('Accélération');
legend("Littérale 1","Littérale 2","Derivation numérique 1","Dérivation numérique 2");

sgtitle("Comparison des signaux de vitesses et d'accélération");

%% Regression Linéaire 

% Matrice qui va contenir les sous matrices du modèle pour chaque
% échantillon.
phi_tot = [];

% Matrice qui va contenir les sous matrices de couples pour chaque
% échantillon.
couple = [];

% Matrice de covariance du bruit de mesure.
R = [];


% Si on souhaite utiliser les vitesses et accélérations de la dérivation
% numérique : décommenter.
% signals.tau_aligned = signals.tau_aligned(3:end,:);
% position = position(3:end,:);
% vitesse = vitesse(2:end,:);


for i = 100:length(acceleration)-100

    % On appelle fonction  modèle dynamique
    q1 = position_filtree(i,1);
    q2 = position_filtree(i,2);

    % Si on souhaite utiliser les vitesses et accélérations de la dérivation
    % numérique : décommenter.
    % dq1 = vitesse_filtree(i,1);
    % dq2 = vitesse_filtree(i,2);
    % 
    % ddq1 = acceleration_filtree(i,1);
    % ddq2 = acceleration_filtree(i,2);

    dq1 = vi(i,1);
    dq2 = vi(i,2);

    ddq1 = ai(i,1);
    ddq2 = ai(i,2);

    % Appel à la fonction mod_dyn qui permet de calculer la matrice phi du
    % modèle linéaire.
    phi = mod_dyn(q1,q2,dq1,dq2,ddq1,ddq2);
    
    % On ajoute phi.
    phi_tot = [phi_tot; phi];
    
    tau_i1 = tau_filtree(i,1);
    tau_i2 = tau_filtree(i,2);

    if (simplifie == 1)
         tau_i1 = tau_filtree(i,1) - (2*(0.2^2) + 0.04)*ddq1;
    end
    
    % On ajoute les deux nouvelles valeurs de couples.
    couple = [couple; tau_i1];
    couple = [couple; tau_i2];

    % Permet d'afficher l'évolution de la boucle dans la commande.
    disp(i);
end



%% Identification

% On forme la matrice de covariance du bruit.
% Si MLE -> 0.36 sinon 1
variance = 1;

R = eye(2*(length(acceleration)-199)).*variance; 
phi = phi_tot;

% On fait le calcul d'identification :  X = inv(phi'*inv(R)*phi)*phi'*inv(R)*couple;

x = R\phi;
x= (phi')*x;

fprintf("Optimalité de la trajectoire : ")
disp(-log(det(inv(x))));
fprintf("Cond de Phi : ")
disp(cond(phi))

x= x\(phi');
x = x/(R);

disp("Paramètres estimés") % x contient l'estimation des paramètres.
x = x*couple



%% Computed Torque issus du modèle identifié 

% Paramètres utilisés dans le rapport. 
% Utilisation de l'estimation courante
%       -> param_ide = x

param_ide =   [
    0.1341 ;
    0.1072 ;
    0.3237 ;
   -0.0150 ;
    0.1216 ;
   -0.0197 ;
   -0.3223 ;
    0.0457 ;
    ];

% Param avec fmincon.
param_ide = [
    0.0637;
    0.0028;
    0.2217;
   -0.0128;
    0.0861;
   -0.0028;
   -0.4176;
   -0.0203];

compute_torque = [];
for i = 100:length(acceleration)-100

    % On appelle fonction modèle dynamique
    q1 = position_filtree(i,1);
    q2 = position_filtree(i,2);

    % Choix d'utiliser les vitesses et accélerations issues de la
    % dérivation numérique ou expression littérale.

    % dq1 = vitesse_filtree(i,1);
    % dq2 = vitesse_filtree(i,2);
    % 
    % ddq1 = acceleration_filtree(i,1);
    % ddq2 = acceleration_filtree(i,2);

    dq1 = vi(i,1);
    dq2 = vi(i,2);

    ddq1 = ai(i,1);
    ddq2 = ai(i,2);

    couple_cal = computed_torque(q1,q2,dq1,dq2,ddq1,ddq2,param_ide);

    if simplifie == 1 
        couple_cal(1) = couple_cal(1) + (2*(0.2^2) + 0.04)*ddq1;
    end
   
    compute_torque = [compute_torque; couple_cal'];

end


%% Computed Torque issus du modèle avec les données du sujet.

param_approx =   [
    0.37 ;
    0.04 ;
    0.2 ;
    0.06 ;
    0.0001;
    0.0001 ;
    0.01 ;
    0.01 ;
    ];

compute_torque_approx = [];

for i = 100:length(acceleration)-100

    % On appelle fonction  modèle dynamique
    q1 = position_filtree(i,1);
    q2 = position_filtree(i,2);

    % dq1 = vitesse_filtree(i,1);
    % dq2 = vitesse_filtree(i,2);
    % 
    % ddq1 = acceleration_filtree(i,1);
    % ddq2 = acceleration_filtree(i,2);

    dq1 = vi(i,1);
    dq2 = vi(i,2);

    ddq1 = ai(i,1);
    ddq2 = ai(i,2);

    couple_cal = computed_torque(q1,q2,dq1,dq2,ddq1,ddq2,param_approx);

    if simplifie == 1 
        couple_cal(1) = couple_cal(1) + (2*(0.2^2) + 0.04)*ddq1;
    end
    compute_torque_approx = [compute_torque_approx; couple_cal'];

end



%% Tracer le signal de couples - Comparaison estimé - réel

figure();

subplot(2,1,1);
plot(signals.tau_aligned(:,1), 'g');
hold on;
plot(compute_torque(:,1), 'b');
hold on;
plot(tau_filtree(:,1), 'r');
hold on;
plot(compute_torque_approx(:,1),'y');
legend('Couple mesuré','Couple prédit','Couple mesuré et filtré','Couple avec données approximatives'); 
title('Comparaison du couple prédit et réel de l''axe 1');
xlabel('Echantillon');
ylabel("Couple (en Nm)");


grid on;

subplot(2,1,2);
plot(signals.tau_aligned(:,2), 'g');
hold on;
plot(compute_torque(:,2), 'b');
hold on;
plot(tau_filtree(:,2), 'r');
hold on;
plot(compute_torque_approx(:,2),'y');
legend('Couple mesuré','Couple prédit','Couple mesuré et filtré','Couple avec données approximatives');
title('Comparaison du couple prédit et réel de l''axe 2');
xlabel('Echantillon');
ylabel("Couple (en Nm)");

grid on;

