% Code écrit par Solal BITTOUN, Lilian DELORY et Maxime LELONG - MASTER SAR
% Dans le cadre du TP/DM du cours d'Estimation et Identification.

% Code qui permet de tester le calcul des vitesses er accélérations
% littérales.

%% Cleanup
clear variables; clc; close all;

%% Full robot 
%Generate measures
% Define a trajectory
TR.T = 10; %s
TR.N = 2; % nb periods
TR.Q = [0; 0];
%TR.C = [1.2 0 1 0.4 0.3; 0 0 0 0 0];
%TR.S = [1.1 0 0.2 pi 0.3; 1.2 0.3 0.7 0.4 0.3];

TR.C = [1 0 0 0 0; 0 0 0 0 0];
TR.S = [0 0 0 0 0; 1 0 0 0 0];

% Run script
option = 'full_robot'; % 'full_robot' or  'only_acquisition_hardware'
[q, tau] = myrobot(TR,option);

params.fe_q = 1e3;
params.fe_tau = 2.5e3;

signals.q = q;
signals.tau = tau;


%% Plot signals 
a = figure();
subplot(211); plot(signals.q); xlabel('Time Sample'); ylabel('q'); axis tight 
legend('Axe 1', 'Axe 2');
subplot(212); plot(signals.tau); xlabel('Time Sample'); ylabel('\tau'); axis tight

%% Filtrage 50 Hz

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

signals.q_no_50 = filtfilt(filt_50hz_1k, q);
signals.tau_no_50 = filtfilt(filt_50hz_2_5K, tau);


%% Resample
signals.tau_1k = resample(signals.tau_no_50,params.fe_q,params.fe_tau);
signals.tau_1k = signals.tau_1k(1:end-5,:);

%% t - alignment

signals.tau_aligned = signals.tau_1k(4:end,:);
signals.q_aligned = signals.q_no_50(1:end-3,:);

%% Tests 

% Paramètres
nh = 5;             % Nombre de termes dans la somme
T = 10;             % Période en secondes
fs = 1000;          % Fréquence d'échantillonnage en Hz
duration = T;       % Durée totale du signal (une période)
t = 0:1/fs:duration-1/fs;  % Vecteur temps

% Coefficients
q0i = zeros(2, 1);       % Valeurs initiales nulles
Cij = TR.C;             % Utilisation des coefficients de TR.C
Sij = TR.S;             % Utilisation des coefficients de TR.S

vi = zeros(length(t), 2);
ai = zeros(length(t), 2);

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
xlabel('Temps (s)');
ylabel('Position');
xlim([0 10000]);
ylim([-2 2]);


subplot(3,1,2);
plot(vi);
title('Vitesses');
xlabel('Temps (s)');
ylabel('Vitesse');

subplot(3,1,3);
plot(ai);
title('Accélérations');
xlabel('Temps (s)');
ylabel('Accélération');

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
vitesse = diff(position_filtree) * params.fe_q; % fs pour compenser la différence causée par diff()

% Filtrage des signaux de position 
vitesse_filtree(:,1) = filtfilt(b, a, vitesse(:,1));
vitesse_filtree(:,2) = filtfilt(b, a, vitesse(:,2));

% Calcul de l'accélération à partir du signal de vitesse
acceleration = diff(vitesse_filtree) * params.fe_q; % fs pour compenser la différence causée par diff()

% Acceleration filtrée 

acceleration_filtree(:,1) = filtfilt(b, a, acceleration(:,1));
acceleration_filtree(:,2) = filtfilt(b, a, acceleration(:,2));

% Plot ou tout autre traitement sur les données de vitesse et d'accélération

% Paramètres du filtre passe-bas
fc = 5;  % Fréquence de coupure

% Création du filtre Butterworth d'ordre 8
[b, a] = butter(8, fc/(params.fe_q/2), 'low');

tau_filtree(:,1) =  filtfilt(b, a, signals.tau_aligned(:,1));
tau_filtree(:,2) =  filtfilt(b, a, signals.tau_aligned(:,2));

temps = 0:1/params.fe_q:((length(position(:,1))-1)/params.fe_q);

% Tracer le signal de position
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

% Tracer le signal de vitesse
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

% Tracer le signal d'accélération
subplot(3,1,3);
plot(temps(1:end-2), acceleration(:,1), 'b');
hold on;
plot(temps(1:end-2), acceleration(:,2), 'r');
title('Signaux d''accélération');
legend('Accélération 1', 'Accélération 2');
xlabel('Temps (s)');
ylabel('Accélération');
grid on;














