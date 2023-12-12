%% Cleanup
clear variables; clc; close all;

%% Full robot 
%Generate measures
% Define a trajectory
TR.T = 2; %s
TR.N = 3; % nb periods
TR.Q = [pi/2; 0];
TR.C = [0 1 0 0; 1 0 0 0];
TR.S = [0 0 0 0; 0 1 0 2]; 

% Run script
option = 'full_robot' ; % 'full_robot' or  'only_acquisition_hardware'
[q, tau] = myrobot(TR,option);

params.fe_q = 1e3;
params.fe_tau = 2.5e3;

signals.q = q;
signals.tau = tau;


%% Plot signals 

a = figure() ;
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

%% FFT des signaux de positions
N = length(signals.q_aligned(:,1));
frequencies = (-N/2:N/2-1) * params.fe_q / N;
Y1 = fftshift(fft(signals.q_aligned(:,1)));

Y2 = fftshift(fft(signals.q_aligned(:,2)));

figure()
plot(frequencies, abs(Y1));
hold on 
plot(frequencies, abs(Y2));
legend("Axe 1", "Axe 2")

%% Derivation numérique

% Paramètres du filtre passe-bas
fc = 50;  % Fréquence de coupure

% Création du filtre Butterworth d'ordre 6
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

% Plot ou tout autre traitement sur les données de vitesse et d'accélération

temps = 0:1/params.fe_q:((length(position(:,1))-1)/params.fe_q);

% Tracer le signal de position
figure;
subplot(3,1,1);
plot(temps, position(:,1), 'b');
hold on;
plot(temps, position(:,2), 'r');
hold on;
plot(temps, position_filtree(:,1), 'g');
hold on;
plot(temps, position_filtree(:,2), 'y');
hold on;

title('Signaux de position');
legend('Position 1', 'Position 2','Position 1 filtree', 'Position 2 filtree');
xlabel('Temps (s)');
ylabel('Position');
grid on;

% Tracer le signal de vitesse
subplot(3,1,2);
plot(temps(1:end-1), vitesse(:,1), 'b');
hold on;
plot(temps(1:end-1), vitesse(:,2), 'r');
hold on;
plot(temps(1:end-1), vitesse_filtree(:,1), 'g');
hold on;
plot(temps(1:end-1), vitesse_filtree(:,2), 'y');

title('Signaux de vitesse');
legend('Vitesse 1', 'Vitesse 2','Vitesse 1 filtree', 'Vitesse 2 filtree' );
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



%% Regression : 

phi_tot = [];
couple = [];
R = [];


signals.tau_aligned = signals.tau_aligned(3:end,:);
position = position(3:end,:);
vitesse = vitesse(2:end,:);


for i = 10:length(acceleration)-10

    % On appelle fonction  modèle dynamique
    q1 = position_filtree(i,1);
    q2 = position_filtree(i,2);

    dq1 = vitesse_filtree(i,1);
    dq2 = vitesse_filtree(i,2);

    ddq1 = acceleration(i,1);
    ddq2 = acceleration(i,2);
 
    phi = mod_dyn(q1,q2,dq1,dq2,ddq1,ddq2);

    phi_tot = [phi_tot; phi];

    tau_i1 = signals.tau_aligned(i,1);
    tau_i2 = signals.tau_aligned(i,2);
     
    couple = [couple; tau_i1];
    couple = [couple; tau_i2];
end


%% Identification
R = eye(2*(length(acceleration)-19)).*0.009578;
phi = phi_tot;

X = inv(phi'*inv(R)*phi)*phi'*inv(R)*couple







