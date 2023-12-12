%% Cleanup
clear variables; clc; close all;

%% Generate measures

% Define a trajectory
TR.T = 1 ; %s
TR.N = 3; % nb periods
% TR.Q = [pi/2; 0];
% TR.C = [0 1 0 0; 1 0 0 0];
% TR.S = [0 0 0 0; 0 1 0 2]; 

TR.Q = [0; 0];
TR.C = [0 0 0 0; 0 0 0 0];
TR.S = [0 0 0 0; 0 0 0 0]; 
    
% Run script
option = 'only_acquisition_hardware' ; % 'full_robot' or  'only_acquisition_hardware'
[q, tau] = myrobot(TR,option);

params.fe_q = 1e3;
params.fe_tau = 2.5e3;

signals.q = q;
signals.tau_1k = tau;

%% Plot results
figure()
subplot(211); plot(q); xlabel('Time Sample'); ylabel('q');
legend('q_1', 'q_2'); title("Mesures avec interface position") ;
subplot(212); plot(tau); xlabel('Time Sample'); ylabel('\tau');
legend('\tau_1', '\tau_2') ; title("Mesures avec interface couple ");

%% Filter 50Hz
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

%% Plot ffts
figure;
[freq_q1, fft_q1] = compute_fft(signals.q(:,2), params.fe_q);
[~, fft_q1_filt] = compute_fft(signals.q_no_50(:,2), params.fe_q);
plot_fft_orig = plot(freq_q1, fft_q1, LineWidth=2); hold on;
plot_fft_orig.Color(4) = 0.30;
plot(freq_q1, fft_q1_filt,LineWidth=2);
title("Comparaison signal avec et sans bruits de secteur");
legend("q_1" ,"q_1 filt");

%% White noise analysis
% With no input to the system => only measurement noise (assumed gaussian)
q1 = q(:,1);
mu = mean(q1);
V = var(q1);
fprintf("bruit -> N(%.3f, %.3f)\n", mu, V);


%% Full robot 
%Generate measures

%close all 
clear all
clc
close all

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

%% Derivation numérique

% Paramètres du filtre passe-bas
fc = 30;  % Fréquence de coupure

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
acceleration = diff(vitesse) * params.fe_q; % fs pour compenser la différence causée par diff()

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
title('Signaux de vitesse');
legend('Vitesse 1', 'Vitesse 2');
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
    q1 = position(i,1);
    q2 = position(i,2);

    dq1 = vitesse(i,1);
    dq2 = vitesse(i,2);

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

R = eye(2*(length(acceleration)-19)).*0.018;


phi = phi_tot;

X = inv(phi'*inv(R)*phi)*phi'*inv(R)*couple







