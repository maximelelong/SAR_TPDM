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
q1 = signals.q_no_50(:,1);
mu = mean(q1);
V = var(q1);
fprintf("bruit -> N(%f, %f)\n", mu, V);
fprintf("Variance tau -> %f \n", var(signals.q_no_50(:,2)));