%% Execute Facial Reanimation Control System (FRCS)

%Run after identifying the EMG Response System (ERS) and inverse Stimulus
%Response System (SRS)^1

%Generates an input desired displacement, the corresponding nerual command
%and runs it through the Hammerstein Model followed by the Inverse of the
%LNL Model to get the stimulus

%% User Input Prompts


tStart = tic;

%% Set initial Parameters

Fs = 1000; 
Nfft = 10000;

EMG_response_model = NHK2;
stimulus_response_model = SRS_inverse;

PRBS_movement = false;
% physiological_movement = true;

%PRBS Signal Parameters
PRBS_movement_time = 180;
variable_amplitude = false;
N = PRBS_movement_time/10;
M = 10000;
PRBS_amplitude = 10;            %mm

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_stimulus_max_amplitude = 0.01;
fr = 0.1;                                       %Frequency distribution mean (Hz)
sig = 0.6;                                      %Std of Frequency Distribution (Hz)
W = 0.55;
nf = physiological_movement_time/10;            %Number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (s)
chance_of_zero = false;


%% Generate the Desired Displacement Signal for ERS Simulation

if PRBS_movement == true

    t_total = 0:0.001:PRBS_movement_time;
    time = PRBS_movement_time;

    A = 0;                                  %Intialize amplitude
    if variable_amplitude == true   
        for k = 1:N
            if k == 1
                R = PRBS_amplitude;
            else
                R = rand(1,1)*PRBS_amplitude;   %Randomly generate a number between 0 and PRBS Amplitude
            end
            for j = 1:M
                A = [A R];
            end
        end
    else
        A = PRBS_amplitude;              %Else set as Constant Amplitude
    end

    Range = [0,0.001]; %Specify what the single-channel PRBS value switches between

    %Specify the clock period of the signal as 1 sample. 
    %That is, the signal value can change at each time step. 
    %For PRBS signals, the clock period is specified in Band = [0 B], 
    %where B is the inverse of the required clock period
    %(Must be less than 1)
    Band = [0 0.01];

    %Generate a nonperiodic PRBS of length 100 samples.
    u = idinput(time*1000+1,'prbs',Band,Range);

    %Create an iddata object from the generated signal. 
    %For this example, specify the sample time as 1 second.
    u = iddata([],u,0.001);

    U = (u.InputData)';
    desired_displacement = A.*U;

else

    t_total = 0:0.001:physiological_movement_time;
    time = physiological_movement_time;

    FR = makedist('Normal','mu',fr,'sigma',sig);
    FrequenciesRandom_max = 1.8;
    FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
    freq_distribution = random(FrequenciesRandom,10000,1);

    AR = makedist('Uniform','lower',0,'upper',physiological_stimulus_max_amplitude);  %Full Amplitude Range
    AmplitudesRandom = AR;
    amp_distribution = random(AmplitudesRandom,10000,1);

    desired_displacement= 0;
    Freq_test = [];
    Pulses_per_interval_test = [];

    for j = 1 : nf    
        t  = 0 : 0.001 : t_interval;         % Time Intervals
        
        if j == 1
            Freq = FrequenciesRandom_max;
            A = physiological_stimulus_max_amplitude;
        else
            Freq = random(FrequenciesRandom,1,1);
            A = random(AmplitudesRandom,1,1);
        end

        if chance_of_zero == true
            nums = randi([0 1], 1, 1);
        else
            nums = 0;
        end

        if nums == 0
            g = 1/Freq;
            D = (1:g:t_interval)';     % pulse delay times
            data = (A*pulstran(t,D,@rectpuls,W))';
            stim_frequency = Freq;
            Pulses_per_interval = t_interval/g;
            data(end) = [];
        else
            data = zeros(20000,1);
            stim_frequency = 0;
            Pulses_per_interval = 0;
        end

        desired_displacement = [desired_displacement; data];
        Freq_test = [Freq_test stim_frequency];
        Pulses_per_interval_test = [Pulses_per_interval_test Pulses_per_interval];
    end

    Pulses_per_interval_total = sum(Pulses_per_interval_test);
    Freq_test_average = sum(Freq_test)/length(Freq_test);

    desired_displacement = desired_displacement';

end

%% Generate Neural Input Command Signal for ERS Simulation

%Neural Input Parameters are based on desired displacement amplitude
Amplitude = desired_displacement*100;  %mV
Frequency = desired_displacement*14000;  %Hz

neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

% figure(figNum)
% figNum = figNum+1;
% plot(t_total,neural)
% title('Neural Command','Fontsize',18)
% xlabel('Time (s)','Fontsize',16)
% ylabel('Amplitude (% of MUs)', 'Fontsize', 16)

neural_simulink = [t_total' neural];


%% Generate EMG Signal use EMG Simulation (Simulink Model)

%Run Simulink;
out = sim('EMG_simulation_simulink',time);

%Get Output from Simulink Model
EMG_simulink = out.EMGout;
t_simulink = out.tout;

EMG_simulink = nldat(EMG_simulink);
set(EMG_simulink, 'domainIncr',0.001);

figure(figNum)
figNum = figNum+1;
plot(EMG_simulink)
title('EMG Input (Simulink)','Fontsize',18)
xlabel('Time (s)','Fontsize',16)
ylabel('Amplitude (V)', 'Fontsize', 16)

%% Simulate Response and Plot Results

%Set domain increments of EMG response model
EMG_response_model_IRF = EMG_response_model{1,2};
set(EMG_response_model_IRF, 'domainIncr', 1.0e-3);
EMG_response_model{1,2} = EMG_response_model_IRF;

% set(EMG_response_model,'domainIncr',0.001)

control_system_displacement = nlsim(EMG_response_model,EMG_simulink);
set(control_system_displacement,'domainIncr',0.001);
%control_system_displacement = control_system_displacement.dataSet;

% figure(figNum)
% figNum = figNum+1;
% plot(t_total,control_system_displacement);
% title('Control System Displacement', 'Fontsize', 18)
% xlabel('Time (s)','Fontsize',16)
% ylabel('Displacement (m)', 'Fontsize', 16)

%Set domain increments of Stimulus response model
% stimulus_response_model_IRF = stimulus_response_model{1,1};
% set(stimulus_response_model_IRF, 'domainIncr', 1.0e-3);
% stimulus_response_model{1,1} = stimulus_response_model_IRF;

control_system_output = nlsim(stimulus_response_model,control_system_displacement);
%control_system_stimulus_output = control_system_stimulus_output.dataSet;

% figure(figNum)
% figNum = figNum+1;
% set(control_system_stimulus_output,'domainIncr',0.001)
% plot(t_total,control_system_stimulus_output);
% title('Control System Stimulus Output', 'Fontsize', 18)
% xlabel('Time (s)','Fontsize',16)
% ylabel('Amplitude (V)', 'Fontsize', 16)

%% Plot all in One Figure

figure(figNum)
figNum = figNum+1;

subplot(3,2,1)
plot(t_total,neural)
title('Neural Command','Fontsize',14)
xlabel('Time (s)','Fontsize',14)
ylabel('Amplitude (% of MUs)', 'Fontsize', 14)
grid on

subplot(3,2,2)
EMG_double = EMG_simulink.dataSet;
plot(t_total,EMG_double)
title('EMG Input (Simulink)','Fontsize',14)
xlabel('Time (s)','Fontsize',14)
ylabel('Amplitude (V)', 'Fontsize', 14)
grid on

subplot(3,2,[3 4])
control_system_displacement_double = control_system_displacement.dataSet;
plot(t_total(1,1:end-1999),control_system_displacement_double(1000:end-1000,1));
title('ERS Output (Healthy Side Predicted Displacement)', 'Fontsize', 14)
% xlabel('Time (s)','Fontsize',14)
ylabel('Displacement (m)', 'Fontsize', 14)
grid on

subplot(3,2,[5 6])
control_system_output_double = control_system_output.dataSet;
plot(t_total(1,1:end-1999),control_system_output_double(1000:end-1000,:));
title('Control System Output (Amplitude Modulation)', 'Fontsize', 14)
xlabel('Time (s)','Fontsize',14)
ylabel('Amplitude', 'Fontsize', 14)
grid on

%% Analysis of Control System Displacement and Control System Output (Inverse SRS Input/Output)

%Spectrum of Inverse Input SRS Input and Output

control_system_displacement_zero = control_system_displacement_double(79002:end-1000,1) - mean(control_system_displacement_double(79002:end-1000,1));
control_system_output_zero = control_system_output_double(79002:end-1000,:) - mean(control_system_output_double(79002:end-1000,:));

Nfft = 10000;

[Pxx,~] = pwelch(control_system_displacement_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);
[Pyy,f] = pwelch(control_system_output_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);

figure(figNum)
figNum = figNum+1;
subplot(2,2,1),plot(t_total(1,1:end-1999),control_system_displacement_double(1000:end-1000,1))
title('ERS Output (Healthy Side Predicted Displacement)','Fontsize',20)
xlabel('Time (s)','Fontsize',20)
ylabel('Displacement (m)','Fontsize',20)
grid on

subplot(2,2,2),plot(t_total(1,1:end-1999),control_system_output_double(1000:end-1000,1))
title('Control System Output (Amplitude Modulation)','Fontsize',20)
xlabel('Time (s)','Fontsize',20)
ylabel('Amplitude','Fontsize',20)
grid on

subplot(2,2,3),semilogy(f(1:200,:),Pxx(1:200,:));
title('ERS Output Spectrum','Fontsize',20);
ylabel('PSD','Fontsize',20);
xlabel('Frequency (Hz)','Fontsize',20);
grid on;

subplot(2,2,4),semilogy(f(1:200,:),Pyy(1:200,:));
title('Control System Output Spectrum','Fontsize',20);
ylabel('PSD','Fontsize',20);
xlabel('Frequency (Hz)','Fontsize',20);
grid on;

%Auto-Correlation and Cross Correlation of Inverse SRS Input and Output
figure(figNum)
figNum = figNum+1;
sgtitle('Correlation of Inverse SRS Input/Output','Fontsize',14)

maxlag = length(control_system_displacement_double(79002:end-1000,1))*0.05;
[Rxx,lags] = xcorr(control_system_displacement_zero,maxlag,'coeff');
lags = lags/Fs;
subplot(3,1,1),plot(lags,Rxx)
title('Auto-Correlation of Inverse SRS Input','Fontsize',16)
% xlabel('Time (s)','Fontsize',16)
ylabel('Correlation','Fontsize',16)
grid on

[Ryy,lags] = xcorr(control_system_output_zero,maxlag,'coeff');
lags = lags/Fs;
subplot(3,1,2),plot(lags,Ryy)
title('Auto-Correlation of Inverse SRS Output','Fontsize',16)
% xlabel('Time (s)','Fontsize',16)
ylabel('Correlation','Fontsize',16)
grid on

[Rxy,lags] = xcorr(control_system_output_zero,control_system_displacement_zero,maxlag,'coeff');
lags = lags/Fs;
subplot(3,1,3),plot(lags,Rxy)
title('Cross-Correlation of Inverse SRS Input/Output','Fontsize',16)
xlabel('Time (s)','Fontsize',16)
ylabel('Correlation','Fontsize',16)
grid on

%Frequency Analysis of Inverse SRS Input/Output

% %Compute Input Spectrum, Output Spectrum, and Input-Output Cross Spectrum
% % Choose FFT size and calculate spectrum
% Nfft = 1000;
% [Sxx,~] = pwelch(Rxx,gausswin(Nfft),Nfft/2,Nfft,Fs);
% [Syy,~] = pwelch(Ryy,gausswin(Nfft),Nfft/2,Nfft,Fs);
% [Sxy,f] = pwelch(Rxy,gausswin(Nfft),Nfft/2,Nfft,Fs);
% 
% % Plot frequency spectrum
% figure(figNum)
% figNum = figNum+1;
% subplot(3,1,1),plot(f,Sxx);
% title('Spectrum of Control System Displacement','Fontsize',20);
% ylabel('PSD','Fontsize',20);
% xlabel('Frequency (Hz)','Fontsize',20);
% grid on;
% 
% subplot(3,1,2),plot(f,Syy);
% title('Spectrum of Control System Stimulus','Fontsize',20);
% ylabel('PSD','Fontsize',20);
% xlabel('Frequency (Hz)','Fontsize',20);
% grid on;
% 
% subplot(3,1,3),plot(f,Sxy);
% title('Cross Spectrum of Inverse SRS Input/Output','Fontsize',20);
% ylabel('PSD','Fontsize',20);
% xlabel('Frequency (Hz)','Fontsize',20);
% grid on;
% 
% %Calculate and Plot Gain and Phase Characterisiics(Stochastic Method)
% G = abs(Sxy./Sxx);
% P = angle(Sxy./Sxx);
% Cxy = mscohere(control_system_displacement_double(1000:end-1000,1),control_system_stimulus_output_double(1000:end-1000,:));
% 
% figure(figNum)
% figNum = figNum+1;
% subplot(3,1,1),plot(G);
% title('Gain','Fontsize',20);
% ylabel('Gain (dB)','Fontsize',20);
% xlabel('Frequency (Hz)','Fontsize',20);
% grid on;
% 
% subplot(3,1,2),plot(P);
% title('Phase Shift','Fontsize',20);
% ylabel('Phase (degrees)','Fontsize',20);
% xlabel('Frequency (Hz)','Fontsize',20);
% grid on;
% 
% subplot(3,1,3),plot(Cxy);
% title('Coherence','Fontsize',20);
% ylabel('Coherence','Fontsize',20);
% xlabel('Frequency (Hz)','Fontsize',20);
% grid on;

%Transfer Function
%Nfft = 1000;
txy = tfestimate(control_system_displacement_double(79002:end-1000,1),control_system_output_double(79002:end-1000,:),hamming(Nfft),Nfft/2,Nfft,Fs);
% txy = tfestimate(control_system_displacement_zero,control_system_stimulus_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);
G = abs(txy);
G = 20*log10(G);
P = angle(txy);
P = P*(180/pi);
[Cxy, f] = mscohere(control_system_displacement_double(79002:end-1000,1),control_system_output_double(79002:end-1000,:),hamming(Nfft),Nfft/2,Nfft,Fs);
% [Cxy, f] = mscohere(control_system_displacement_zero,control_system_stimulus_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);

figure(figNum)
figNum = figNum+1;
sgtitle('Frequency Response of Inverse SRS Input/Output','Fontsize',14)

subplot(3,1,1),plot(f(1:31,:),G(1:31,:));
title('Gain','Fontsize',16);
ylabel('Gain (dB)','Fontsize',16);
% xlabel('Frequency (Hz)','Fontsize',16);
grid on;

subplot(3,1,2),plot(f(1:31,:),P(1:31,:));
title('Phase Shift','Fontsize',16);
ylabel('Phase (degrees)','Fontsize',16);
% xlabel('Frequency (Hz)','Fontsize',16);
grid on;

subplot(3,1,3),plot(f(1:31,:),Cxy(1:31,:));
title('Coherence','Fontsize',16);
ylabel('Coherence','Fontsize',16);
xlabel('Frequency (Hz)','Fontsize',16);
grid on;

%% Run the Control System Output (Amplitude Modulation) through the forward model (Paralyzed Side Simulation Model)
% And perform some analysis on the results

predicted_healthy_displacement = control_system_displacement_double(1000:end-1000,1);

%Run Paralyzed Side Simulation Model
%stimulus_simulink = [(t_total(1,1:end-1999))' control_system_output_double(1000:end-1000,:)];
amplitude_modulation = control_system_output_double;
amplitude_modulation_simulink = [(t_total(1,1:end-1999))' amplitude_modulation(1000:end-1000,:)];

% set_param('Paralyzed_Model_Simulink_ControlSystemTest/Output Noise','Cov','set_output_noise_power')
set_param('Paralyzed_Model_Simulink/Output Noise','Cov','set_output_noise_power')
%output_noise_power = [output_noise_power set_output_noise_power];

%Run Simulink;
time = length(t_total(1,1:end-2000))/1000;
% out = sim('Paralyzed_Model_Simulink_ControlSystemTest',time);
out = sim('Paralyzed_Model_Simulink',time);

input_stimulus = out.Paralyzed_Model_Stimulus;
predicted_paralyzed_displacement = out.Paralyzed_Model_Displacement;
t_simulink = out.tout;

%plot the stimulus that is created based on the amplitude modulation
figure(figNum)
figNum = figNum+1;
plot(t_total(1,1:end-1999),input_stimulus)
title('Electrical Stimulus Signal based on Control System Amplitude Modulation Output')
ylabel('Voltage (V)')
xlabel('Time (s)')
grid on

%Compare the predicted healthy displacement and predicted paralyzed
%displacement with the desired displacement
figure(figNum)
figNum = figNum+1;
subplot(4,1,1)
plot(t_total(1,1:end-1999),desired_displacement(:,1000:end-1000))
ax = gca;
ax.FontSize = 10;
title('(a) Desired Displacement', 'Fontsize', 12)
% xlabel('Time (s)', 'Fontsize', 12)
ylabel('Displacement (m)', 'Fontsize', 12)
grid on

subplot(4,1,2)
plot(t_total(1,1:end-1999),predicted_healthy_displacement)
ax = gca;
ax.FontSize = 10;
title('(b) Predicted Healthy Displacement, Pos_H(t)', 'Fontsize', 12)
% xlabel('Time (s)', 'Fontsize', 12)
ylabel('Displacement (m)', 'Fontsize', 12)
grid on

subplot(4,1,3)
plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
ax = gca;
ax.FontSize = 10;
title('(c) Predicted Paralyzed Displacement, Pos_P(t)', 'Fontsize', 12)
% xlabel('Time (s)', 'Fontsize', 12)
ylabel('Displacement (m)', 'Fontsize', 12)
grid on

subplot(4,1,4)
Variance = vaf(predicted_healthy_displacement,predicted_paralyzed_displacement);
hold on
plot(t_total(1,1:end-1999),predicted_healthy_displacement)
plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
hold off
ax = gca;
ax.FontSize = 10;
title(['(d) Superimposed, VAF = ' num2str(round(Variance,1)) '%'], 'Fontsize', 12)
xlabel('Time (s)', 'Fontsize', 12)
ylabel('Displacement (m)', 'Fontsize', 12)
legend('Pos_H(t)', 'Pos_P(t)','Fontsize',10)
grid on

% Compares the desired displacement with the predicted paralyzed
% displacement
figure(figNum)
figNum = figNum+1;
subplot(3,1,1)
plot(t_total(1,1:end-1999),desired_displacement(:,1000:end-1000))
title('(a) Desired Displacement', 'Fontsize', 12)
% xlabel('Time (s)', 'Fontsize', 12)
ylabel('Displacement (m)', 'Fontsize', 12)
grid on

subplot(3,1,2)
plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
title('(b) Predicted Paralyzed Displacement, Pos_P(t)', 'Fontsize', 12)
% xlabel('Time (s)', 'Fontsize', 12)
ylabel('Displacement (m)', 'Fontsize', 12)
grid on

subplot(3,1,3)
Variance2 = vaf((desired_displacement(:,1000:end-1000))',predicted_paralyzed_displacement);
hold on
plot(t_total(1,1:end-1999),desired_displacement(:,1000:end-1000))
plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
title(['(c) Superimposed, VAF = ' num2str(round(Variance2,1)) '%'], 'Fontsize', 12)
xlabel('Time (s)', 'Fontsize', 12)
ylabel('Displacement (m)', 'Fontsize', 12)
legend('Desired Displacement', 'Predicted Paralyzed Displacement, Pos_P(t)')
grid on

%PDF and Spectrum of Residuals (between predicted healthy displacement and
%predicted paralyzed displacement)
residuals = predicted_healthy_displacement - predicted_paralyzed_displacement;

Nfft = 10000;
residuals_zero = residuals - mean(residuals);
%[Prr,f] = pwelch(residuals_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);
[Prr,f] = pwelch(residuals_zero,Nfft,[],Nfft,Fs);

figure(figNum)
figNum = figNum+1;
sgtitle('Residuals between Predicted Healthy Displacement, Pos_H(t), and Predicted Paralyzed Displacement, Pos_P(t)','Fontsize',14)

subplot(2,2,[1 2])
plot(t_total(1,1:end-3998),residuals(1000:end-1000,:))
ax = gca;
ax.FontSize = 14;
title('(a) Residuals','Fontsize',16)
ylabel('Displacement (m)','Fontsize',16)
xlabel('Time (s)','Fontsize',16)
grid on

subplot(2,2,3)
histogram(residuals(1000:end-1000,:))
ax = gca;
ax.FontSize = 14;
title('(b) Residual Distribution','Fontsize',16)
xlabel('Displacement (m)','Fontsize',16)
ylabel('Density','Fontsize',16)
grid on

subplot(2,2,4)
semilogy(f(1:1300,:),Prr(1:1300,:))
ax = gca;
ax.FontSize = 14;
title('(c) Spectrum of Residuals','Fontsize',16)
ylabel('PSD (log)','Fontsize',16);
xlabel('Frequency (Hz)','Fontsize',16);
grid on;

figure(figNum)
figNum = figNum+1;
sgtitle('Residuals between Predicted Healthy Displacement, Pos_H(t), and Predicted Paralyzed Displacement, Pos_P(t)','Fontsize',14)

subplot(2,2,[1 2])
plot(t_total(1,1:end-3998),residuals(1000:end-1000,:))
ax = gca;
ax.FontSize = 14;
title('(a) Residuals','Fontsize',16)
ylabel('Displacement (m)','Fontsize',16)
xlabel('Time (s)','Fontsize',16)
grid on

subplot(2,2,3)
histogram(residuals(1000:end-1000,:))
ax = gca;
ax.FontSize = 14;
title('(b) Residual Distribution','Fontsize',16)
xlabel('Displacement (m)','Fontsize',16)
ylabel('Density','Fontsize',16)
grid on

subplot(2,2,4)
semilogy(f(1:500,:),Prr(1:500,:))
ax = gca;
ax.FontSize = 14;
title('(c) Spectrum of Residuals','Fontsize',16)
ylabel('PSD (log)','Fontsize',16);
xlabel('Frequency (Hz)','Fontsize',16);
grid on;

%Spectrum of Predicted Healthy Displacement and Predicted Paralyzed
%Displacement
predicted_healthy_displacement_zero = predicted_healthy_displacement - mean(predicted_healthy_displacement);
predicted_paralyzed_displacement_zero = predicted_paralyzed_displacement - mean(predicted_paralyzed_displacement);

Nfft = 10000;
%[Pxx,~] = pwelch(predicted_healthy_displacement_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);
%[Pyy,f] = pwelch(predicted_paralyzed_displacement_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);
[Pxx,~] = pwelch(predicted_healthy_displacement_zero,Nfft,[],Nfft,Fs);
[Pyy,f] = pwelch(predicted_paralyzed_displacement_zero,Nfft,[],Nfft,Fs);

figure(figNum)
figNum = figNum+1;
subplot(2,2,1),plot(t_total(1,1:end-1999),predicted_healthy_displacement)
title('Predicted Healthy Displacement, Pos_H(t)','Fontsize',14)
xlabel('Time (s)','Fontsize',16)
ylabel('Displacement (m)','Fontsize',16)
grid on

subplot(2,2,2),plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
title('Predicted Paralyzed Displacement, Pos_P(t)','Fontsize',14)
xlabel('Time (s)','Fontsize',16)
ylabel('Displacement (m)','Fontsize',16)
grid on

subplot(2,2,3),semilogy(f(1:500,:),Pxx(1:500,:));
title('Predicted Healthy Displacement Spectrum','Fontsize',14);
ylabel('PSD (log scale)','Fontsize',16);
xlabel('Frequency (Hz)','Fontsize',16);
grid on;

subplot(2,2,4),semilogy(f(1:500,:),Pyy(1:500,:));
title('Predicted Paralyzed Displacement Spectrum','Fontsize',14);
ylabel('PSD','Fontsize',16);
xlabel('Frequency (Hz)','Fontsize',16);
grid on;

%Auto-Correlation and Cross Correlation of Predicted Healthy Displacement
%and Predicted Paralyzed Displacement
figure(figNum)
figNum = figNum+1;
sgtitle('Correlation of Predicted Healthy Displacement, Pos_H(t), and Predicted Paralyzed Displacement, Pos_P(t)','Fontsize',14)

maxlag = floor(length(predicted_healthy_displacement_zero)*0.05);
[Rxx,lags] = xcorr(predicted_healthy_displacement_zero,maxlag,'coeff');
lags = lags/Fs;
subplot(3,1,1),plot(lags,Rxx)
title('(a) Auto-Correlation of Predicted Healthy Displacement','Fontsize',16)
% xlabel('Time (s)','Fontsize',16)
ylabel('Correlation','Fontsize',16)
grid on

[Ryy,lags] = xcorr(predicted_paralyzed_displacement_zero,maxlag,'coeff');
lags = lags/Fs;
subplot(3,1,2),plot(lags,Ryy)
title('(b) Auto-Correlation of Predicted Paralyzed Displacement','Fontsize',16)
% xlabel('Time (s)','Fontsize',16)
ylabel('Correlation','Fontsize',16)
grid on

[Rxy,lags] = xcorr(predicted_paralyzed_displacement_zero,predicted_healthy_displacement_zero,maxlag,'coeff');
lags = lags/Fs;
subplot(3,1,3),plot(lags,Rxy)
title('(c) Cross-Correlation','Fontsize',16)
xlabel('Time (s)','Fontsize',16)
ylabel('Correlation','Fontsize',16)
grid on

% Frequency Response between Predicted Healthy Side Displacement and
% Predicted Paralyzed Displacement
%Transfer Function
Nfft = 2000;
txy = tfestimate(predicted_healthy_displacement,predicted_paralyzed_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);
%txy = tfestimate(predicted_healthy_displacement,predicted_paralyzed_displacement,[],[],[],Fs);
G = abs(txy);
G = 20*log10(G);
P = angle(txy);
P = P*(180/pi);
[Cxy, f] = mscohere(predicted_healthy_displacement,predicted_paralyzed_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);
%[Cxy, f] = mscohere(predicted_healthy_displacement,predicted_paralyzed_displacement,[],[],[],Fs);

figure(figNum)
figNum = figNum+1;
sgtitle('Frequency Response of Predicted Healthy Displacement, Pos_H(t), and Predicted Paralyzed Displacement, Pos_P(t)','Fontsize',14)

subplot(3,1,1),plot(f(1:600,:),G(1:600,:));
ax = gca;
ax.FontSize = 14;
title('(a) Gain','Fontsize',16);
ylabel('Gain (dB)','Fontsize',16);
% xlabel('Frequency (Hz)','Fontsize',20);
grid on;

subplot(3,1,2),plot(f(1:600,:),P(1:600,:));
ax = gca;
ax.FontSize = 14;
title('(b) Phase Shift','Fontsize',16);
ylabel('Phase (degrees)','Fontsize',16);
% xlabel('Frequency (Hz)','Fontsize',20);
grid on;

subplot(3,1,3),plot(f(1:600,:),Cxy(1:600,:));
ax = gca;
ax.FontSize = 14;
title('(c) Coherence','Fontsize',16);
ylabel('Coherence','Fontsize',16);
xlabel('Frequency (Hz)','Fontsize',16);
grid on;

figure(figNum)
figNum = figNum+1;
plot(f(1:100,:),Cxy(1:100,:));
ax = gca;
ax.FontSize = 16;
title('Coherence between Pos_H(t) and Pos_P(t)','Fontsize',20);
ylabel('Coherence','Fontsize',20);
xlabel('Frequency (Hz)','Fontsize',20);
grid on;

% Frequency Response between Desired Displacement and
% Predicted Paralyzed Displacement
%Transfer Function
Nfft = 1000;
txy = tfestimate(desired_displacement(:,1000:end-1000),predicted_paralyzed_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);
G = abs(txy);
G = 20*log10(G);
P = angle(txy);
P = P*(180/pi);
[Cxy, f] = mscohere(desired_displacement(:,1000:end-1000),predicted_paralyzed_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);

figure(figNum)
figNum = figNum+1;
sgtitle('Frequency Response of Desired Displacement and Predicted Paralyzed Displacement','Fontsize',14)

subplot(3,1,1),plot(f(1:31,:),G(1:31,:));
ax = gca;
ax.FontSize = 14;
title('(a) Gain','Fontsize',16);
ylabel('Gain (dB)','Fontsize',16);
% xlabel('Frequency (Hz)','Fontsize',20);
grid on;

subplot(3,1,2),plot(f(1:31,:),P(1:31,:));
ax = gca;
ax.FontSize = 14;
title('(b) Phase Shift','Fontsize',16);
ylabel('Phase (degrees)','Fontsize',16);
% xlabel('Frequency (Hz)','Fontsize',20);
grid on;

subplot(3,1,3),plot(f(1:31,:),Cxy(1:31,:));
ax = gca;
ax.FontSize = 14;
title('(c) Coherence','Fontsize',16);
ylabel('Coherence','Fontsize',16);
xlabel('Frequency (Hz)','Fontsize',16);
grid on;

% Plotting the residuals as a funtion of desired position
figure(figNum)
figNum = figNum+1;
plot(residuals, desired_displacement(:,1000:end-1000))
title('Desired Displacement vs Residuals','Fontsize',20);
ylabel('Displacement (m)','Fontsize',18);
xlabel('Residual (m)','Fontsize',18);
grid on;

figure(figNum)
figNum = figNum+1;
plot(residuals, predicted_healthy_displacement)
title('Predicted Healthy Displacement vs Residuals','Fontsize',20);
ylabel('Displacement (m)','Fontsize',18);
xlabel('Residual (m)','Fontsize',18);
grid on;

figure(figNum)
figNum = figNum+1;
plot(residuals, predicted_paralyzed_displacement)
title('Predicted Paralyzed Displacement vs Residuals','Fontsize',20);
ylabel('Displacement (m)','Fontsize',18);
xlabel('Residual (m)','Fontsize',18);
grid on;

%% Analysis of Output of ERS and and Output of SRS

% figure(figNum)
% figNum = figNum+1;
% subplot(2,1,1)
% plot(t_total, control_system_displacement_double)
% title('Output of ERS')
% xlabel('Time (s)')
% ylabel('Displacement (m)')
% grid on
% 
% subplot(2,1,2)
% plot(t_total, simulated_outputs(:,1))
% title('Output of SRS')
% xlabel('Time (s)')
% ylabel('Displacement (m)')
% grid on

% figure(figNum)
% figNum = figNum+1;
% hold on
% plot(t_total, control_system_displacement_double)
% plot(t_total, simulated_outputs(:,1))
% hold off
% title('Superimposed ERS and SRS Output')
% xlabel('Time (s)')
% ylabel('Displacement (m)')
% legend('ERS Output','SRS Output')
% grid on
% 
% %Frequency of ERS Output and SRS Output
% 
% control_system_displacement_zero = control_system_displacement_double(1000:end-1000,1) - mean(control_system_displacement_double(1000:end-1000,1));
% simulated_output_zero = simulated_outputs(1000:end-1000,1) - mean(simulated_outputs(1000:end-1000,1));
% 
% [Pxx,f] = pwelch(control_system_displacement_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);
% [Pyy,~] = pwelch(simulated_output_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);
% 
% figure(figNum)
% figNum = figNum+1;
% subplot(2,2,1),plot(t_total(1,1:end-1999),control_system_displacement_double(1000:end-1000,1))
% title('ERS Output','Fontsize',20)
% xlabel('Time (s)','Fontsize',20)
% ylabel('Displacement (m)','Fontsize',20)
% grid on
% 
% subplot(2,2,2),plot(t_total(1,1:end-1999),simulated_outputs(1000:end-1000,1))
% title('SRS Output','Fontsize',20)
% xlabel('Time (s)','Fontsize',20)
% ylabel('Displacement (m)','Fontsize',20)
% grid on
% 
% subplot(2,2,3)
% plot(f(1:20,:),Pxx(1:20,:));
% ax = gca;
% ax.FontSize = 14;
% title('Power Spectrum of ERS Output','Fontsize',20);
% ylabel('PSD','Fontsize',20); 
% xlabel('Frequency (Hz)','Fontsize',20);
% grid on;
% 
% subplot(2,2,4)
% plot(f(1:20,:),Pyy(1:20,:));
% ax = gca;
% ax.FontSize = 14;
% title('Power Spectrum of SRS Output','Fontsize',20);
% ylabel('PSD','Fontsize',20); 
% xlabel('Frequency (Hz)','Fontsize',20);
% grid on;
% 
% %Auto-Correlation and Cross Correlation of Inverse SRS Input and Output
% figure(figNum)
% figNum = figNum+1;
% 
% maxlag = length(control_system_displacement_double(1002:end-1000,1))*0.05;
% [Rxx,lags] = xcorr(control_system_displacement_double(1000:end-1000,1),maxlag,'coeff');
% lags = lags/Fs;
% subplot(3,1,1),plot(lags,Rxx)
% title('Auto-Correlation of ERS Output','Fontsize',20)
% xlabel('Time (s)','Fontsize',20)
% ylabel('Correlation','Fontsize',20)
% grid on
% 
% [Ryy,lags] = xcorr(simulated_outputs(1000:end-1000,1),maxlag,'coeff');
% lags = lags/Fs;
% subplot(3,1,2),plot(lags,Ryy)
% title('Auto-Correlation of SRS Output','Fontsize',20)
% xlabel('Time (s)','Fontsize',20)
% ylabel('Correlation','Fontsize',20)
% grid on
% 
% [Rxy,lags] = xcorr(simulated_outputs(1000:end-1000,1),control_system_displacement_double(1000:end-1000,1),maxlag,'coeff');
% lags = lags/Fs;
% subplot(3,1,3),plot(lags,Rxy)
% title('Cross-Correlation of ERS Output and SRS Output','Fontsize',20)
% xlabel('Time (s)','Fontsize',20)
% ylabel('Correlation','Fontsize',20)
% grid on
% 
% %Transfer Function
% Nfft = 1000;
% txy = tfestimate(control_system_displacement_double(1000:end-1000,1),simulated_outputs(1000:end-1000,1),gausswin(Nfft),Nfft/2,Nfft,Fs);
% G = abs(txy);
% G = 20*log10(G);
% P = angle(txy);
% P = P*(180/pi);
% [Cxy, f] = mscohere(control_system_displacement_double(1000:end-1000,1),simulated_outputs(1000:end-1000,1),gausswin(Nfft),Nfft/2,Nfft,Fs);
% 
% figure(figNum)
% figNum = figNum+1;
% % sgtitle('Frequency Response of ERS Output and SRS Output','Fontsize',20)
% 
% subplot(3,1,1),plot(f,G);
% title('Gain','Fontsize',20);
% ylabel('Gain (dB)','Fontsize',20);
% xlabel('Frequency (Hz)','Fontsize',20);
% grid on;
% 
% subplot(3,1,2),plot(f,P);
% title('Phase Shift','Fontsize',20);
% ylabel('Phase (degrees)','Fontsize',20);
% xlabel('Frequency (Hz)','Fontsize',20);
% grid on;
% 
% subplot(3,1,3),plot(f,Cxy);
% title('Coherence','Fontsize',20);
% ylabel('Coherence','Fontsize',20);
% xlabel('Frequency (Hz)','Fontsize',20);
% grid on;
% 
% %%
% %Spectrum of Control System Displacement and Control System Stimulus\
% 
% [Pxx1,f1] = pwelch(control_system_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);
% 
% figure(figNum)
% figNum = figNum+1;
% 
% subplot(2,1,1)
% plot(t_total,control_system_displacement);
% ax = gca;
% ax.FontSize = 14;
% title('(a) Control System Displacement','Fontsize',18);
% ylabel('PSD (log)','Fontsize',18); 
% xlabel('Frequency (Hz)','Fontsize',18);
% grid on;
% 
% subplot(2,1,2)
% plot(f1,Pxx1);
% ax = gca;
% ax.FontSize = 14;
% title('(b) Power Spectrum of Control System Displacement','Fontsize',18);
% ylabel('PSD (log)','Fontsize',18); 
% xlabel('Frequency (Hz)','Fontsize',18);
% grid on;
% 
% [Pxx2,f2] = pwelch(control_system_stimulus_output,gausswin(Nfft),Nfft/2,Nfft,Fs);
% 
% figure(figNum)
% figNum = figNum+1;
% 
% subplot(2,1,1)
% plot(t_total,control_system_stimulus_output);
% ax = gca;
% ax.FontSize = 14;
% title('(a) Control System Stimulus Output','Fontsize',18);
% ylabel('PSD (log)','Fontsize',18); 
% xlabel('Frequency (Hz)','Fontsize',18);
% grid on;
% 
% subplot(2,1,2)
% plot(f2,Pxx2);
% ax = gca;
% ax.FontSize = 14;
% title('(b) Power Spectrum of Control System Stimulus Output','Fontsize',18);
% ylabel('PSD (log)','Fontsize',18); 
% xlabel('Frequency (Hz)','Fontsize',18);
% grid on;
% 


%%
% control_system_stimulus_output_truncated = max(min(control_system_stimulus_output,10e24),-10e24);
% figure(figNum)
% figNum = figNum+1;
% % set(control_system_stimulus_output,'domainIncr',0.001)
% plot(t_total,control_system_stimulus_output_truncated);
% title('Control System Stimulus Displacement (Truncated)')

%%
tEnd = toc(tStart)/60
