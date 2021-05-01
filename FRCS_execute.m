%% Execute Facial Reanimation Control System (FRCS)

%INSTRUCTIONS:
%Run after identifying the EMG Response System (ERS) and inverse Stimulus
%Response System (SRS)^1 (First run FRCS_identify_models followed by
%FRCS_identify_inverse_SRS)

%Generates a desired displacement and the corresponding input nerual command
%and runs the neural command through the ERS Simulation to get an EMG Signal. 
%This EMG signal is then used as input to the FRCS (ERS followed by inverse SRS).
%The output of the FRCS is a stimulus amplitude modulation signal required
%to duplicate the displacement on the paralyzed side. The amplitude
%modulation signal is passed into the SRS Simulation to get a paralyzed
%displacement that is compared to the predicted healthy displacement.

%When running the script, you need to provide the following input:
% 1. Type of FRCS Input? PRBS/Physiological
%       This determines which type of input is used for the FRCS. Default
%       is Physiological as this is closer to a real world signal.

%% User Input Prompts

prompt1 = 'Type of FRCS Input? PRBS/Physiological(Phys) [Phys]: ';
str1 = input(prompt1,'s');
if ~strcmp(str1,'PRBS') & ~strcmp(str1,'Phys') & ~isempty(str1)
    disp('Invalid Input')
    return
elseif isempty(str1)
    str1 = 'Phys';
end

if strcmp(str1,'PRBS')
    PRBS_signal = true;
elseif strcmp(str1,'Phys')
    PRBS_signal = false;
end 
FRCS_input_type = str1;

%Number of Trials

%Number of Models

tStart = tic;

%% Set initial Parameters

Fs = 1000; 
Nfft = 10000;
figNum = 200;

%ERS and SRS Model
if strcmp(SRS_inverse_input_type,'PRBS')
    EMG_response_system = ERS_all(1);
elseif strcmp(SRS_inverse_input_type,'Phys')
    EMG_response_system = ERS_all(2);
end

%PRBS Signal Parameters
PRBS_movement_time = 180;
variable_amplitude = false;
N = PRBS_movement_time/10;
M = 10000;
PRBS_amplitude = 10;                            %PRBS Amplitude (mm)

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

%This is used to get an EMG Signal for Input to the FRCS

if PRBS_signal == true

    t_total = 0:0.001:PRBS_movement_time;
    time = PRBS_movement_time;

    A = 0;                                      %Intialize amplitude
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

    %Generate a nonperiodic PRBS of length time
    u = idinput(time*1000+1,'prbs',Band,Range);

    %Create an iddata object from the generated signal. 
    %For this example, specify the sample time as 0.001 second.
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
neural_simulink = [t_total' neural];

%% Generate EMG Signal use EMG Simulation (Simulink Model)

%Run Simulink;
out = sim('ERS_simulation',time);

%Get Output from Simulink Model (EMG)
EMG_simulink = out.ERS_Simulation_EMG;
t_simulink = out.tout;

EMG_simulink = nldat(EMG_simulink);
set(EMG_simulink, 'domainIncr',0.001);

%% Simulate Response and Plot Results

%Set domain increments of EMG response model
EMG_response_system_IRF = EMG_response_system{1,2};
set(EMG_response_system_IRF, 'domainIncr', 1.0e-3);
EMG_response_system{1,2} = EMG_response_system_IRF;

%Simulate Response of ERS to EMG Input
control_system_displacement = nlsim(EMG_response_system,EMG_simulink);
set(control_system_displacement,'domainIncr',0.001);

%Simulate Response of Inverse SRS to Healthy Displacement (Output of ERS)
control_system_output = nlsim(SRS_inverse,control_system_displacement);
set(control_system_output,'domainIncr',0.001);

%% Plot Inputs and Output of Facial Reanimation Control System

figure(figNum)
figNum = figNum+1;
sgtitle('Inputs/Outputs of FRCS','Fontsize',14)

subplot(3,2,1)
plot(t_total,neural)
title('Neural Command Input for ERS Simulation','Fontsize',14)
xlabel('Time (s)','Fontsize',14)
ylabel('Amplitude (% of MUs)', 'Fontsize', 14)
grid on

subplot(3,2,2)
EMG_double = EMG_simulink.dataSet;
plot(t_total,EMG_double)
title('EMG Input from ERS Simulation','Fontsize',14)
xlabel('Time (s)','Fontsize',14)
ylabel('Amplitude (V)', 'Fontsize', 14)
grid on

subplot(3,2,[3 4])
control_system_displacement_double = control_system_displacement.dataSet;
plot(t_total(1,1:end-1999),control_system_displacement_double(1000:end-1000,1));
title('ERS Output (Healthy Side Predicted Displacement)', 'Fontsize', 14)
ylabel('Displacement (m)', 'Fontsize', 14)
grid on

subplot(3,2,[5 6])
control_system_output_double = control_system_output.dataSet;
plot(t_total(1,1:end-1999),control_system_output_double(1000:end-1000,:));
title('FRCS Output (Amplitude Modulation)', 'Fontsize', 14)
xlabel('Time (s)','Fontsize',14)
ylabel('Amplitude', 'Fontsize', 14)
grid on

%% Run the Control System Output (Amplitude Modulation) through SRS Simulation

predicted_healthy_displacement = control_system_displacement_double(1000:end-1000,1);

%Execute SRS Simulation with amplitude modulation input
amplitude_modulation = control_system_output_double;
amplitude_modulation_simulink = [(t_total(1,1:end-1999))' amplitude_modulation(1000:end-1000,:)];

%Set Output Noise as zero
set_output_noise_power = 0;
set_param('SRS_simulation/Output Noise','Cov','set_output_noise_power')

%Run Simulink;
time = length(t_total(1,1:end-2000))/1000;
out = sim('SRS_simulation',time);

%Get Output signals from SRS Simulation
electrical_stimulus = out.SRS_Simulation_Stimulus;
predicted_paralyzed_displacement = out.SRS_Simulation_Displacement;
t_simulink = out.tout;

%plot the stimulus that is created based on the amplitude modulation
figure(figNum)
figNum = figNum+1;
plot(t_total(1,1:end-1999),electrical_stimulus)
title('Electrical Stimulus Signal based on FRCS Amplitude Modulation Output')
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
ylabel('Displacement (m)', 'Fontsize', 12)
grid on

subplot(4,1,2)
plot(t_total(1,1:end-1999),predicted_healthy_displacement)
ax = gca;
ax.FontSize = 10;
title('(b) Predicted Healthy Displacement, Pos_H(t)', 'Fontsize', 12)
ylabel('Displacement (m)', 'Fontsize', 12)
grid on

subplot(4,1,3)
plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
ax = gca;
ax.FontSize = 10;
title('(c) Predicted Paralyzed Displacement, Pos_P(t)', 'Fontsize', 12)
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

% % Compares the desired displacement with the predicted paralyzed
% % displacement
% figure(figNum)
% figNum = figNum+1;
% subplot(3,1,1)
% plot(t_total(1,1:end-1999),desired_displacement(:,1000:end-1000))
% title('(a) Desired Displacement', 'Fontsize', 12)
% ylabel('Displacement (m)', 'Fontsize', 12)
% grid on
% 
% subplot(3,1,2)
% plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
% title('(b) Predicted Paralyzed Displacement, Pos_P(t)', 'Fontsize', 12)
% ylabel('Displacement (m)', 'Fontsize', 12)
% grid on
% 
% subplot(3,1,3)
% Variance2 = vaf((desired_displacement(:,1000:end-1000))',predicted_paralyzed_displacement);
% hold on
% plot(t_total(1,1:end-1999),desired_displacement(:,1000:end-1000))
% plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
% title(['(c) Superimposed, VAF = ' num2str(round(Variance2,1)) '%'], 'Fontsize', 12)
% xlabel('Time (s)', 'Fontsize', 12)
% ylabel('Displacement (m)', 'Fontsize', 12)
% legend('Desired Displacement', 'Predicted Paralyzed Displacement, Pos_P(t)')
% grid on

%PDF and Spectrum of Residuals (between predicted healthy displacement and
%predicted paralyzed displacement)
residuals = predicted_healthy_displacement - predicted_paralyzed_displacement;

Nfft = 10000;
residuals_zero = residuals - mean(residuals);
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
semilogy(f(1:200,:),Prr(1:200,:))
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

subplot(2,2,3),semilogy(f(1:200,:),Pxx(1:200,:));
title('Predicted Healthy Displacement Spectrum','Fontsize',14);
ylabel('PSD (log scale)','Fontsize',16);
xlabel('Frequency (Hz)','Fontsize',16);
grid on;

subplot(2,2,4),semilogy(f(1:200,:),Pyy(1:200,:));
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
ylabel('Correlation','Fontsize',16)
grid on

[Ryy,lags] = xcorr(predicted_paralyzed_displacement_zero,maxlag,'coeff');
lags = lags/Fs;
subplot(3,1,2),plot(lags,Ryy)
title('(b) Auto-Correlation of Predicted Paralyzed Displacement','Fontsize',16)
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
txy = tfestimate(predicted_healthy_displacement,predicted_paralyzed_displacement,Nfft,[],Nfft,Fs);
G = abs(txy);
G = 20*log10(G);
P = angle(txy);
P = P*(180/pi);
[Cxy, f] = mscohere(predicted_healthy_displacement,predicted_paralyzed_displacement,Nfft,[],Nfft,Fs);

figure(figNum)
figNum = figNum+1;
sgtitle('Frequency Response of Predicted Healthy Displacement, Pos_H(t), and Predicted Paralyzed Displacement, Pos_P(t)','Fontsize',14)

subplot(3,1,1),plot(f(1:20,:),G(1:20,:));
ax = gca;
ax.FontSize = 14;
title('(a) Gain','Fontsize',16);
ylabel('Gain (dB)','Fontsize',16);
grid on;

subplot(3,1,2),plot(f(1:20,:),P(1:20,:));
ax = gca;
ax.FontSize = 14;
title('(b) Phase Shift','Fontsize',16);
ylabel('Phase (degrees)','Fontsize',16);
grid on;

subplot(3,1,3),plot(f(1:20,:),Cxy(1:20,:));
ax = gca;
ax.FontSize = 14;
title('(c) Coherence','Fontsize',16);
ylabel('Coherence','Fontsize',16);
xlabel('Frequency (Hz)','Fontsize',16);
grid on;

figure(figNum)
figNum = figNum+1;
plot(f(1:20,:),Cxy(1:20,:));
ax = gca;
ax.FontSize = 16;
title('Coherence between Pos_H(t) and Pos_P(t)','Fontsize',20);
ylabel('Coherence','Fontsize',20);
xlabel('Frequency (Hz)','Fontsize',20);
grid on;

% Frequency Response between Desired Displacement and
% Predicted Paralyzed Displacement
%Transfer Function
% Nfft = 1000;
% txy = tfestimate(desired_displacement(:,1000:end-1000),predicted_paralyzed_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);
% G = abs(txy);
% G = 20*log10(G);
% P = angle(txy);
% P = P*(180/pi);
% [Cxy, f] = mscohere(desired_displacement(:,1000:end-1000),predicted_paralyzed_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);
% 
% figure(figNum)
% figNum = figNum+1;
% sgtitle('Frequency Response of Desired Displacement and Predicted Paralyzed Displacement','Fontsize',14)
% 
% subplot(3,1,1),plot(f(1:31,:),G(1:31,:));
% ax = gca;
% ax.FontSize = 14;
% title('(a) Gain','Fontsize',16);
% ylabel('Gain (dB)','Fontsize',16);
% grid on;
% 
% subplot(3,1,2),plot(f(1:31,:),P(1:31,:));
% ax = gca;
% ax.FontSize = 14;
% title('(b) Phase Shift','Fontsize',16);
% ylabel('Phase (degrees)','Fontsize',16);
% grid on;
% 
% subplot(3,1,3),plot(f(1:31,:),Cxy(1:31,:));
% ax = gca;
% ax.FontSize = 14;
% title('(c) Coherence','Fontsize',16);
% ylabel('Coherence','Fontsize',16);
% xlabel('Frequency (Hz)','Fontsize',16);
% grid on;

% %Plotting the residuals as a funtion of desired displacement
% figure(figNum)
% figNum = figNum+1;
% plot(residuals, desired_displacement(:,1000:end-1000))
% title('Desired Displacement vs Residuals','Fontsize',20);
% ylabel('Displacement (m)','Fontsize',18);
% xlabel('Residual (m)','Fontsize',18);
% grid on;
% 
% %Plotting the residuals as a funtion of healthy displacement
% figure(figNum)
% figNum = figNum+1;
% plot(residuals, predicted_healthy_displacement)
% title('Predicted Healthy Displacement vs Residuals','Fontsize',20);
% ylabel('Displacement (m)','Fontsize',18);
% xlabel('Residual (m)','Fontsize',18);
% grid on;
% 
% %Plotting the residuals as a funtion of paralyzed displacement
% figure(figNum)
% figNum = figNum+1;
% plot(residuals, predicted_paralyzed_displacement)
% title('Predicted Paralyzed Displacement vs Residuals','Fontsize',20);
% ylabel('Displacement (m)','Fontsize',18);
% xlabel('Residual (m)','Fontsize',18);
% grid on;

tEnd = toc(tStart)/60
