%% ERS Model Test (MODEL VALIDATION)

%Validates the estimated models from ERS_model by using them to predict the
%movements generated in response to input realizations which are different 
%than that used to identify the models. Will validate the specified model 
%with a PRBS input signal and then a "Physiological" input signal. Plots the 
%superimposed predicted and observed displacements, residual distribution, 
%and residual spectrum 

%INSTRUCTIONS: RUN ERS_model FIRST and THEN RUN THIS

%When running the script, you need to provide the following input:
%If both PRBS and Physiological Models were Identified in ERS_model,
% 1. Which Identified Model Type do you want to Validate?
%       Choose either the PRBS model or Phys model for validation
%If not, will skip question 1 and default to the only model identified.
% 2. Number of Validation Trials?
%       Number of validation trials used to calculate VAF mean and std

%% User Input Prompts

if compare_two_models == true
    prompt1 = 'Which Identified Model do you want to Validate? PRBS/Phys [PRBS]: ';
    str1 = input(prompt1,'s');
    if ~strcmp(str1,'PRBS') & ~strcmp(str1,'Phys') & ~isempty(str1)
        disp('Invalid Input')
        return
    elseif isempty(str1)
        str1 = 'PRBS';
    end
    
    if strcmp(str1,'PRBS')
        NHK = NHK1;
        disp('Validating Model Identified from PRBS Input')
    elseif strcmp(str1,'Phys')
        NHK = NHK2;
        disp('Validating Model Identified from Physiological Input')
    end
else
    NHK = NHK_all(1);
    disp(['Validating Model Identified from ' input_type ' Input'])
end

prompt2 = 'Number of Validation Trials? 1-50 [1]: ';
str2 = input(prompt2);
if str2<1 | str2>50
    disp('Invalid Input')
    return
elseif isempty(str2)
    str2 = 1;
end

tStart = tic;

%% Set initial Parameters

%Noise Parameters
noise_snr = [];
set_output_noise_power = 0;
output_noise_power = [];
figNum = 100;

%For Power Spectrums
Fs = 1000; 
Nfft = 10000;

%PRBS Signal Parameters
PRBS_movement_time = 180;
variable_amplitude = true;      %PRBS can either be constant amplitude or variable amplitude
N = PRBS_movement_time/10;      %Number of times the amplitude randomly changes (For Variable Amplitude Only)
M = 10000;                      %Number of each random value (For Variable Amplitude Only)
PRBS_amplitude = 10;            %PRBS Amplitude (mm)

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_movement_max_amplitude = 0.01;    %"Physioloigcal" Amplitude (m)
fr = 0.1;                                       %Frequency distribution mean (Hz) (Max is 1.8 Hz)
sig = 0.6;                                      %Std of Frequency Distribution (Hz)
W = 0.55;                                       %Width of signal pulse (seconds) 
nf = physiological_movement_time/10;            %Number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (seconds)
chance_of_zero = false;

%Set number of Validation Trials (To get a mean and std VAF)
num_trials = str2;
validation_accuracy = [];
validation_accuracy_all = [];
Zcur_all_val = [];

%for validation we test with PRBS signal first, and then with Physiological signal
PRBS_movement = [true false];

%% Generate Desired Displacement Signals for Model Validation

for signal = 1:length(PRBS_movement)
    
    for trial = 1:num_trials

        if PRBS_movement(signal) == true
            t_total = 0:0.001:PRBS_movement_time;
            time = PRBS_movement_time;

            A = 0;                                      %Intialize amplitude
            if variable_amplitude == true   
                for k = 1:N
                    if k == 1
                        R = PRBS_amplitude;             %First interval is at max PRBS Amplitude
                    else
                        R = rand(1,1)*PRBS_amplitude;   %Randomly generate a number between 0 and PRBS Amplitude
                    end
                    for j = 1:M
                        A = [A R];
                    end
                end
            else
                A = PRBS_amplitude;                      %Else set as Constant Amplitude
            end

            Range = [0,0.001]; %Specify what the single-channel PRBS value switches between

            %Specify the clock period of the signal as 1 sample. 
            %That is, the signal value can change at each time step. 
            %For PRBS signals, the clock period is specified in Band = [0 B], 
            %where B is the inverse of the required clock period
            %(Must be less than 1)
            Band = [0 0.01];

            %Generate a nonperiodic PRBS of length time.
            u = idinput(time*1000+1,'prbs',Band,Range);

            %Create an iddata object from the generated signal. 
            %For this example, specify the sample time as 0.001 second.
            u = iddata([],u,0.001);

            U = (u.InputData)';
            desired_displacement = A.*U;

            %Power Pectrum of Desired Displacement
            desired_displacement_zero = desired_displacement - mean(desired_displacement);
            [Pxx1,f1] = pwelch(desired_displacement_zero,Nfft,[],Nfft,Fs);

        else

            % "Physiological" Movements
            t_total = 0:0.001:physiological_movement_time;
            time = physiological_movement_time;

            FR = makedist('Normal','mu',fr,'sigma',sig);
            FrequenciesRandom_max = 1.8;
            FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
            freq_distribution = random(FrequenciesRandom,10000,1);

            AR = makedist('Uniform','lower',0,'upper',physiological_movement_max_amplitude);  %Full Amplitude Range
            AmplitudesRandom = AR;
            amp_distribution = random(AmplitudesRandom,10000,1);

            desired_displacement= 0;
            Freq_test = [];
            Pulses_per_interval_test = [];

            for j = 1 : nf    
                t  = 0 : 0.001 : t_interval;         % Length of Time Interval
                if j == 1
                    Freq = FrequenciesRandom_max;
                    A = physiological_movement_max_amplitude;
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
                    movement_frequency = Freq;
                    Pulses_per_interval = t_interval/g;
                    data(end) = [];
                else
                    data = zeros(20000,1);
                    movement_frequency = 0;
                    Pulses_per_interval = 0;
                end

                %data = A*square(2*pi*f/Fs*t)';        % Generate Square Wave
                desired_displacement = [desired_displacement; data];
                Freq_test = [Freq_test movement_frequency];
                Pulses_per_interval_test = [Pulses_per_interval_test Pulses_per_interval];
            end

            Pulses_per_interval_total = sum(Pulses_per_interval_test);
            Freq_test_average = sum(Freq_test)/length(Freq_test);

            %Power Pectrum of Desired Displacement
            desired_displacement_zero = desired_displacement - mean(desired_displacement);
            [Pxx1,f1] = pwelch(desired_displacement_zero,Nfft,[],Nfft,Fs);

            desired_displacement = desired_displacement';

        end

        %% Create Neural Input for ERS Simulation

        %Create Frequency and Amplitude Parameters of the Neural Input based on
        %Desired Displacement Amplitude
        Amplitude = desired_displacement*100;    %mV
        Frequency = desired_displacement*14000;  %Hz

        %Generate Neural Command Signal with Frequency and Amplitdue Parameters
        neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

        neural_simulink = [t_total' neural]; %Input for the Simulink Model

        %Power Pectrum of Neural Input
        neural_zero = neural - mean(neural);
        [Pxx2,f2] = pwelch(neural_zero,Nfft,[],Nfft,Fs);

        %% Execute ERS Simulation Simulink Model

        %Set Output Noise (set as zero)
        set_param('ERS_simulation/Output Noise','Cov','set_output_noise_power')
        output_noise_power = [output_noise_power set_output_noise_power];

        %Run Simulink Model
        out = sim('ERS_simulation',time);

        %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
        %set_output_noise = set_output_noise_power;

        %% Get Output Signals from ERS Simulation Model

        %EMG, Muscle Force, and Healthy Displacement Output
        emg_simulink = out.ERS_Simulation_EMG;
        force_simulink = out.ERS_Siomulation_Force;
        output_displacement_simulink = out.ERS_Simulation_Displacement;
        t_simulink = out.tout;

        %% Input/Output for Model Validation

        Zcur = [emg_simulink,output_displacement_simulink];
        Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});

        %% Calculate Signal to Noise Ratio
        output_noise_simulink = out.Output_Noise;

        signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
        noise_snr = [noise_snr signal_to_noise];

        %% Model Validation
        
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

        %Plot Results of Validation
        figure(figNum)
        [R, V, yp] = nlid_resid(NHK,Zcur);
        
        %validation accuracy for all trials of current signal type
        validation_accuracy = [validation_accuracy V];
        
        Zcur_all_val = [Zcur_all_val Zcur]; 

    end
    
    %Store the validation accuracy for the signal type and reset for next
    %signal iteration
    validation_accuracy_all(signal,:) = validation_accuracy;
    validation_accuracy = [];

    %% Plot the Validation results of last validation trial

    %Plots the Superimposed Predicted and Observed Healthy Displacement with %VAF in the title
    figNum = figNum+1;
    figure(figNum);
    figNum = figNum+1;
    pred = double(yp);
    plot(t_total,pred);
    hold on
    plot(t_total, output_displacement_simulink)
    ax = gca;
    ax.FontSize = 14;
    hold off
    title(['Superimposed, VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 24)
    xlabel('Time (s)', 'Fontsize', 18)
    ylabel('Displacement, Pos_H(t) (m)', 'Fontsize', 18)
    legend('Predicted', 'Observed','Fontsize',14)

    %Zero mean of Residuals and Spectrum of Residuals
    R_zero = R - mean(R);
    [PxxR,fR] = pwelch(double(R_zero),[],[],[],Fs);

    %Plots the Superimposed Predicted and Observed Outputs, Residual
    %Distribution, and Residual Spectrum
    figure(figNum)
    figNum = figNum+1;
    subplot(2,2,[1 2])
    plot(t_total,pred);
    hold on
    plot(t_total, output_displacement_simulink)
    ax = gca;
    ax.FontSize = 12;
    hold off
    title(['(a) Superimposed, VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 18)
    xlabel('Time (s)', 'Fontsize', 16)
    ylabel('Displacement, Pos_H(t) (m)', 'Fontsize', 16)
    legend('Predicted', 'Observed','Fontsize', 14)

    subplot(2,2,3)
    plot(p)
    ax = gca;
    ax.FontSize = 12;
    xlabel('Displacement (m)','Fontsize',16)
    ylabel('Density','Fontsize',16)
    title('(b) Residual Distribution','Fontsize',18)
    grid on

    subplot(2,2,4)
    semilogy(fR(1:650,:),PxxR(1:650,:),'LineWidth',1.5);
    ax = gca;
    ax.FontSize = 12;
    title('(c) Power Spectrum of Residuals','Fontsize',18);
    ylabel('PSD (log)','Fontsize',16); 
    xlabel('Frequency (Hz)','Fontsize',16);
    grid on

    %Plots the Superimposed Predicted and Observed Outputs, Residual
    %Distribution, and Residual Auto-Correlation
    figure(figNum)
    figNum = figNum+1;
    subplot(2,2,[1 2])
    plot(t_total,pred);
    hold on
    plot(t_total, output_displacement_simulink)
    ax = gca;
    ax.FontSize = 12;
    hold off
    title(['(a) Superimposed, VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 18)
    xlabel('Time (s)', 'Fontsize', 16)
    ylabel('Displacement, Pos_H(t) (m)', 'Fontsize', 16)
    legend('Predicted', 'Observed','Fontsize', 14)

    subplot(2,2,3)
    plot(p)
    ax = gca;
    ax.FontSize = 12;
    xlabel('Displacement (m)','Fontsize',16)
    ylabel('Density','Fontsize',16)
    title('(b) Residual Distribution','Fontsize',18)
    grid on

    subplot(2,2,4)
    R_double = double(R_zero);
    L = length(R_double)-1;
    maxlag = L*0.10;
    [res_corr,lags] = xcorr(R_double,maxlag,'normalized');
    lags = lags/Fs;
    subplot(2,2,4),plot(lags,res_corr,'LineWidth',1.5)
    ax = gca;
    ax.FontSize = 12;
    title('(c) Auto-Correlation Residuals','Fontsize',18)
    xlabel('Time (s)','Fontsize',16)
    ylabel('Correlation','Fontsize',16)
    grid on

    %Plots the Residuals vs Output Displacement
%     figure(figNum)
%     figNum = figNum+1;
%     Rx = double(R);
%     plot(Rx, output_displacement_simulink)
%     ax = gca;
%     ax.FontSize = 14;
%     title('Residuals vs Displacement Amplitude','Fontsize',24);
%     ylabel('Displacement Amplitude (m)','Fontsize',18); 
%     xlabel('Residuals (m)','Fontsize',18);
%     grid on

end

%% Calculate Validation Mean and Std for both Input Types
validation_accuracy_mean_PRBS = mean(validation_accuracy_all(1,:),2);
disp(['Validation Accuracy Mean for PRBS Input: ' num2str(round(validation_accuracy_mean_PRBS,1)) '%'])

validation_accuracy_std_PRBS = std(validation_accuracy_all(1,:),0,2);
disp(['Validation Accuracy Std for PRBS Input: ' num2str(round(validation_accuracy_std_PRBS,1)) '%'])

validation_accuracy_mean_phys = mean(validation_accuracy_all(2,:),2);
disp(['Validation Accuracy Mean for Physiological Input: ' num2str(round(validation_accuracy_mean_phys,1)) '%'])

validation_accuracy_std_phys = std(validation_accuracy_all(2,:),0,2);
disp(['Validation Accuracy Std for Physiological Input: ' num2str(round(validation_accuracy_std_phys,1)) '%'])

tEnd = toc(tStart)/60