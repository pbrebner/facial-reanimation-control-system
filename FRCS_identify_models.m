%% Facial Reanimation Control System (FRCS) - Identifiy ERS Model and SRS Model

% Identifies a Hammerstein Model for the healthy side and a Model for
% the paralyzed side. Can use both types of inputs (PRBS and Physiological)
% and have the option to compare the resulting model from each type of
% movement

%Since the physiologically based signal needs to be slightly different for
%each model type, the desired displacement signals are generated for each
%model type

clc
clear all

%% User Input Prompts

prompt1 = 'Type of SRS Model Structure? LNL/Hammerstein(Hamm)/Wiener/IRF [Wiener]: ';
str1 = input(prompt1,'s');
if ~strcmp(str1,'LNL') & ~strcmp(str1,'Hamm') & ~strcmp(str1,'Wiener') & ~strcmp(str1,'IRF') & ~isempty(str1)
    disp('Invalid Input')
    return
elseif isempty(str1)
    str1 = 'Wiener';
end
model_type = str1;

prompt2 = 'Compare two models (PRBS & Physiological)? Y/N [Y]: ';
str2 = input(prompt2,'s');
if ~strcmp(str2,'Y') & ~strcmp(str2,'N') & ~isempty(str2)
    disp('Invalid Input')
    return
elseif isempty(str2)
    str2 = 'Y';
end

if strcmp(str2,'Y')
    str2 = true;
elseif strcmp(str2,'N')
    str2 = false;
end

if str2 == false
    prompt3 = 'Type of Input? PRBS/Physiological(Phys) [PRBS]: ';
    str3 = input(prompt3,'s');
    if ~strcmp(str3,'PRBS') & ~strcmp(str3,'Phys') & ~isempty(str3)
        disp('Invalid Input')
        return
    elseif isempty(str3)
        str3 = 'PRBS';
    end
    
    if strcmp(str3,'PRBS')
        PRBS_signal = true;
    elseif strcmp(str3,'Phys')
        PRBS_signal = false;
    end 
    input_type = str3;
end

tStart = tic;

%% Set initial Parameters

%Noise Parameters
set_output_noise_power = 0;
noise_snr_ERS = [];
noise_snr_SRS = [];
output_noise_power = [];
figNum = 1;

%For Power Spectrums or FFT
Fs = 1000; 
Nfft = 10000;

%PRBS Signal Parameters
PRBS_movement_time = 180;
PRBS_stimulus_time = 480;
variable_amplitude_ERS = true;
variable_amplitude_SRS = true;
PRBS_amplitude_ERS = 10; %mm
PRBS_amplitude_SRS = 20; %mm

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_stimulus_time = 480;
physiological_stimulus_max_amplitude_ERS = 0.01;
physiological_stimulus_max_amplitude_SRS = 0.02;
fr = 0.1;                                         %Frequency distribution mean (Hz)
sig_ERS = 0.6;                                    %Std of Frequency Distribution (Hz)
sig_SRS = 0.8;                                    %Std of Frequency Distribution (Hz)
W_ERS = 0.55;
W_SRS = 0.45;
chance_of_zero = false;

%Initialize Models
ERS_all = [];
SRS_all = [];
%Identification Accuracy
accuracy_ERS = [];
accuracy_SRS = [];
%Predicted Displacement
pred_ERS = [];
pred_SRS = [];
%Input/Output for Identification
Zcur_ERS = [];
Zcur_SRS = [];

%Compare PRBS and Physiological Models
compare_two_models = str2;

if compare_two_models == true
    PRBS_signal = [true false];
end

%SRS Model Structure
LNL_model = false;
Hammerstein_model = false;
Weiner_model = false;
Linear_IRF_model = false;

if strcmp(str1,'LNL')
    LNL_model = true;
elseif strcmp(str1,'Hamm')
    Hammerstein_model = true;
elseif strcmp(str1,'Wiener')
    Weiner_model = true;
elseif strcmp(str1,'IRF')
    Linear_IRF_model = true;
end

%% Generate Desired Displacement Signals for ERS Model

for num_signals = 1:length(PRBS_signal)
    
    if PRBS_signal(num_signals) == true

        t_total = 0:0.001:PRBS_movement_time;
        time = PRBS_movement_time;
        
        N = PRBS_movement_time/10;
        M = 10000;

        A = [0];                            %Intialize Amplitude
        if variable_amplitude_ERS == true   
            for k = 1:N
                if k == 1
                    R = PRBS_amplitude_ERS;
                else
                    R = rand(1,1)*PRBS_amplitude_ERS;   %Randomly generate a number between 0 and max PRBS Amplitude
                end
                for j = 1:M
                    A = [A R];
                end
            end
        else
            A = PRBS_amplitude_ERS;              % Else set as Constant Amplitude
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

        FR = makedist('Normal','mu',fr,'sigma',sig_ERS);
        FrequenciesRandom_max = 1.8;
        FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
        freq_distribution = random(FrequenciesRandom,10000,1);

        AR = makedist('Uniform','lower',0,'upper',physiological_stimulus_max_amplitude_ERS);  %Full Amplitude Range
        AmplitudesRandom = AR;
        amp_distribution = random(AmplitudesRandom,10000,1);

        desired_displacement= 0;
        Freq_test = [];
        Pulses_per_interval_test = [];
        
        nf = time/10;                                   %number of random signal changes
        t_interval = physiological_movement_time/nf;    %Length of random interval (s)

        for j = 1 : nf    
            t  = 0 : 0.001 : t_interval;                %Time Intervals
            
            if j == 1
                Freq = FrequenciesRandom_max;
                A = physiological_stimulus_max_amplitude_ERS;
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
                data = (A*pulstran(t,D,@rectpuls,W_ERS))';
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
    
    %Neural Input Parameters Based on Desired Displacement Amplitude
    Amplitude = desired_displacement*100;  %mV
    Frequency = desired_displacement*14000;  %Hz
    
    neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';
    
    neural_simulink = [t_total' neural];
    
    figure(figNum)
    figNum = figNum+1;
    plot(t_total,neural)
    title('Neural Command Input')
    
    %% Execute ERS Simulation (Simulink Model)
    
    %Set Output Noise as zero
    set_param('ERS_simulation/Output Noise','Cov','set_output_noise_power')
    output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink;
    out = sim('ERS_simulation',time);

    %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
    %set_output_noise = set_output_noise_power;
    
    %% Get Output Signals from ERS Simulation
    
    %EMG, Muscle Force, and Healthy Displacement Output
    emg_simulink = out.ERS_Simulation_EMG;
    force_simulink = out.ERS_Simulation_Force;
    output_displacement_simulink = out.ERS_Simulation_Displacement;
    t_simulink = out.tout;
    
%     %Plot Healthy Displacement Output and Displacement Spectrum
%     output_displacement_simulink_zero = output_displacement_simulink - mean(output_displacement_simulink);
%     [Pxx,f] = pwelch(output_displacement_simulink_zero,[],[],[],Fs);
%     figure(figNum)
%     figNum = figNum+1;
%     subplot(2,1,1),plot(t_total,output_displacement_simulink)
%     title('Simulated Healthy Displacement, Pos_H(t)','Fontsize',20)
%     xlabel('Time (s)','Fontsize',20)
%     ylabel('Displacement (m)','Fontsize',20)
%     grid on
% 
%     subplot(2,1,2),semilogy(f(1:200,:),Pxx(1:200,:));
%     title('Healthy Displacement Spectrum','Fontsize',20);
%     ylabel('PSD','Fontsize',20);
%     xlabel('Frequency (Hz)','Fontsize',20);
%     grid on;
    
    %% Input/Output for ERS Identification
    
    Zcur = [emg_simulink,output_displacement_simulink];
    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});
    
    %% Calculate Signal to Noise Ratio
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    noise_snr_ERS = [noise_snr_ERS signal_to_noise];
    
    %% ERS Model Identification
    
    %Hammerstein Model(HK Method)
    set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
    NHK=nlbl;
    set(NHK,'idMethod','hk','displayFlag',true,'threshNSE',.001);

    % Set Number of lags in IRF
    I=NHK{1,2};
    set(I,'nLags',2400); 
    NHK{1,2}=I;
    
    %Identify and Normalize Linear Element
    NHK=nlident(NHK,Zcur);
    NHK = normCoefLE(NHK);
    
    figure(figNum);
    figNum = figNum+1;
    [R, V, yp] = nlid_resid(NHK,Zcur);
    
    accuracy_ERS = [accuracy_ERS V];
    pred_ERS = [pred_ERS double(yp)];

    ERS_all = [ERS_all NHK];
    Zcur_ERS = [Zcur_ERS Zcur];
    
end

% %Spectrum of Predicted Displacement
% pred_zero = double(yp) - mean(double(yp));
% [Pyy,f] = pwelch(pred_zero,[],[],[],Fs);
% figure(figNum)
% figNum = figNum+1;
% subplot(2,1,1),plot(t_total,double(yp))
% title('Predicted Healthy Displacement for ERS Identification, Pos_H(t)','Fontsize',20)
% xlabel('Time (s)','Fontsize',20)
% ylabel('Displacement (m)','Fontsize',20)
% grid on
% 
% subplot(2,1,2),semilogy(f(1:200,:),Pyy(1:200,:));
% title('Predicted Healthy Displacement Spectrum','Fontsize',20);
% ylabel('PSD','Fontsize',20);
% xlabel('Frequency (Hz)','Fontsize',20);
% grid on;

%% Generate Desired Displacement Signals for SRS Model Identification

for num_signals = 1:length(PRBS_signal)
    
    if PRBS_signal(num_signals) == true

        t_total = 0:0.001:PRBS_stimulus_time;
        time = PRBS_stimulus_time;
        
        N = PRBS_stimulus_time/10;
        M = 10000;

        A = [0];                                        %Intialize amplitude
        if variable_amplitude_SRS == true   
            for k = 1:N
                if k == 1
                    R = PRBS_amplitude_SRS;
                else
                    R = rand(1,1)*PRBS_amplitude_SRS;   %Randomly generate a number between 0 and PRBS Amplitude
                end
                for j = 1:M
                    A = [A R];
                end
            end
        else
            A = PRBS_amplitude_SRS;              %Else Set as Constant Amplitude
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
        
        t_total = 0:0.001:physiological_stimulus_time;
        time = physiological_stimulus_time;
        
        FR = makedist('Normal','mu',fr,'sigma',sig_SRS);
        FrequenciesRandom_max = 2.21;
        FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
        freq_distribution = random(FrequenciesRandom,10000,1);

        AR = makedist('Uniform','lower',0,'upper',physiological_stimulus_max_amplitude_SRS);  %Full Amplitude Range
        AmplitudesRandom = AR;
        amp_distribution = random(AmplitudesRandom,10000,1);

        desired_displacement= 0;
        Freq_test = [];
        Pulses_per_interval_test = [];
        
        nf = time/10;                                   %number of random signal changes
        t_interval = physiological_stimulus_time/nf;    %Length of random interval (s)

        for j = 1 : nf    
            t  = 0 : 0.001 : t_interval;         %Time Intervals
            
            if j == 1
                Freq = FrequenciesRandom_max;
                A = physiological_stimulus_max_amplitude_SRS;
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
                data = (A*pulstran(t,D,@rectpuls,W_SRS))';
                stim_frequency = Freq;
                Pulses_per_interval = t_interval/g;
                data(end) = [];
            else
                data = zeros(20000,1);
                stim_frequency = 0;
                Pulses_per_interval = 0;
            end

            %data = A*square(2*pi*f/Fs*t)';        % Generate Square Wave
            desired_displacement = [desired_displacement; data];
            Freq_test = [Freq_test stim_frequency];
            Pulses_per_interval_test = [Pulses_per_interval_test Pulses_per_interval];
        end

        Pulses_per_interval_total = sum(Pulses_per_interval_test);
        Freq_test_average = sum(Freq_test)/length(Freq_test);

        desired_displacement = desired_displacement';

    end
    
    %% Generate Amplitude Modulation Input for SRS Simulation
    
    stim_frequency = 50;
    stim_amplitude = desired_displacement*170;

    %Theoretical Electrical Stimulus
    input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);
    
    %Amplitude Modulation Input
    amplitude_modulation = stim_amplitude;
    amplitude_modulation_simulink = [t_total' amplitude_modulation'];
    
    figure(figNum)
    figNum = figNum+1;
    plot(t_total,amplitude_modulation)
    title('Amplitude Modulation Input')
    
    %% Execute SRS Simulation (Simulink Model)
    
    %Set Output Noise as zero
    set_param('SRS_simulation/Output Noise','Cov','set_output_noise_power')
    output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink;
    out = sim('SRS_simulation',time);

    %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
    %set_output_noise = set_output_noise_power;
    
    %% Get Output Signals from SRS Simulation

    %Muscle Force, Electrical Stimulus, and Paralyzed Displacement Output
    force_simulink = out.SRS_Simulation_Force;
    input_stimulus = out.SRS_Simulation_Stimulus;
    output_displacement_simulink = out.SRS_Simulation_Displacement;
    t_simulink = out.tout;
    
%     %Plot Paralyzed Displacement Output and Spectrum
%     output_displacement_simulink_zero = output_displacement_simulink - mean(output_displacement_simulink);
%     [Pxx,f] = pwelch(output_displacement_simulink_zero,[],[],[],Fs);
%     figure(figNum)
%     figNum = figNum+1;
%     subplot(2,1,1),plot(t_total,output_displacement_simulink)
%     title('Simulated Paralyzed Displacement, Pos_P(t)','Fontsize',20)
%     xlabel('Time (s)','Fontsize',20)
%     ylabel('Displacement (m)','Fontsize',20)
%     grid on
% 
%     subplot(2,1,2),semilogy(f(1:300,:),Pxx(1:300,:));
%     title('Paralyzed Displacement Spectrum','Fontsize',20);
%     ylabel('PSD','Fontsize',20);
%     xlabel('Frequency (Hz)','Fontsize',20);
%     grid on;

    %% Input/Output for SRS Identification
    
    Zcur = [amplitude_modulation',output_displacement_simulink];
    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Stimulus, Output Displacement','chanNames', {'Stimulus (V)' 'Displacement (m)'});
    
    %% Calculate Signal to Noise Ratio
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    noise_snr_SRS = [noise_snr_SRS signal_to_noise];
    
    %% SRS Model Identification (LNL, Hammerstein, Wiener, or Linear IRF)
    
    set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
    
    if LNL_model == true
        LNL=lnlbl;
        set(LNL,'idMethod','hk','nLags1',750,'polyOrderMax',2,'nLags2',750,... 
        'hkTolerance', 0.1, 'nhkMaxIts', 4, 'nhkMaxInner', 4);

        LNL=nlident(LNL,Zcur);

        figure(figNum)
        figNum = figNum+1;
        [R, V, yp] = nlid_resid(LNL,Zcur);

        accuracy_SRS = [accuracy_SRS V];
        pred_SRS = [pred_SRS double(yp)];

        SRS_all = [SRS_all LNL];
        Zcur_SRS = [Zcur_SRS Zcur];
        
    elseif Hammerstein_model == true
        
        Hammerstein=nlbl;
        set(Hammerstein,'idMethod','hk','displayFlag',true,'threshNSE',.001);
        I=Hammerstein{1,2};
        set(I,'nLags',2400,'nSides',1);
        Hammerstein{1,2}=I;
        
        Hammerstein=nlident(Hammerstein,Zcur);
        
        accuracy_SRS = [accuracy_SRS V];
        pred_SRS = [pred_SRS double(yp)];
        
        SRS_all = [SRS_all Hammerstein];
        Zcur_SRS = [Zcur_SRS Zcur];
        
    elseif Weiner_model == true
        
        Weiner = lnbl; %Wiener
        set(Weiner,'idMethod','hk');
        I = Weiner{1,1};
        set(I,'nLags',850,'nSides',1); % Set Number of lags and Sides in IRF
        Weiner{1,1} = I;
        
        Weiner=nlident(Weiner,Zcur);
        
        figure(figNum)
        figNum = figNum+1;
        [R, V, yp] = nlid_resid(Weiner,Zcur);
        
        accuracy_SRS = [accuracy_SRS V];
        pred_SRS = [pred_SRS double(yp)];
        
        SRS_all = [SRS_all Weiner];
        Zcur_SRS = [Zcur_SRS Zcur];
        
    elseif Linear_IRF_model == true
        
        IRF_model = irf(Zcur,'nLags',1000,'nSides',1);

        figure(figNum)
        figNum = figNum+1;
        [R, V, yp] = nlid_resid(IRF_model,Zcur);
        
        accuracy_SRS = [accuracy_SRS V];
        pred_SRS = [pred_SRS double(yp)];

        if compare_two_models == true && PRBS_signal(num_signals) == true
            IRF1 = IRF_model;
        elseif compare_two_models == true && PRBS_signal(num_signals) == false
            IRF2 = IRF_model;
        end
        
        Zcur_SRS = [Zcur_SRS Zcur];
        
    end
    
%     %Spectrum of Predicted Displacement
%     pred_zero = double(yp) - mean(double(yp));
%     [Pyy,f] = pwelch(pred_zero,[],[],[],Fs);
%     figure(figNum)
%     figNum = figNum+1;
%     subplot(2,1,1),plot(t_total,double(yp))
%     title('Predicted Displacement for SRS Identification','Fontsize',20)
%     xlabel('Time (s)','Fontsize',20)
%     ylabel('Displacement (m)','Fontsize',20)
%     grid on
% 
%     subplot(2,1,2),semilogy(f(1:200,:),Pyy(1:200,:));
%     title('Predicted Displacement Spectrum','Fontsize',20);
%     ylabel('PSD','Fontsize',20);
%     xlabel('Frequency (Hz)','Fontsize',20);
%     grid on;
    
end
%% Plots the Models (Compares the PRBS and Physiological Models)

%For ERS Model Non-linear Element (Limits to plot only "High Quality"
%Values)
upper_limit = 0.000586;
lower_limit = -0.00063;
    
if compare_two_models == true && LNL_model == true
    
    NHK1 = ERS_all(1);
    NHK2 = ERS_all(2);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(NHK1{1,1})
    hold on
    plot(NHK2{1,1})
    ax = gca;
    ax.FontSize = 15;
    xlim([lower_limit upper_limit])
    title('(a) Static Nonlinearities','Fontsize', 20)
    xlabel('EMG Input (V)','Fontsize',18)
    ylabel('Displacement Output (m)','Fontsize',18)
    grid on
    hold off

    subplot(1,2,2)
    plot(NHK1{1,2})
    hold on
    plot(NHK2{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Linear Elements','Fontsize', 20)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    legend('PRBS', 'Physiological','Fontsize',14)
    grid on
    hold off
    
    LNL1 = SRS_all(1);
    LNL2 = SRS_all(2);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,3,1)
    hold on
    plot(LNL1{1,1}) 
    plot(LNL2{1,1})
    hold off
    ax = gca;
    ax.FontSize = 13;
    title('(a) Input Linear Element','Fontsize', 16)
    xlabel('Lags (s)','Fontsize',16)
    ylabel('X1','Fontsize',16)
    grid on

    subplot(1,3,2)
    hold on
    plot(LNL1{1,2})
    plot(LNL2{1,2})
    hold off
    ax = gca;
    ax.FontSize = 13;
    title('(b) Static Nonlinearity','Fontsize', 16)
    xlabel('Amplitude Modulation Input','Fontsize',16)
    ylabel('Displacement Output (m)','Fontsize',16)
    grid on

    subplot(1,3,3)
    hold on
    plot(LNL1{1,3}); 
    plot(LNL2{1,3});
    hold off
    ax = gca;
    ax.FontSize = 13;
    title('(c) Output Linear Element','Fontsize', 16)
    xlabel('Lags (s)','Fontsize',16)
    ylabel('X1','Fontsize',16)
    legend('PRBS', 'Physiological','Fontsize',14)
    grid on
    
elseif compare_two_models == true && Hammerstein_model == true
    
    NHK1 = ERS_all(1);
    NHK2 = ERS_all(2);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(NHK1{1,1})
    hold on
    plot(NHK2{1,1})
    ax = gca;
    ax.FontSize = 15;
    xlim([lower_limit upper_limit])
    title('(a) Static Nonlinearities','Fontsize', 20)
    xlabel('EMG Input (V)','Fontsize',18)
    ylabel('Displacement Output (m)','Fontsize',18)
    grid on
    hold off

    subplot(1,2,2)
    plot(NHK1{1,2})
    hold on
    plot(NHK2{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Linear Elements','Fontsize', 20)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    legend('PRBS', 'Physiological','Fontsize',14)
    grid on
    hold off
    
    Hammerstein1 = SRS_all(1);
    Hammerstein2 = SRS_all(2);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(Hammerstein1{1,1})
    hold on
    plot(Hammerstein2{1,1})
    ax = gca;
    ax.FontSize = 15;
    title('(a) Static Nonlinearities','Fontsize', 20)
    xlabel('Amplitude Modulation Input','Fontsize',18)
    ylabel('Displacement Output (m)','Fontsize',18)
    grid on
    hold off

    subplot(1,2,2)
    plot(Hammerstein1{1,2})
    hold on
    plot(Hammerstein2{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Linear Elements','Fontsize', 20)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    legend('PRBS', 'Physiological','Fontsize', 14)
    grid on
    hold off
    
elseif compare_two_models == true && Weiner_model == true
    
    NHK1 = ERS_all(1);
    NHK2 = ERS_all(2);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(NHK1{1,1})
    hold on
    plot(NHK2{1,1})
    ax = gca;
    ax.FontSize = 15;
    xlim([lower_limit upper_limit])
    title('(a) Static Nonlinearities','Fontsize', 20)
    xlabel('EMG Input, E(t) (V)','Fontsize',18)
    ylabel('Transformed EMG (V)','Fontsize',18)
    grid on
    hold off

    subplot(1,2,2)
    plot(NHK1{1,2})
    hold on
    plot(NHK2{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Linear Elements','Fontsize', 20)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    legend('PRBS', 'Physiological','Fontsize',14)
    grid on
    hold off
    
    Weiner1 = SRS_all(1);
    Weiner2 = SRS_all(2);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(Weiner1{1,1})
    hold on
    plot(Weiner2{1,1})
    ax = gca;
    ax.FontSize = 15;
    title('(a) Linear Elements','Fontsize', 20)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    grid on
    hold off

    subplot(1,2,2)
    plot(Weiner1{1,2})
    hold on
    plot(Weiner2{1,2})
    ax = gca;
    ax.FontSize = 15;
    %xlim([lower_limit upper_limit])
    title('(b) Static Nonlinearities','Fontsize', 20)
    xlabel('Transformed Displacement Input (m)','Fontsize',18)
    ylabel('Paralyzed Displacement Output, Pos_P(t) (m)','Fontsize',18)
    legend('PRBS', 'Physiological','Fontsize', 14)
    grid on
    hold off
    
elseif compare_two_models == true && Linear_IRF_model == true
    
    NHK1 = ERS_all(1);
    NHK2 = ERS_all(2);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(NHK1{1,1})
    hold on
    plot(NHK2{1,1})
    ax = gca;
    ax.FontSize = 15;
    xlim([lower_limit upper_limit])
    title('(a) Static Nonlinearities','Fontsize', 20)
    xlabel('EMG Input','Fontsize',18)
    ylabel('Displacement Output (m)','Fontsize',18)
    grid on
    hold off

    subplot(1,2,2)
    plot(NHK1{1,2})
    hold on
    plot(NHK2{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Linear Elements','Fontsize', 20)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    legend('PRBS', 'Physiological','Fontsize',14)
    grid on
    hold off
    
    figure(figNum)
    figNum = figNum+1;
    hold on
    plot(IRF1) 
    plot(IRF2)
    hold off
    ax = gca;
    ax.FontSize = 14;
    title('Linear IRF Models','Fontsize', 20)
    xlabel('Lags (s)','Fontsize',20)
    ylabel('X1','Fontsize',20)
    legend('PRBS', 'Physiological','Fontsize',14)
    grid on

    
elseif compare_two_models == false && LNL_model == true
    
    NHK1 = ERS_all(1);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(NHK1{1,1})
    ax = gca;
    ax.FontSize = 15;
    title('(a) Static Nonlinearities','Fontsize', 20)
    xlabel('EMG Input (V)','Fontsize',18)
    ylabel('Displacement Output (m)','Fontsize',18)
    grid on

    subplot(1,2,2)
    plot(NHK1{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Linear Elements','Fontsize', 20)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    grid on
    
    LNL1 = SRS_all(1);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,3,1)
    plot(LNL1{1,1}) 
    ax = gca;
    ax.FontSize = 13;
    title('(a) Input Linear Element','Fontsize', 16)
    xlabel('Lags (s)','Fontsize',16)
    ylabel('X1','Fontsize',16)
    grid on

    subplot(1,3,2)
    plot(LNL1{1,2})
    ax = gca;
    ax.FontSize = 13;
    title('(b) Static Nonlinearity','Fontsize', 16)
    xlabel('Amplitude Modulation Input','Fontsize',16)
    ylabel('Displacement Output (m)','Fontsize',16)
    grid on

    subplot(1,3,3)
    plot(LNL1{1,3}); 
    ax = gca;
    ax.FontSize = 13;
    title('(c) Output Linear Element','Fontsize', 16)
    xlabel('Lags (s)','Fontsize',16)
    ylabel('X1','Fontsize',16)
    grid on
    
elseif compare_two_models == false && Weiner_model == true
    
    NHK1 = ERS_all(1);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(NHK1{1,1})
    ax = gca;
    ax.FontSize = 15;
    title('(a) Static Nonlinearities','Fontsize', 20)
    xlabel('EMG Input (V)','Fontsize',18)
    ylabel('Displacement Output (m)','Fontsize',18)
    grid on

    subplot(1,2,2)
    plot(NHK1{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Linear Elements','Fontsize', 20)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    grid on
    
    Weiner1 = SRS_all(1);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(Weiner1{1,1})
    ax = gca;
    ax.FontSize = 15;
    title('(a) Linear Elements','Fontsize', 20)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    grid on

    subplot(1,2,2)
    plot(Weiner1{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Static Nonlinearities','Fontsize', 20)
    xlabel('Amplitude Modulation Input','Fontsize',18)
    ylabel('Displacement Output (m)','Fontsize',18)
    legend('PRBS', 'Physiological','Fontsize', 14)
    grid on
    
else
    
    NHK1 = ERS_all(1);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(NHK1{1,1})
    ax = gca;
    ax.FontSize = 15;
    title('(a) Static Nonlinearities','Fontsize', 20)
    xlabel('EMG Input (V)','Fontsize',18)
    ylabel('Displacement Output (m)','Fontsize',18)
    grid on

    subplot(1,2,2)
    plot(NHK1{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Linear Elements','Fontsize', 20)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    grid on
    
    figure(figNum)
    figNum = figNum+1;
    plot(IRF_model) 
    ax = gca;
    ax.FontSize = 14;
    title('Linear IRF Model','Fontsize',20)
    xlabel('Lags (s)','Fontsize',20)
    ylabel('X1','Fontsize',20)
    grid on
    
end

tEnd = toc(tStart)/60