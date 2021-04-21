%EMG Model Test

%Validate the estimated model by using it to predict the movements
%generated in response to an input realization which is different than that
%used to identify the model

%Tests the Model based on the number of pulses of the Physiological Input

%Identifies models with one set of signals and then validates with multiple
%different signals

clc
clear all

tStart = tic;

%%
%Set initial Parameters
set_output_noise_power = 0;
%noise_multiplier = 10;
noise_snr = [];
output_noise_power = [];
% figNum = 1;

%For Power Spectrums
Fs = 1000; 
Nfft = 1000;

%Desired Displacement Signal Type (Simple, PRBS or Physiologically Based)
% simple_movement = false;
% PRBS_movement = true;
%physiological_movement = true;

%Signal Amplitude (PRBS Signal)
% PRBS_movement_time = 180;
% variable_amplitude = false;
% N = 18;                     %Number of times the amplitude randomly changes (For Variable Amplitude Only)
% M = 10000;                  %Number of each random value (For Variable Amplitude Only)
% PRBS_amplitude = 10;             %Amplitude (For Constant Amplitude Only)

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_movement_max_amplitude = 0.01;
fr = 0.1;                 %Frequency distribution mean (Hz) (Max is 1.8 Hz)
sig = 0.6;                %Std of Frequency Distribution (Hz)
W = 0.55;                  
nf = 18;                  %number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (s)
chance_of_zero = false;

NHK_all_temp = [];
Zcur_all_temp = [];
NHK_all = [];
Zcur_all = [];

identification_accuracy = [];
validation_accuracy = [];

%%
use_fr = false; %Use frequency instead of number of pulses for Physiological Input
if use_fr == true
    fr = 0:0.1:0.3;
    sig = 0.6;
end

num_pulses = 1:1:20;
num_trials_ident = 10;
num_trials_val = 30;
variable_time = false;

if use_fr == true
    num_models_for_ident = length(fr);
else
    num_models_for_ident = length(num_pulses);
end

validation_frequency = [];
validation_pulses_total = [];

%%
for trial = 1:num_trials_ident
    
    for model = 1:num_models_for_ident
        
        if variable_time == true 
            t_total = 0:0.001:((physiological_movement_time/num_models_for_ident)*model);
            time = ((physiological_movement_time/num_models_for_ident)*model);
        else
            t_total = 0:0.001:physiological_movement_time;
            time = physiological_movement_time;
        end

        if use_fr == true
            FR = makedist('Normal','mu',fr(model),'sigma',sig);
            FrequenciesRandom_max = 1.8;
            FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
            freq_distribution = random(FrequenciesRandom,10000,1);
        %     figure(figNum)
        %     figNum = figNum+1;
        %     histogram(freq_distribution,100)
        %     title('Frequency Distribution')
        %     xlabel('Frequency (Hz)')
        end

        AR = makedist('Uniform','lower',0,'upper',physiological_movement_max_amplitude);  %Full Amplitude Range
        AmplitudesRandom = AR;
        amp_distribution = random(AmplitudesRandom,10000,1);
        %figure(figNum)
        %figNum = figNum+1;
        %histogram(amp_distribution,100)
        %title('Amplitude Distribution')
        %xlabel('Amplitude (mm)')

        %PWR = makedist('Normal','mu',0.001,'sigma',0.005);
        %PulseWidthRandom = truncate(PWR,0.001,0.010);   % 1ms - 10 ms
        %pw_distribution = random(PulseWidthRandom,10000,1);
        %figure(107)
        %histogram(pw_distribution,100)
        %title('Pulse Width Distribution')

        desired_displacement= 0;
        Freq_test = [];
        Pulses_per_interval_test = [];

        if use_fr == true

            for j = 1 : nf    
                t  = 0 : 0.001 : t_interval;         % Time Samples

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

        else
            
            t  = 0 : 0.001 : time;

            %A = random(AmplitudesRandom,1,1);
            A = 0.01;
            g = time/num_pulses(model);

            DS = makedist('Uniform','lower',0,'upper',g);
            DelayStart = DS;
            delay_start_distribution = random(DelayStart,10000,1);

            d = random(DelayStart,1,1);

            Delay = d;

            for hh = 1:model-1

                %normal distribution from 0 - g
                DR = makedist('Uniform','lower',0,'upper',g);
                DelaysRandom = DR;
                delay_random_distribution = random(DelaysRandom,10000,1);
    %             figure(figNum)
    %             figNum = figNum+1;
    %             histogram(delay_random_distribution,100)

                random_delay = random(DelaysRandom,1,1);

                Delay = [Delay hh*g+random_delay];
                %Delay = [Delay Delay(end)+random_delay];

            end

            %Delay = (d:g:time)';

            data = (A*pulstran(t,Delay,@rectpuls,W))';
            desired_displacement = data;

            Pulses_per_interval_total = num_pulses(model);
            Freq_test_average = [];

            %Power Pectrum of Desired Displacement
            [Pxx1,f1] = pwelch(desired_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);

            desired_displacement = desired_displacement';

        end

        %%
        %Create Frequency and Amplitude Parameters based on Analog Signal
        Amplitude = desired_displacement*100;  %mV
        Frequency = desired_displacement*14000;  %Hz

        %%
        %Generate Neural Command Signal with Frequency and Amplitdue Parameters
        neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

    %     figure(figNum)
    %     figNum = figNum+1;
    %     plot(t_total,neural)
    %     ax = gca;
    %     ax.FontSize = 26;
    %     xlabel('Time (s)','Fontsize',32)
    %     ylabel('Amplitude (% of MUs)','Fontsize',32);
    %     title('Neural Command Input','Fontsize',36)
    %     grid on

        neural_simulink = [t_total' neural]; %Input for the Simulink Model

        %Power Pectrum of Neural Input
        [Pxx2,f2] = pwelch(neural,gausswin(Nfft),Nfft/2,Nfft,Fs);


        %%
        %Execute the Simulink Model with Set Output Noise

        %Set Output Noise
        set_param('EMG_Model_Simulink/Output Noise','Cov','set_output_noise_power')
        output_noise_power = [output_noise_power set_output_noise_power];

        %Run Simulink;
        out = sim('EMG_Model_Simulink',time);

        %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
        %set_output_noise = set_output_noise_power;

        %%
        %Get Output Signals from Simulink
        emg_simulink = out.EMGout;

        %PDF of EMG
    %     figure(figNum)
    %     figNum = figNum+1;
    %     emg_pdf = emg_simulink;
    %     emg_pdf(emg_pdf==0) = []; %Removing the zero values
    %     histogram(emg_pdf,'Normalization','pdf')
    %     ax = gca;
    %     ax.FontSize = 15;
    %     xlabel('Volts (V)','Fontsize',18)
    %     ylabel('Density','Fontsize',18)
    %     title('EMG Distribution','Fontsize',24)
    %     grid on

        %Power Pectrum of EMG
        [Pxx_emg,f_emg] = pwelch(emg_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);

        force_simulink = out.EMG_Model_Force;

        %Power Pectrum of FOrce
        [Pxx_force,f_force] = pwelch(force_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);

        output_displacement_simulink = out.EMG_Model_Displacement;
        t_simulink = out.tout;

        Zcur = [emg_simulink,output_displacement_simulink];

        Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});

    %     figure(figNum)
    %     figNum = figNum+1;
    %     subplot(2,1,1)
    %     plot(t_total,emg_simulink)
    %     ax = gca;
    %     ax.FontSize = 15;
    %     xlabel('Time (s)','Fontsize',18)
    %     ylabel('EMG (V)','Fontsize',18);
    %     title('(a) Output EMG','Fontsize',24)
    %     grid on
    %     
    %     subplot(2,1,2)
    %     plot(t_total,output_displacement_simulink)
    %     ax = gca;
    %     ax.FontSize = 15;
    %     xlabel('Time (s)','Fontsize',18)
    %     ylabel('Displacement (m)','Fontsize',18);
    %     title('(b) Output Displacement','Fontsize',24)
    %     grid on

        %%
        %Calculate Signal to Noise Ratio
        output_noise_simulink = out.Output_Noise;

        signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
        noise_snr = [noise_snr signal_to_noise];

        %%
        %Hammerstein System Identification (EMG-Movement Model)
        %HK Method
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

        NHK=nlbl;
        set(NHK,'idMethod','hk','displayFlag',true,'threshNSE',.001);

        % Set umber of lags in IRF
        I=NHK{1,2};
        set(I,'nLags',2400); %(accuracy increases if nlags is increased)
        NHK{1,2}=I;

        NHK=nlident(NHK,Zcur);

        figure(model);
        [R, V, yp] = nlid_resid(NHK,Zcur);

        identification_accuracy(trial,model) = V;
        
        if use_fr == true
            validation_frequency(trial,model) = Freq_test_average;
            validation_pulses_total(trial,model) = Pulses_per_interval_total;
        else
            validation_pulses_total(trial,model) = Pulses_per_interval_total;
        end
        
        NHK_all_temp = [NHK_all_temp NHK];
        NHK = [];
        Zcur_all_temp = [Zcur_all_temp Zcur];

    end
    %%
    NHK_all = [NHK_all;NHK_all_temp];
    NHK_all_temp = [];
    Zcur_all = [Zcur_all;Zcur_all_temp];
    Zcur_all_temp = [];
    
end

%% Model Validation
%Set Initial Parmeters
noise_snr = [];
set_output_noise_power = 0;
%noise_multiplier = 10;
output_noise_power = [];

%For Power Spectrums
Fs = 1000; 
Nfft = 1000;

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_movement_max_amplitude = 0.01;
fr = 0.1;                 %Frequency distribution mean (Hz) (Max is 1.8 Hz)
sig = 0.6;                %Std of Frequency Distribution (Hz)
W = 0.55;                  
nf = 18;                  %number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (s)
chance_of_zero = false;


%%
for trial = 1:num_trials_val
    
    %(Physiologically Based Movment - Square pulses with random amplitude and random frequency)
    t_total = 0:0.001:physiological_movement_time;             %Total time
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
        t  = 0 : 0.001 : t_interval;         % Time Samples

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
    [Pxx1,f1] = pwelch(desired_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);

    desired_displacement = desired_displacement';
    
    %%
    %Create Frequency and Amplitude Parameters based on Analog Signal
    Amplitude = desired_displacement*100;  %mV
    Frequency = desired_displacement*14000;  %Hz

    %%
    %Generate Neural Command Signal with Frequency and Amplitdue Parameters
    neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

    neural_simulink = [t_total' neural]; %Input for the Simulink Model

    %Power Pectrum of Neural Input
    [Pxx2,f2] = pwelch(neural,gausswin(Nfft),Nfft/2,Nfft,Fs);
    
    %%
    %Execute Simulink Model with Set Output Noise

    %Set Output Noise
    set_param('EMG_Model_Simulink/Output Noise','Cov','set_output_noise_power')
    output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink Model
    out = sim('EMG_Model_Simulink',time);

    %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
    %set_output_noise = set_output_noise_power;
    
    %%
    %Get Output Signals from Simulink
    emg_simulink = out.EMGout;
    
    emg_pdf = emg_simulink;
    emg_pdf(emg_pdf==0) = []; %Removing the zero values
    
    %Power Pectrum of EMG
    [Pxx_emg,f_emg] = pwelch(emg_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);
    
    force_simulink = out.EMG_Model_Force;
    
    output_displacement_simulink = out.EMG_Model_Displacement;
    t_simulink = out.tout;

    Zcur = [emg_simulink,output_displacement_simulink];

    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});

    %%
    %Calculate Signal to Noise Ratio
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    noise_snr = [noise_snr signal_to_noise];
    
    %%
    %Model Validation
    
    set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
    
    for model = 1:num_models_for_ident
        
        figure(model)
        [R, V, yp] = nlid_resid(NHK_all(1,model),Zcur);
        
        validation_accuracy(trial,model) = V;
        
    end

end
%% Plot Accuracy (%VAF) vs Average Frequency of identification signal

figNum = 100;

if use_fr == true
    figure(figNum)
    figNum = figNum+1;
    plot(fr,validation_accuracy,'LineWidth', 3)
    ax = gca;
    ax.FontSize = 26;
    title('Expected Average Frequency vs Validation Accuracy','Fontsize', 36)
    xlabel('Frequency (Hz)','Fontsize', 32)
    ylabel('Accuracy (%VAF)','Fontsize', 32)
    grid on

    figure(figNum)
    figNum = figNum+1;
    plot(validation_frequency,validation_accuracy,'LineWidth', 3)
    ax = gca;
    ax.FontSize = 26;
    title('Observed Average Frequency vs Validation Accuracy','Fontsize', 36)
    xlabel('Frequency (Hz)','Fontsize', 32)
    ylabel('Accuracy (%VAF)','Fontsize', 32)
    grid on
    
else
    
    
    identification_accuracy = max(0, identification_accuracy);
    identification_accuracy_mean = mean(identification_accuracy);
    identification_accuracy_var = var(identification_accuracy);
    identification_accuracy_std = std(identification_accuracy);
    
    figure(figNum)
    figNum = figNum+1;
    hold on
    plot(validation_pulses_total(1,:),min(identification_accuracy_mean+identification_accuracy_std,100),'LineStyle','none','LineWidth',2)
    plot(validation_pulses_total(1,:),max(identification_accuracy_mean-identification_accuracy_std,0),'LineStyle','none','LineWidth',2)
    patch([validation_pulses_total(1,:) fliplr(validation_pulses_total(1,:))], [min(identification_accuracy_mean+identification_accuracy_std,100) fliplr(max(identification_accuracy_mean-identification_accuracy_std,0))], 'b','FaceAlpha',0.2)
    plot(validation_pulses_total(1,:),identification_accuracy_mean,'LineWidth', 3)
    hold off
    ax = gca;
    ax.FontSize = 14;
    title('Total Pulses vs Identifcation Accuracy','Fontsize', 24) 
    xlabel('Number of Signal Pulses','Fontsize', 18)
    ylabel('Accuracy (%VAF)','Fontsize', 18)
    grid on
    
    validation_accuracy = max(0, validation_accuracy);
    validation_accuracy_mean = mean(validation_accuracy);
    validation_accuracy_var = var(validation_accuracy);
    validation_accuracy_std = std(validation_accuracy);
    
    figure(figNum)
    figNum = figNum+1;
    hold on
    % errorbar(validation_pulses_total,validation_accuracy_mean,validation_accuracy_std,'CapSize',14,'LineWidth',3)
    % plot(validation_pulses_total,min(validation_accuracy_mean+validation_accuracy_std,100),'--r','LineWidth',2)
    % plot(validation_pulses_total,max(validation_accuracy_mean-validation_accuracy_std,0),'--r','LineWidth',2)
    plot(validation_pulses_total(1,:),min(validation_accuracy_mean+validation_accuracy_std,100),'LineStyle','none','LineWidth',2)
    plot(validation_pulses_total(1,:),max(validation_accuracy_mean-validation_accuracy_std,0),'LineStyle','none','LineWidth',2)
    patch([validation_pulses_total(1,:) fliplr(validation_pulses_total(1,:))], [min(validation_accuracy_mean+validation_accuracy_std,100) fliplr(max(validation_accuracy_mean-validation_accuracy_std,0))], 'b','FaceAlpha',0.2)
    plot(validation_pulses_total(1,:),validation_accuracy_mean,'LineWidth',3)
    hold off
    ax = gca;
    ax.FontSize = 14;
    title('Total Pulses vs Validation Accuracy','Fontsize', 24) 
    xlabel('Number of Signal Pulses','Fontsize', 18)
    ylabel('Accuracy (%VAF)','Fontsize', 18)
    grid on
    
end

tEnd = toc(tStart)/60


