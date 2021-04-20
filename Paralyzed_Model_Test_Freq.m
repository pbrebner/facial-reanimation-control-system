%Paralyzed Model Test

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
%set_output_noise_power = 1e-10;        %Output Noise Power
set_output_noise_power = 0;
% time = 180;                             %Total time in seconds
%noise_multiplier = 10;
% accuracy = [];
% record_length = [];
noise_snr = [];
output_noise_power = [];

%For Power Spectrums
Fs = 1000; 
Nfft = 1000;

PRBS_stimulus = false;
% physiological_stimulus = true;

%PRBS Stimulus
% PRBS_stimulus_time = 180;
% variable_amplitude = false;
% N = PRBS_stimulus_time/10;
% M = 1000;
% PRBS_amplitude = 10; %mm

%Set Physiological Signal Parameters
physiological_stimulus_time = 180;
physiological_stimulus_max_amplitude = 0.01;
fr = 0.1;                 %Frequency distribution mean (Hz)
sig = 0.8;                %Std of Frequency Distribution (Hz)
W = 0.45;                  
nf = physiological_stimulus_time/10;                  %number of random signal changes
t_interval = physiological_stimulus_time/nf;    %Length of random interval (s)
chance_of_zero = false;

LNL_model = false;
Hammerstein_model = false;
Weiner_model = true;
Linear_IRF_model = true;

LNL_all_temp = [];
LNL_all = [];
Hammerstein_all_temp = [];
Hammerstein_all = [];
Weiner_all_temp = [];
Weiner_all = [];
IRF_model_all_temp = [];
IRF_model_all = [];
NLN_all_temp = [];
NLN_all = [];
Zcur_all_temp = [];
Zcur_all = [];

identification_accuracy = [];
validation_accuracy = [];

%%
use_fr = false; %Use frequency instead of number of pulses for Physiological Input
if use_fr == true
    fr = 0:0.1:0.3;
    sig = 0.6;
end

num_pulses = 1:1:30;
num_trials_ident = 1;
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
            t_total = 0:0.001:((physiological_stimulus_time/num_models_for_ident)*model);
            time = ((physiological_stimulus_time/num_models_for_ident)*model);
        else
            t_total = 0:0.001:physiological_stimulus_time;
            time = physiological_stimulus_time;
        end

        if use_fr == true
            FR = makedist('Normal','mu',fr(model),'sigma',sig);
            FrequenciesRandom_max = 2.1;
            FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
            freq_distribution = random(FrequenciesRandom,10000,1);
        %     figure(figNum)
        %     figNum = figNum+1;
        %     histogram(freq_distribution,100)
        %     title('Frequency Distribution')
        %     xlabel('Frequency (Hz)')
        end

        AR = makedist('Uniform','lower',0,'upper',physiological_stimulus_max_amplitude);  %Full Amplitude Range
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
            A = physiological_stimulus_max_amplitude;
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
    %                         figure(figNum)
    %                         figNum = figNum+1;
    %                         histogram(delay_random_distribution,100)

                random_delay = random(DelaysRandom,1,1);

                Delay = [Delay hh*g+random_delay];
                %Delay = [Delay Delay(end)+random_delay];

            end

            %Delay = (d:g:time)';

            data = (A*pulstran(t,Delay,@rectpuls,W))';
            desired_displacement = data;

            Pulses_per_interval_total = num_pulses(model);
            Freq_test_average = [];

        end


    %     figure(figNum)
    %     figNum = figNum+1;
    %     plot(t_total,desired_displacement)
    %     ax = gca;
    %     ax.FontSize = 15;
    %     xlabel('Time (s)','Fontsize',18)
    %     ylabel('Desired Displacement (m)','Fontsize',18)
    %     title('Analog Signal of Desired Displacement','Fontsize',24)
    %     grid on

        %FFT of Desired Displacement
        L = length(desired_displacement);

        Y = fft(desired_displacement);
        P1 = abs(Y/L);
        Pxx1 = P1(1:L/2+1);
        Pxx1(2:end-1) = 2*Pxx1(2:end-1);
    %     Pxx1 = Pxx1';

        f1 = (Fs*(0:(L/2))/L)';
    %     figure(figNum)
    %     figNum = figNum+1;
    %     plot(f1(1:5000,:),Pxx1(1:5000,:)) 
    %     title('FFT of Desired Displacement')
    %     xlabel('f (Hz)')
    %     ylabel('|Pxx1(f)|')

        %Power Pectrum of Desired Displacement
        %[Pxx1,f1] = pwelch(desired_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);

        % figure(figNum)
        % figNum = figNum+1;
        % subplot(2,1,1)
        % plot(t_total,desired_displacement)
        % ax = gca;
        % ax.FontSize = 14;
        % xlabel('Time (s)','Fontsize',18)
        % ylabel('Desired Displacement (m)','Fontsize',18)
        % title('Analog Signal of Desired Movement','Fontsize',24)
        % grid on

        %Plot frequency spectrum
        % subplot(2,1,2)
        % semilogy(f1(1:6,1),Pxx1(1:6,1));
        % ax = gca;
        % ax.FontSize = 14;
        % title('Power Spectrum of Desired Displacement','Fontsize',24);
        % ylabel('PSD (log)','Fontsize',18); 
        % xlabel('Frequency (Hz)','Fontsize',18);
        % grid on;

        desired_displacement = desired_displacement';

        stim_amplitude = desired_displacement*170;  %mV
        stim_frequency = 50;

    %     input_stimulus = max(stim_amplitude.*square(2*pi*stim_frequency.*t_total),0);
        input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);

        amplitude_modulation = stim_amplitude;
        amplitude_modulation_simulink = [t_total' amplitude_modulation'];


        %FFT of Input Stimulus
        L = length(input_stimulus);

        Y = fft(input_stimulus);
        P2 = abs(Y/L);
        Pxx2 = P2(1:L/2+1);
        Pxx2(2:end-1) = 2*Pxx2(2:end-1);
        Pxx2 = Pxx2';

        f2 = (Fs*(0:(L/2))/L)';
    %     figure(figNum)
    %     figNum = figNum +1;
    %     plot(f2(1:15000,:),Pxx2(1:15000,:)) 
    %     title('FFT of Input Stimulus')
    %     xlabel('f (Hz)')
    %     ylabel('|Pxx2(f)|')

        %Power Pectrum of Stimulus Input
        %[Pxx2,f2] = pwelch(input_stimulus,gausswin(Nfft),Nfft/2,Nfft,Fs);

        stimulus_simulink = [t_total' input_stimulus'];


        %%
        %Execute Simulink Model with Set Output Noise

        set_param('Paralyzed_Model_Simulink/Output Noise','Cov','set_output_noise_power')
        output_noise_power = [output_noise_power set_output_noise_power];

        %Run Simulink;
        out = sim('Paralyzed_Model_Simulink',time);

        %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
        %set_output_noise = set_output_noise_power;

        %%
        %Get Output Signals from Simulink

        %Muscle Force
        force_simulink = out.Paralyzed_Model_Force;

        % figure(figNum)
        % figNum = figNum+1;
        % plot(t_total,force_simulink)
        % ax = gca;
        % ax.FontSize = 15;
        % xlabel('Time (s)','Fontsize',18)
        % ylabel('Force (N)','Fontsize',18);
        % title('Muscle Force','Fontsize',24)
        % grid on

        %Power Pectrum of Force
        [Pxx_force,f_force] = pwelch(force_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);


        %%
        %Input Stimulus and Output Displacement
        input_stimulus = out.Paralyzed_Model_Stimulus;
        output_displacement_simulink = out.Paralyzed_Model_Displacement;
        t_simulink = out.tout;

        Zcur = [amplitude_modulation',output_displacement_simulink];

        Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});

        %%
        %Calculate Signal to Noise Ratio
        output_noise_simulink = out.Output_Noise;

        signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
        noise_snr = [noise_snr signal_to_noise];

        %%
        %Model Identification (LNL Model or NLN Model)

        if LNL_model == true

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

            LNL=lnlbl;
            set(LNL,'idMethod','hk','nLags1',750,'polyOrderMax',2,'nLags2',750,... 
            'hkTolerance', 0.1, 'nhkMaxIts', 10, 'nhkMaxInner', 5);

            LNL=nlident(LNL,Zcur);

            figure(model)
            [R, V, yp] = nlid_resid(LNL,Zcur);

            identification_accuracy(trial,model) = V;

            if use_fr == true
                validation_frequency(trial,model) = Freq_test_average;
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            else
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            end

            LNL_all_temp = [LNL_all_temp LNL];
            LNL = [];
            Zcur_all_temp = [Zcur_all_temp Zcur];
        
        elseif Hammerstein_model == true
        
            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

            Hammerstein=nlbl;
            set(Hammerstein,'idMethod','hk','displayFlag',true,'threshNSE',.001);
            I=Hammerstein{1,2};
            set(I,'nLags',2400,'nSides',1); %(accuracy increases if nlags is increased)
            Hammerstein{1,2}=I;

            Hammerstein=nlident(Hammerstein,Zcur);

            figure(model)
            [R, V, yp] = nlid_resid(Hammerstein,Zcur);

            identification_accuracy(trial,model) = V;

            if use_fr == true
                validation_frequency(trial,model) = Freq_test_average;
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            else
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            end

            Hammerstein_all_temp = [Hammerstein_all_temp Hammerstein];
            Hammerstein = [];
            Zcur_all_temp = [Zcur_all_temp Zcur];

        elseif Weiner_model == true

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

            Weiner = lnbl; %Wiener
            set(Weiner,'idMethod','hk');
            I = Weiner{1,1};
            set(I,'nLags',850,'nSides',1); % Set Number of lags and Sides in IRF
            Weiner{1,1} = I;

            Weiner=nlident(Weiner,Zcur);

            figure(model)
            [R, V, yp] = nlid_resid(Weiner,Zcur);

            identification_accuracy(trial,model) = V;

            if use_fr == true
                validation_frequency(trial,model) = Freq_test_average;
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            else
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            end

            Weiner_all_temp = [Weiner_all_temp Weiner];
            Weiner = [];
            Zcur_all_temp = [Zcur_all_temp Zcur];

        elseif Linear_IRF_model == true
            % Identify a two-sided IRF Model

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

            IRF_model = irf(Zcur,'nLags',1200,'nSides',1);

            figure(model)
            [R, V, yp] = nlid_resid(IRF_model,Zcur);

            identification_accuracy(trial,model) = V;

            if use_fr == true
                validation_frequency(trial,model) = Freq_test_average;
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            else
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            end

            IRF_model_all_temp = [IRF_model_all_temp IRF_model];
            IRF_model = [];
            Zcur_all_temp = [Zcur_all_temp Zcur];
            
        else
            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
            Zcur_NLN = iddata(output_displacement_simulink,input_stimulus',Fs);
            %'domainIncr',0.001,'comment','Input Stimulus, Output Displacement','chanNames', {'Stimulus (V)' 'Displacement (m)'});

            set(Zcur_NLN,'InputName','Input Stimulus','OutputName','Output Displacement');

            NLN = nlhw(Zcur_NLN,[2 2 2]);

            [y,fit,ic] = compare(Zcur_NLN,NLN);

            pred = y.y(:,1);
            V = vaf(output_displacement_simulink,pred);
            identification_accuracy(trial,model) = V; 

            R = output_displacement_simulink - pred;

            %Plots just the superimposed part with %VAF in the title
            figure(model);
            plot(t_total,pred);
            hold on
            plot(t_total, output_displacement_simulink)
            ax = gca;
            ax.FontSize = 18;
            hold off
            title(['Superimposed, VAF = ' num2str(V) '%'], 'Fontsize', 28)
            xlabel('Time (s)', 'Fontsize', 24)
            ylabel('Displacement (m)', 'Fontsize', 24)
            legend('Predicted', 'Observed', 'Fontsize', 20)
            grid on

            if use_fr == true
                validation_frequency(trial,model) = Freq_test_average;
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            else
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            end

            NLN_all_temp = [NLN_all_temp NLN];
            NLN = [];
            Zcur_all_temp = [Zcur_all_temp Zcur];

        end

    end

    LNL_all = [LNL_all;LNL_all_temp];
    LNL_all_temp = [];
    
    Hammerstein_all = [Hammerstein_all;Hammerstein_all_temp];
    Hammerstein_all_temp = [];
    
    Weiner_all = [Weiner_all;Weiner_all_temp];
    Weiner_all_temp = [];
    
    IRF_model_all = [IRF_model_all;IRF_model_all_temp];
    IRF_model_all_temp = [];

    NLN_all = [NLN_all;NLN_all_temp];
    NLN_all_temp = [];

    Zcur_all = [Zcur_all;Zcur_all_temp];
    Zcur_all_temp = [];

end


%% Model Validation
%Set initial Parameters
%set_output_noise_power = 1e-10;        %Output Noise Power
set_output_noise_power = 0;
%noise_multiplier = 10;
noise_snr = [];
output_noise_power = [];

%For Power Spectrums
Fs = 1000; 
Nfft = 1000;

%Set Physiological Signal Parameters
physiological_stimulus_time = 180;
physiological_stimulus_max_amplitude = 0.010;
fr = 0.1;                 %Frequency distribution mean (Hz)
sig = 0.8;                %Std of Frequency Distribution (Hz)
W = 0.45;                  
nf = physiological_stimulus_time/10;                  %number of random signal changes
t_interval = physiological_stimulus_time/nf;    %Length of random interval (s)
chance_of_zero = false;


%%
for trial = 1:num_trials_val

    t_total = 0:0.001:physiological_stimulus_time;
    time = physiological_stimulus_time;

    FR = makedist('Normal','mu',fr,'sigma',sig);
    FrequenciesRandom_max = 2.21;
    FrequenciesRandom = truncate(FR,0,1.8);
    freq_distribution = random(FrequenciesRandom,10000,1);

    AR = makedist('Uniform','lower',0,'upper',physiological_stimulus_max_amplitude);  %Full Amplitude Range
    AmplitudesRandom = AR;
    amp_distribution = random(AmplitudesRandom,10000,1);

    desired_displacement= 0;
    Freq_test = [];
    Pulses_per_interval_test = [];

    for j = 1 : nf    
        t  = 0 : 0.001 : t_interval;         % Time Samples

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

        %data = A*square(2*pi*f/Fs*t)';        % Generate Square Wave
        desired_displacement = [desired_displacement; data];
        Freq_test = [Freq_test stim_frequency];
        Pulses_per_interval_test = [Pulses_per_interval_test Pulses_per_interval];
    end

    Pulses_per_interval_total = sum(Pulses_per_interval_test);
    Freq_test_average = sum(Freq_test)/length(Freq_test);
    
    %FFT of Desired Displacement
    L = length(desired_displacement);
    
    Y = fft(desired_displacement);
    P1 = abs(Y/L);
    Pxx1 = P1(1:L/2+1);
    Pxx1(2:end-1) = 2*Pxx1(2:end-1);
%     Pxx1 = Pxx1';
    
    f1 = (Fs*(0:(L/2))/L)';
%     figure(figNum)
%     figNum = figNum+1;
%     plot(f1(1:5000,:),Pxx1(1:5000,:)) 
%     title('FFT of Desired Displacement')
%     xlabel('f (Hz)')
%     ylabel('|Pxx1(f)|')
    
    %Power Pectrum of Desired Displacement
    %[Pxx1,f1] = pwelch(desired_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);

    desired_displacement = desired_displacement';

    stim_amplitude = desired_displacement*170;  %mV
    stim_frequency = 50;

    %input_stimulus = max(stim_amplitude.*square(2*pi*stim_frequency.*t_total),0);
    input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);
    
    amplitude_modulation = stim_amplitude;
    amplitude_modulation_simulink = [t_total' amplitude_modulation'];

    %FFT of Input Stimulus
    L = length(input_stimulus);
    
    Y = fft(input_stimulus);
    P2 = abs(Y/L);
    Pxx2 = P2(1:L/2+1);
    Pxx2(2:end-1) = 2*Pxx2(2:end-1);
    Pxx2 = Pxx2';
    
    f2 = (Fs*(0:(L/2))/L)';
%     figure(figNum)
%     figNum = figNum +1;
%     plot(f2(1:15000,:),Pxx2(1:15000,:)) 
%     title('FFT of Input Stimulus')
%     xlabel('f (Hz)')
%     ylabel('|Pxx2(f)|')
    
    %Power Pectrum of Stimulus Input
    %[Pxx2,f2] = pwelch(input_stimulus,gausswin(Nfft),Nfft/2,Nfft,Fs);

    stimulus_simulink = [t_total' input_stimulus'];
    
    %%
    %Execute Simulink Model with Set Output Noise

    set_param('Paralyzed_Model_Simulink/Output Noise','Cov','set_output_noise_power')
    output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink;
    out = sim('Paralyzed_Model_Simulink',time);

    %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
    %set_output_noise = set_output_noise_power;

    %%
    %Get Output Signals from Simulink
    %Muscle Force
    force_simulink = out.Paralyzed_Model_Force;

    % figure(figNum)
    % figNum = figNum+1;
    % plot(t_total,force_simulink)
    % ax = gca;
    % ax.FontSize = 15;
    % xlabel('Time (s)','Fontsize',18)
    % ylabel('Force (N)','Fontsize',18);
    % title('Muscle Force','Fontsize',24)
    % grid on

    %Power Pectrum of Force
    [Pxx_force,f_force] = pwelch(force_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);


    %%
    %Input Stimulus and Output Displacement
    input_stimulus = out.Paralyzed_Model_Stimulus;
    output_displacement_simulink = out.Paralyzed_Model_Displacement;
    t_simulink = out.tout;

    %Zcur = [input_stimulus',output_displacement_simulink];
    Zcur = [amplitude_modulation',output_displacement_simulink];

    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});

    
    %%
    %Calculate Signal to Noise Ratio
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    noise_snr = [noise_snr signal_to_noise];

    
    %%
    %Model Validation

    if LNL_model == true

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for model = 1:num_models_for_ident
            
            figure(model)
            [R, V, yp] = nlid_resid(LNL_all(1,model),Zcur);

            validation_accuracy(trial,model) = V;
        end
        
    elseif Hammerstein_model == true
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for model = 1:num_models_for_ident
            
            figure(model)
            [R, V, yp] = nlid_resid(Hammerstein_all(1,model),Zcur);

            validation_accuracy(trial,model) = V;
        end
        
    elseif Weiner_model == true
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for model = 1:num_models_for_ident
            
            figure(model)
            [R, V, yp] = nlid_resid(Weiner_all(1,model),Zcur);

            validation_accuracy(trial,model) = V;
        end
        
    elseif Linear_IRF_model == true
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for model = 1:num_models_for_ident
            
            figure(model)
            [R, V, yp] = nlid_resid(IRF_model_all(1,model),Zcur);

            validation_accuracy(trial,model) = V;
        end

    else

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

        Zcur_NLN = iddata(output_displacement_simulink,input_stimulus',Fs);
        set(Zcur_NLN,'InputName','Input Stimulus','OutputName','Output Displacement');
        
        for model = 1:num_models_for_ident

            [y,fit,ic] = compare(Zcur_NLN,NLN_all(1,model));

            pred = y.y(:,1);
            V = vaf(output_displacement_simulink,pred);
            validation_accuracy(trial,model) = V; 

            R = output_displacement_simulink - pred;

            %Plots just the superimposed part with %VAF in the title
            figure(model);
            plot(t_total,pred);
            hold on
            plot(t_total, output_displacement_simulink)
            ax = gca;
            ax.FontSize = 18;
            hold off
            title(['Superimposed, VAF = ' num2str(V) '%'], 'Fontsize', 28)
            xlabel('Time (s)', 'Fontsize', 24)
            ylabel('Displacement (m)', 'Fontsize', 24)
            legend('Predicted', 'Observed', 'Fontsize', 20)
            grid on

        end

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
    xlabel('Number of Movement Pulses in Identification Signal','Fontsize', 18)
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
    xlabel('Number of Movement Pulses in Identification Signal','Fontsize', 18)
    ylabel('Accuracy (%VAF)','Fontsize', 18)
    grid on
    
end

tEnd = toc(tStart)/60

