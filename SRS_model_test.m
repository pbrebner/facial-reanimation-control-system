%% SRS Model Test (Model Validation)

%Validate the estimated model by using it to predict the movements
%generated in response to an input realization which is different than that
%used to identify the model. Will validate the specified model 
%with a PRBS input signal and then a "Physiological" input signal. Plots the 
%superimposed predicted and observed displacements, residual distribution, 
%and residual spectrum 

%INSTRUCTIONS: RUN Paralyzed_Model_Updated FIRST and THEN RUN THIS

%When running the script, you need to provide the following input:
%If both PRBS and Physiological Models were Identified in SRS_model,
% 1. Which Identified Model Type do you want to Validate? PRBS/Phys
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
        SRS_model = SRS_models(1);
        disp(['Validating ' model_type ' Model Identified from PRBS Input'])
    elseif strcmp(str1,'Phys')
        SRS_model = SRS_models(2);
        disp(['Validating ' model_type ' Model Identified from Physiological Input'])
    end
else
    SRS_model = SRS_models(1);
    disp(['Validating ' model_type ' Model Identified from ' input_type ' Input'])
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
set_output_noise_power = 0;
noise_snr = [];
output_noise_power = [];
figNum = 100;

%For Power Spectrums
Fs = 1000; 
Nfft = 200000;

%PRBS Signal Parameters
PRBS_stimulus_time = 480;
variable_amplitude = true;     %PRBS can either be constant amplitude or variable amplitude
N = PRBS_stimulus_time/10;     %Number of times the amplitude randomly changes (For Variable Amplitude Only)
M = 10000;                     %Number of each random value (For Variable Amplitude Only)
PRBS_amplitude = 20;           %Amplitude of PRBS Signal (mm)

%Set Physiological Signal Parameters
physiological_stimulus_time = 480;
physiological_stimulus_max_amplitude = 0.02;    %"Physiological" Amplitude (m)
fr = 0.1;                                       %Frequency distribution mean (Hz)
sig = 0.8;                                      %Std of Frequency Distribution (Hz)
W = 0.45;                                       %Width of signal pulse (seconds) 
nf = physiological_stimulus_time/10;            %number of random signal changes
t_interval = physiological_stimulus_time/nf;    %Length of random interval (s)
chance_of_zero = false;

% %Set which model you want to use
% if compare_two_models == true
%     if LNL_model == true
%         LNL = LNL1;
%     elseif Hammerstein_model == true
%         Hammerstein = Hammerstein1;
%     elseif Weiner_model == true
%         Weiner = Weiner2;
%     elseif Linear_IRF_model == true
%         IRF_model = IRF1;
%     end
% end

%Number of Validation Trials
num_trials = str2;
validation_accuracy = [];
validation_accuracy_all = [];

%for validation we test with PRBS signal first, and then with Physiological signal
PRBS_stimulus = [true false];

%% Generate Desired Displacement and Amplitude Modulation Signals for Model Validation
for signal = 1:length(PRBS_stimulus)
    
    for trial = 1:num_trials

        if PRBS_stimulus(signal) == true

            t_total = 0:0.001:PRBS_stimulus_time;
            time = PRBS_stimulus_time;

            A = [0];                    %Intialize amplitude
            if variable_amplitude == true   
                for k = 1:N
                    if k == 1
                        R = PRBS_amplitude;
                    else
                        R = rand(1,1)*PRBS_amplitude;   %Randomly generate a number between 0 and 10
                    end

                    for j = 1:M
                        A = [A R];
                    end
                end
            else
                A = PRBS_amplitude;              %Constant Amplitude
            end

            Range = [0,0.001]; %Specify that the single-channel PRBS value switches between -2 and 2

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

            %FFT of Desired Displacement
            L = length(desired_displacement);

            Y = fft(desired_displacement);
            P1 = abs(Y/L);
            Pxx1 = P1(1:L/2+1);
            Pxx1(2:end-1) = 2*Pxx1(2:end-1);
            Pxx1 = Pxx1';
            f1 = (Fs*(0:(L/2))/L)';

            stim_frequency = 50;
            stim_amplitude = desired_displacement*170;

            input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);

            amplitude_modulation = stim_amplitude;
            amplitude_modulation_simulink = [t_total' amplitude_modulation'];

        else

            t_total = 0:0.001:physiological_stimulus_time;
            time = physiological_stimulus_time;

            FR = makedist('Normal','mu',fr,'sigma',sig);
            FrequenciesRandom_max = 2.21;
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
            f1 = (Fs*(0:(L/2))/L)';

            desired_displacement = desired_displacement';

            stim_amplitude = desired_displacement*170;  %mV
            stim_frequency = 50;

            input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);

            amplitude_modulation = stim_amplitude;
            amplitude_modulation_simulink = [t_total' amplitude_modulation'];

        end

        %% Execute SRS Simulation (Simulink Model)

        %Set Output Noise as zero
        set_param('SRS_simulation/Output Noise','Cov','set_output_noise_power')
        output_noise_power = [output_noise_power set_output_noise_power];

        %Run Simulink;
        out = sim('SRS_simulation',time);

        %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
        %set_output_noise = set_output_noise_power;

        %% Get Output Signals from SRS Simulation (Simulink Model)

        %Muscle Force
        force_simulink = out.Paralyzed_Model_Force;

        %Electrical Stimulus
        input_stimulus = out.Paralyzed_Model_Stimulus;

        %FFT of Electrical Stimulus
        L = length(input_stimulus);

        Y = fft(input_stimulus);
        P2 = abs(Y/L);
        Pxx2 = P2(1:L/2+1);
        Pxx2(2:end-1) = 2*Pxx2(2:end-1);
        f2 = (Fs*(0:(L/2))/L)';

        %Paralyzed Displacement Output
        output_displacement_simulink = out.Paralyzed_Model_Displacement;
        t_simulink = out.tout;

        %% Input/Output for Model Validation

        Zcur = [amplitude_modulation',output_displacement_simulink];
        Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});

        %% Calculate Signal to Noise Ratio

        %Get noise from SRS Simulation
        output_noise_simulink = out.Output_Noise;

        signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
        noise_snr = [noise_snr signal_to_noise];

        %% Model Validation (LNL, Hammerstein, Wiener, or Linear IRF)

        if LNL_model == true

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

            figure(figNum)
            [R, V, yp] = nlid_resid(SRS_model,Zcur);

            validation_accuracy = [validation_accuracy V];

        elseif Hammerstein_model == true

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

            figure(figNum)
            [R, V, yp] = nlid_resid(SRS_model,Zcur);

            validation_accuracy = [validation_accuracy V];

        elseif Weiner_model == true

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

            figure(figNum)
            [R, V, yp] = nlid_resid(SRS_model,Zcur);

            validation_accuracy = [validation_accuracy V];    

        elseif Linear_IRF_model == true

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

            figure(figNum)
            [R, V, yp] = nlid_resid(SRS_model,Zcur);

            validation_accuracy = [validation_accuracy V];   

        end

    end
    
    %Store the validation accuracy for the signal type and reset for next
    %signal iteration
    validation_accuracy_all(signal,:) = validation_accuracy;
    validation_accuracy = [];
    
    %% Plot the Validation results of last validation trial
    
    pred = double(yp);

    %Plot Sumperimposed Accuracy, Residual PDF and Spectrum
    figure(figNum)
    figNum = figNum+1;
    subplot(2,2,[1 2])
    plot(t_total,pred);
    hold on
    plot(t_total, output_displacement_simulink)
    ax = gca;
    ax.FontSize = 15;
    hold off
    title(['Superimposed, VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 20)
    xlabel('Time (s)', 'Fontsize', 16)
    ylabel('Displacement, Pos_P(t) (m)', 'Fontsize', 16)
    legend('Predicted', 'Observed', 'Fontsize', 15)

    subplot(2,2,3)
    histogram(double(R))
    ax = gca;
    ax.FontSize = 15;
    xlabel('Displacement (m)','Fontsize',16)
    ylabel('Density','Fontsize',16)
    title('(b) Residual Distribution','Fontsize',20)
    grid on

    R_zero = R - mean(R);

    [PxxR,fR] = pwelch(double(R_zero),[],[],[],Fs);
    subplot(2,2,4)
    semilogy(fR(1:1300),PxxR(1:1300,:),'LineWidth',1);
    ax = gca;
    ax.FontSize = 15;
    title('(c) Residual Power Spectrum','Fontsize',20);
    ylabel('PSD (log)','Fontsize',16); 
    xlabel('Frequency (Hz)','Fontsize',16);
    grid on
    
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