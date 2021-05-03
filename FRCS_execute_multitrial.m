%% Execute Control System for x number of trials
% Does not display figures 

%Run after identifying the EMG Response System (ERS) and inverse Stimulus
%Response System (SRS)^1

%Generates an input desired displacement, the corresponding nerual command
%and runs it through the Hammerstein Model followed by the Inverse of the
%LNL Model to get the stimulus

tStart = tic;

%% Set initial Parameters

time = 180;
figNum = 300;

Fs = 1000;
Nfft = 1000;

EMG_response_model = NHK2;
% stimulus_response_model = SRS_inverse;

PRBS_movement = false;
% physiological_movement = true;

%PRBS Signal Parameters
PRBS_movement_time = 180;
variable_amplitude = false;
N = PRBS_movement_time/10;
M = 10000;
PRBS_amplitude = 10; %mm

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_stimulus_max_amplitude = 0.01;
fr = 0.1;                                       %Frequency distribution mean (Hz)
sig = 0.6;                                      %Std of Frequency Distribution (Hz)
W = 0.55;
nf = 18;                                        %Number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (s)
chance_of_zero = false;

%Select which Model to Use
if compare_two_models == true && LNL_model == true
    SRS_model = LNL2;
elseif compare_two_models == true && Hammerstein_model == true
    SRS_model = Hammerstein2;
elseif compare_two_models == true && Weiner_model == true
    SRS_model = Weiner2;
elseif compare_two_models == true && Linear_IRF_model == true
    SRS_model = IRF2;
elseif compare_two_models == false && LNL_model == true
    SRS_model = LNL;
elseif compare_two_models == false && Weiner_model == true
    SRS_model = Weiner;
else
    SRS_model = IRF_model;
end

num_trials = 100;
num_models = 3;

%% Simulated Input Test 2 (PRBS Input)
%PRBS Stimulus
simulated_input_PRBS = [];

num_signals = num_models;

PRBS_stimulus_time = 480;
variable_amplitude = true;
N = PRBS_stimulus_time/10;
M = 10000;
PRBS_amplitude = 20; %mm

for signal = 1:num_signals

    t_total = 0:0.001:PRBS_stimulus_time;
    time = PRBS_stimulus_time;

    A = [0];                    %Intialize amplitude
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

    stim_frequency = 50;
    stim_amplitude = desired_displacement*170;

    input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);

    simulated_input_PRBS(:,signal) = stim_amplitude;

end

%% Simulate the reponse of the SRS model
simulated_input = simulated_input_PRBS;

simulated_output = [];
simulated_outputs = [];

for signal = 1:num_signals

    simulated_output = nlsim(SRS_model,simulated_input(:,signal));
    set(simulated_output, 'domainIncr',0.001)

    simulated_output = max(simulated_output,0);

    simulated_outputs(:,signal) = simulated_output.dataSet;

    Zcur_simulated(:,:,signal) = [simulated_outputs(1000:end,signal) simulated_input(1000:end,signal)];
    Zcur_simulated(:,:,signal) = nldat(Zcur_simulated(:,:,signal),'domainIncr',0.001,'comment',...
        'Simulated Output (as Input), Simulated Input (as Output)',...
        'chanNames', {'Simulated Output (as Input)' 'Simulated Input (as Output)'});

end
 
%%

SRS_inverse_all = [];
control_system_displacement_double = [];
control_system_output_double = [];
predicted_healthy_displacement = [];
amplitude_modulation = [];
predicted_paralyzed_displacement = [];
Variance = [];
control_system_validation_accuracy = [];

for model = 1:num_models
    nLags = 400;
    
    %LNL Model
    SRS_inverse = lnlbl;  %LNL Model
    set(SRS_inverse,'idMethod','hk','hkTolerance', 0.5,...
        'nhkMaxIts', 3, 'nhkMaxInner', 3);
    
    I1 = irf;
    set(I1,'nLags',1,'nSides',2,'domainIncr',0.001);
    P = polynom;
    I3 = irf;
    set(I3,'nLags',nLags,'nSides', 2,'domainIncr',0.001);
    SRS_inverse.elements = {I1 P I3};
    
%     SRS_inverse = nlbl;  %Hammerstein
%     set(SRS_inverse,'idMethod','hk','displayFlag',true,'threshNSE',.001);
%     I2 = irf;
%     set(I2,'nLags',nLags, 'nSides', 2); % Set number of lags and Sides in IRF
%     SRS_inverse{1,2} = I2;
    
    SRS_inverse  = nlident(SRS_inverse,Zcur_simulated(:,:,model));
    
    SRS_inverse_all = [SRS_inverse_all SRS_inverse];
    
    figure(figNum)
    figNum = figNum+1;
    [R, V, yp] = nlid_resid(SRS_inverse_all(model),Zcur_simulated(:,:,model));
    
    figure(figNum)
    figNum = figNum+1;
    plot(SRS_inverse_all(model))

    %%
    
    for trial = 1:num_trials

        %Generate the Desired Displacement Signal
        if PRBS_movement == true

            t_total = 0:0.001:PRBS_movement_time;
            time = PRBS_movement_time;

            A = 0;                    %Intialize amplitude
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

                desired_displacement = [desired_displacement; data];
                Freq_test = [Freq_test stim_frequency];
                Pulses_per_interval_test = [Pulses_per_interval_test Pulses_per_interval];
            end

            Pulses_per_interval_total = sum(Pulses_per_interval_test);
            Freq_test_average = sum(Freq_test)/length(Freq_test);

            desired_displacement = desired_displacement';

        end

        %% Generate Neural Input Command Signal
        Amplitude = desired_displacement*100;  %mV
        Frequency = desired_displacement*14000;  %Hz

        neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

        neural_simulink = [t_total' neural];

        %% Generate EMG Signal (Simulink)
        
        %Run Simulink;
        out = sim('EMG_Generation_Simulink',time);

        %Get Output from Simulink
        EMG_simulink = out.EMGout;
        t_simulink = out.tout;

        EMG_simulink = nldat(EMG_simulink);
        set(EMG_simulink, 'domainIncr',0.001);

        %% Simulate Response and Plot Results

        %Set domain increments of EMG response model
        EMG_response_model_IRF = EMG_response_model{1,2};
        set(EMG_response_model_IRF, 'domainIncr', 1.0e-3);
        EMG_response_model{1,2} = EMG_response_model_IRF;

        control_system_displacement = nlsim(EMG_response_model,EMG_simulink);
        set(control_system_displacement,'domainIncr',0.001);
        %control_system_displacement = control_system_displacement.dataSet;

        control_system_output = nlsim(SRS_inverse_all(:,model),control_system_displacement);
        %control_system_stimulus_output = control_system_stimulus_output.dataSet;

        control_system_displacement_double(:,model) = control_system_displacement.dataSet;

        control_system_output_double(:,model) = control_system_output.dataSet;

        %% Run the Control System Stimulus through the forward model (Paralyzed Side Simulation Model)

        predicted_healthy_displacement(:,model) = control_system_displacement_double(1000:end-1000,model);

        %Run Paralyzed Side Simulation Model
        amplitude_modulation(:,model) = control_system_output_double(:,model);
        amplitude_modulation_simulink = [(t_total(1,1:end-1999))' amplitude_modulation(1000:end-1000,model)];

        set_param('Paralyzed_Model_Simulink/Output Noise','Cov','set_output_noise_power')

        %Run Simulink;
        time = length(t_total(1,1:end-2000))/1000;
        out = sim('Paralyzed_Model_Simulink',time);

        input_stimulus = out.Paralyzed_Model_Stimulus;
        predicted_paralyzed_displacement(:,model) = out.Paralyzed_Model_Displacement;
        t_simulink = out.tout;

        Variance(trial,model) = vaf(predicted_healthy_displacement(:,model),predicted_paralyzed_displacement(:,model));
        control_system_validation_accuracy(trial,model) = Variance(trial,model);
        

    end

end

%% Calculate Mean and Std VAF
validation_accuracy_mean = [];
validation_accuracy_std = [];

for model = 1:num_models
    validation_accuracy_mean(:,model) = mean(control_system_validation_accuracy(:,model));
    validation_accuracy_std(:,model) = std(control_system_validation_accuracy(:,model));
end

validation_accuracy_mean
validation_accuracy_std

tEnd = toc(tStart)/60