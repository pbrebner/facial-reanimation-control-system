% Execute Control System for x number of trials
% Does not display figures 

%Run after identifying the EMG Response System (ERS) and inverse Stimulus
%Response System (SRS)^1

%Generates an input desired displacement, the corresponding nerual command
%and runs it through the Hammerstein Model followed by the Inverse of the
%LNL Model to get the stimulus

tStart = tic;

%%
%Set initial Parameters
% set_output_noise_power = 0;
% noise_snr_NHK = [];
% noise_snr_LNL = [];
% output_noise_power = [];
%figNum = 850;
time = 180;
figNum = 300;

Fs = 1000;
Nfft = 1000;

EMG_response_model = NHK2;
% stimulus_response_model = SRS_inverse;

PRBS_movement = false;
% physiological_movement = true;

%PRBS Stimulus
PRBS_movement_time = 180;
variable_amplitude = false;
N = PRBS_movement_time/10;
M = 10000;
PRBS_amplitude = 10; %mm

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_stimulus_max_amplitude = 0.01;
fr = 0.1;                 %Frequency distribution mean (Hz)
sig = 0.6;                %Std of Frequency Distribution (Hz)
W = 0.55;
nf = 18;                  %number of random signal changes
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

%%
% Simulated Input Test 2 (PRBS Input)
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

    stim_frequency = 50;
    stim_amplitude = desired_displacement*170;

    %input_stimulus = max(stim_amplitude.*square(2*pi*stim_frequency.*t_total),0);
    input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);

    %simulated_input_PRBS(:,signal) = input_stimulus';
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
% Identify a model between the simulated output (as an input) and the
% simulated input (as an output)

%     Zcur_simulated = [simulated_outputs(1000:end,1) simulated_input(1000:end,1)];
%     Zcur_simulated = nldat(Zcur_simulated,'domainIncr',0.001,'comment','Simulated Output (as Input), Simulated Input (as Output)','chanNames', {'Simulated Output (as Input)' 'Simulated Input (as Output)'});

%     test_output = [];

%     for signal = 2:num_signals
% 
%         Test_Output = nldat(simulated_outputs(1000:end,signal));
%         set(Test_Output,'domainIncr',0.001);
% 
%         test_output(:,signal-1) = Test_Output.dataSet;
% 
%     end

    
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
    % set(P,'polyType','tcheb');
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

        %%
        %Generate Neural Input Command Signal
        Amplitude = desired_displacement*100;  %mV
        Frequency = desired_displacement*14000;  %Hz

        neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

        neural_simulink = [t_total' neural];

        %%
        %Generate EMG Signal (Simulink)
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

        % set(EMG_response_model,'domainIncr',0.001)

        control_system_displacement = nlsim(EMG_response_model,EMG_simulink);
        set(control_system_displacement,'domainIncr',0.001);
        %control_system_displacement = control_system_displacement.dataSet;

        %Set domain increments of Stimulus response model
        % stimulus_response_model_IRF = stimulus_response_model{1,1};
        % set(stimulus_response_model_IRF, 'domainIncr', 1.0e-3);
        % stimulus_response_model{1,1} = stimulus_response_model_IRF;

        control_system_output = nlsim(SRS_inverse_all(:,model),control_system_displacement);
        %control_system_stimulus_output = control_system_stimulus_output.dataSet;


        %% Plot all in One Figure

    %     figure(figNum)
    %     figNum = figNum+1;
    % 
    %     subplot(3,2,1)
    %     plot(t_total,neural)
    %     title('Neural Command','Fontsize',14)
    %     xlabel('Time (s)','Fontsize',14)
    %     ylabel('Amplitude (% of MUs)', 'Fontsize', 14)
    %     grid on
    % 
    %     subplot(3,2,2)
    %     EMG_double = EMG_simulink.dataSet;
    %     plot(t_total,EMG_double)
    %     title('EMG Input (Simulink)','Fontsize',14)
    %     xlabel('Time (s)','Fontsize',14)
    %     ylabel('Amplitude (V)', 'Fontsize', 14)
    %     grid on
    % 
    %     subplot(3,2,[3 4])
        control_system_displacement_double(:,model) = control_system_displacement.dataSet;
    %     plot(t_total(1,1:end-1999),control_system_displacement_double(1000:end-1000,1));
    %     title('ERS Output (Healthy Side Predicted Displacement)', 'Fontsize', 14)
    %     % xlabel('Time (s)','Fontsize',14)
    %     ylabel('Displacement (m)', 'Fontsize', 14)
    %     grid on
    % 
    %     subplot(3,2,[5 6])
        control_system_output_double(:,model) = control_system_output.dataSet;
    %     plot(t_total(1,1:end-1999),control_system_output_double(1000:end-1000,:));
    %     title('Control System Output (Amplitude Modulation)', 'Fontsize', 14)
    %     xlabel('Time (s)','Fontsize',14)
    %     ylabel('Amplitude', 'Fontsize', 14)
    %     grid on


        %% Run the Control System Stimulus through the forward model (Paralyzed Side Simulation Model)
        % And perform some analysis on the results

        predicted_healthy_displacement(:,model) = control_system_displacement_double(1000:end-1000,model);

        %Run Paralyzed Side Simulation Model
        %stimulus_simulink = [(t_total(1,1:end-1999))' control_system_output_double(1000:end-1000,:)];
        amplitude_modulation(:,model) = control_system_output_double(:,model);
        amplitude_modulation_simulink = [(t_total(1,1:end-1999))' amplitude_modulation(1000:end-1000,model)];

        % set_param('Paralyzed_Model_Simulink_ControlSystemTest/Output Noise','Cov','set_output_noise_power')
        set_param('Paralyzed_Model_Simulink/Output Noise','Cov','set_output_noise_power')
        %output_noise_power = [output_noise_power set_output_noise_power];

        %Run Simulink;
        time = length(t_total(1,1:end-2000))/1000;
        % out = sim('Paralyzed_Model_Simulink_ControlSystemTest',time);
        out = sim('Paralyzed_Model_Simulink',time);

        input_stimulus = out.Paralyzed_Model_Stimulus;
        predicted_paralyzed_displacement(:,model) = out.Paralyzed_Model_Displacement;
        t_simulink = out.tout;

        %plot the stimulus that is created based on the amplitude modulation
        % figure(figNum)
        % figNum = figNum+1;
        % plot(t_total(1,1:end-1999),input_stimulus)
        % title('Stimulus Signal based on Control System Amplitude Modulation Output')
        % ylabel('Voltage (V)')
        % xlabel('Time (s)')
        % grid on

        %Compare the predicted healthy displacement and predicted paralyzed
        %displacement with the desired displacement
        % figure(figNum)
        % figNum = figNum+1;
        % subplot(4,1,1)
        % plot(t_total(1,1:end-1999),desired_displacement(:,1000:end-1000))
        % ax = gca;
        % ax.FontSize = 10;
        % title('(a) Desired Displacement', 'Fontsize', 12)
        % % xlabel('Time (s)', 'Fontsize', 12)
        % ylabel('Displacement (m)', 'Fontsize', 12)
        % grid on
        % 
        % subplot(4,1,2)
        % plot(t_total(1,1:end-1999),predicted_healthy_displacement)
        % ax = gca;
        % ax.FontSize = 10;
        % title('(b) Predicted Healthy Displacement (Output of ERS)', 'Fontsize', 12)
        % % xlabel('Time (s)', 'Fontsize', 12)
        % ylabel('Displacement (m)', 'Fontsize', 12)
        % grid on
        % 
        % subplot(4,1,3)
        % plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
        % ax = gca;
        % ax.FontSize = 10;
        % title('(c) Predicted Paralyzed Displacement (Based on Control System Amplitude Modulation Output)', 'Fontsize', 12)
        % % xlabel('Time (s)', 'Fontsize', 12)
        % ylabel('Displacement (m)', 'Fontsize', 12)
        % grid on
        % 
        % subplot(4,1,4)
        Variance(trial,model) = vaf(predicted_healthy_displacement(:,model),predicted_paralyzed_displacement(:,model));
        control_system_validation_accuracy(trial,model) = Variance(trial,model);
        % hold on
        % plot(t_total(1,1:end-1999),predicted_healthy_displacement)
        % plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
        % hold off
        % ax = gca;
        % ax.FontSize = 10;
        % title(['(d) Superimposed, VAF = ' num2str(Variance) '%'], 'Fontsize', 12)
        % xlabel('Time (s)', 'Fontsize', 12)
        % ylabel('Displacement (m)', 'Fontsize', 12)
        % legend('Predicted Healthy Displacement', 'Predicted Paralyzed Displacement','Fontsize',10)
        % grid on

        % % Compares the desired displacement with the predicted paralyzed
        % % displacement
        % figure(figNum)
        % figNum = figNum+1;
        % subplot(3,1,1)
        % plot(t_total(1,1:end-1999),desired_displacement(:,1000:end-1000))
        % title('(a) Desired Displacement', 'Fontsize', 12)
        % % xlabel('Time (s)', 'Fontsize', 12)
        % ylabel('Displacement (m)', 'Fontsize', 12)
        % grid on
        % 
        % subplot(3,1,2)
        % plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
        % title('(b) Predicted Paralyzed Displacement (Based on Control System Amplitude Modulation Output)', 'Fontsize', 12)
        % % xlabel('Time (s)', 'Fontsize', 12)
        % ylabel('Displacement (m)', 'Fontsize', 12)
        % grid on
        % 
        % subplot(3,1,3)
        % Variance2 = vaf((desired_displacement(:,1000:end-1000))',predicted_paralyzed_displacement);
        % hold on
        % plot(t_total(1,1:end-1999),desired_displacement(:,1000:end-1000))
        % plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
        % title(['(c) Superimposed, VAF = ' num2str(Variance2) '%'], 'Fontsize', 12)
        % xlabel('Time (s)', 'Fontsize', 12)
        % ylabel('Displacement (m)', 'Fontsize', 12)
        % legend('Desired Displacement', 'Predicted Paralyzed Displacement')
        % grid on

        % %PDF and Spectrum of Residuals (between predicted healthy displacement and
        % %predicted paralyzed displacement)
        % residuals = predicted_healthy_displacement - predicted_paralyzed_displacement;
        % 
        % Nfft = 1000;
        % residuals_zero = residuals - mean(residuals);
        % [Prr,f] = pwelch(residuals_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);
        % 
        % figure(figNum)
        % figNum = figNum+1;
        % sgtitle('Residuals between Predicted Healthy Side Displacement and Predicted Paralyzed Displacement','Fontsize',14)
        % 
        % subplot(2,2,[1 2])
        % plot(t_total(1,1:end-3998),residuals(1000:end-1000,:))
        % ax = gca;
        % ax.FontSize = 14;
        % title('(a) Residuals','Fontsize',16)
        % ylabel('Displacement (m)','Fontsize',16)
        % xlabel('Time (s)','Fontsize',16)
        % grid on
        % 
        % subplot(2,2,3)
        % histogram(residuals(1000:end-1000,:))
        % ax = gca;
        % ax.FontSize = 14;
        % title('(b) Residual Distribution','Fontsize',16)
        % xlabel('Displacement (m)','Fontsize',16)
        % ylabel('Density','Fontsize',16)
        % grid on
        % 
        % subplot(2,2,4)
        % semilogy(f(1:200,:),Prr(1:200,:))
        % ax = gca;
        % ax.FontSize = 14;
        % title('(c) Spectrum of Residuals','Fontsize',16)
        % ylabel('PSD (log scale)','Fontsize',16);
        % xlabel('Frequency (Hz)','Fontsize',16);
        % grid on;
        % 
        % figure(figNum)
        % figNum = figNum+1;
        % sgtitle('Residuals between Predicted Healthy Side Displacement and Predicted Paralyzed Displacement','Fontsize',14)
        % 
        % subplot(2,2,[1 2])
        % plot(t_total(1,1:end-3998),residuals(1000:end-1000,:))
        % ax = gca;
        % ax.FontSize = 14;
        % title('(a) Residuals','Fontsize',16)
        % ylabel('Displacement (m)','Fontsize',16)
        % xlabel('Time (s)','Fontsize',16)
        % grid on
        % 
        % subplot(2,2,3)
        % histogram(residuals(1000:end-1000,:))
        % ax = gca;
        % ax.FontSize = 14;
        % title('(b) Residual Distribution','Fontsize',16)
        % xlabel('Displacement (m)','Fontsize',16)
        % ylabel('Density','Fontsize',16)
        % grid on
        % 
        % subplot(2,2,4)
        % semilogy(f(1:61,:),Prr(1:61,:))
        % ax = gca;
        % ax.FontSize = 14;
        % title('(c) Spectrum of Residuals','Fontsize',16)
        % ylabel('PSD (log scale)','Fontsize',16);
        % xlabel('Frequency (Hz)','Fontsize',16);
        % grid on;

        % %Spectrum of Predicted Healthy Displacement and Predicted Paralyzed
        % %Displacement
        % predicted_healthy_displacement_zero = predicted_healthy_displacement - mean(predicted_healthy_displacement);
        % predicted_paralyzed_displacement_zero = predicted_paralyzed_displacement - mean(predicted_paralyzed_displacement);
        % 
        % Nfft = 1000;
        % [Pxx,~] = pwelch(predicted_healthy_displacement_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);
        % [Pyy,f] = pwelch(predicted_paralyzed_displacement_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);
        % 
        % figure(figNum)
        % figNum = figNum+1;
        % subplot(2,2,1),plot(t_total(1,1:end-1999),predicted_healthy_displacement)
        % title('Predicted Healthy Displacement','Fontsize',16)
        % xlabel('Time (s)','Fontsize',16)
        % ylabel('Displacement (m)','Fontsize',16)
        % grid on
        % 
        % subplot(2,2,2),plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
        % title('Predicted Paralyzed Displacement','Fontsize',16)
        % xlabel('Time (s)','Fontsize',16)
        % ylabel('Displacement (m)','Fontsize',16)
        % grid on
        % 
        % subplot(2,2,3),semilogy(f(1:200,:),Pxx(1:200,:));
        % title('Predicted Healthy Displacement Spectrum','Fontsize',16);
        % ylabel('PSD (log scale)','Fontsize',16);
        % xlabel('Frequency (Hz)','Fontsize',16);
        % grid on;
        % 
        % subplot(2,2,4),semilogy(f(1:200,:),Pyy(1:200,:));
        % title('Predicted Paralyzed Displacement Spectrum','Fontsize',16);
        % ylabel('PSD','Fontsize',16);
        % xlabel('Frequency (Hz)','Fontsize',16);
        % grid on;

        % %Auto-Correlation and Cross Correlation of Predicted Healthy Displacement
        % %and Predicted Paralyzed Displacement
        % figure(figNum)
        % figNum = figNum+1;
        % sgtitle('Correlation of Predicted Healthy Displacement and Predicted Paralyzed Displacement','Fontsize',14)
        % 
        % maxlag = floor(length(predicted_healthy_displacement_zero)*0.05);
        % [Rxx,lags] = xcorr(predicted_healthy_displacement_zero,maxlag,'coeff');
        % lags = lags/Fs;
        % subplot(3,1,1),plot(lags,Rxx)
        % title('(a) Auto-Correlation of Predicted Healthy Displacement','Fontsize',16)
        % % xlabel('Time (s)','Fontsize',16)
        % ylabel('Correlation','Fontsize',16)
        % grid on
        % 
        % [Ryy,lags] = xcorr(predicted_paralyzed_displacement_zero,maxlag,'coeff');
        % lags = lags/Fs;
        % subplot(3,1,2),plot(lags,Ryy)
        % title('(b) Auto-Correlation of Predicted Paralyzed Displacement','Fontsize',16)
        % % xlabel('Time (s)','Fontsize',16)
        % ylabel('Correlation','Fontsize',16)
        % grid on
        % 
        % [Rxy,lags] = xcorr(predicted_paralyzed_displacement_zero,predicted_healthy_displacement_zero,maxlag,'coeff');
        % lags = lags/Fs;
        % subplot(3,1,3),plot(lags,Rxy)
        % title('(c) Cross-Correlation','Fontsize',16)
        % xlabel('Time (s)','Fontsize',16)
        % ylabel('Correlation','Fontsize',16)
        % grid on
        % 
        % % Frequency Response between Predicted Healthy Side Displacement and
        % % Predicted Paralyzed Displacement
        % %Transfer Function
        % Nfft = 1000;
        % txy = tfestimate(predicted_healthy_displacement,predicted_paralyzed_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);
        % G = abs(txy);
        % G = 20*log10(G);
        % P = angle(txy);
        % P = P*(180/pi);
        % [Cxy, f] = mscohere(predicted_healthy_displacement,predicted_paralyzed_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);
        % 
        % figure(figNum)
        % figNum = figNum+1;
        % sgtitle('Frequency Response of Predicted Healthy Displacement and Predicted Paralyzed Displacement','Fontsize',14)
        % 
        % subplot(3,1,1),plot(f(1:31,:),G(1:31,:));
        % ax = gca;
        % ax.FontSize = 14;
        % title('(a) Gain','Fontsize',16);
        % ylabel('Gain (dB)','Fontsize',16);
        % % xlabel('Frequency (Hz)','Fontsize',20);
        % grid on;
        % 
        % subplot(3,1,2),plot(f(1:31,:),P(1:31,:));
        % ax = gca;
        % ax.FontSize = 14;
        % title('(b) Phase Shift','Fontsize',16);
        % ylabel('Phase (degrees)','Fontsize',16);
        % % xlabel('Frequency (Hz)','Fontsize',20);
        % grid on;
        % 
        % subplot(3,1,3),plot(f(1:31,:),Cxy(1:31,:));
        % ax = gca;
        % ax.FontSize = 14;
        % title('(c) Coherence','Fontsize',16);
        % ylabel('Coherence','Fontsize',16);
        % xlabel('Frequency (Hz)','Fontsize',16);
        % grid on;
        % 
        % figure(figNum)
        % figNum = figNum+1;
        % plot(f(1:11,:),Cxy(1:11,:));
        % ax = gca;
        % ax.FontSize = 16;
        % title('Coherence','Fontsize',20);
        % ylabel('Coherence','Fontsize',20);
        % xlabel('Frequency (Hz)','Fontsize',20);
        % grid on;
        % 
        % % Frequency Response between Desired Displacement and
        % % Predicted Paralyzed Displacement
        % %Transfer Function
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
        % % xlabel('Frequency (Hz)','Fontsize',20);
        % grid on;
        % 
        % subplot(3,1,2),plot(f(1:31,:),P(1:31,:));
        % ax = gca;
        % ax.FontSize = 14;
        % title('(b) Phase Shift','Fontsize',16);
        % ylabel('Phase (degrees)','Fontsize',16);
        % % xlabel('Frequency (Hz)','Fontsize',20);
        % grid on;
        % 
        % subplot(3,1,3),plot(f(1:31,:),Cxy(1:31,:));
        % ax = gca;
        % ax.FontSize = 14;
        % title('(c) Coherence','Fontsize',16);
        % ylabel('Coherence','Fontsize',16);
        % xlabel('Frequency (Hz)','Fontsize',16);
        % grid on;
        % 
        % % Plotting the residuals as a funtion of desired position
        % figure(figNum)
        % figNum = figNum+1;
        % plot(residuals, desired_displacement(:,1000:end-1000))
        % title('Desired Displacement vs Residuals','Fontsize',20);
        % ylabel('Displacement (m)','Fontsize',18);
        % xlabel('Residual (m)','Fontsize',18);
        % grid on;
        % 
        % figure(figNum)
        % figNum = figNum+1;
        % plot(residuals, predicted_healthy_displacement)
        % title('Predicted Healthy Displacement vs Residuals','Fontsize',20);
        % ylabel('Displacement (m)','Fontsize',18);
        % xlabel('Residual (m)','Fontsize',18);
        % grid on;
        % 
        % figure(figNum)
        % figNum = figNum+1;
        % plot(residuals, predicted_paralyzed_displacement)
        % title('Predicted Paralyzed Displacement vs Residuals','Fontsize',20);
        % ylabel('Displacement (m)','Fontsize',18);
        % xlabel('Residual (m)','Fontsize',18);
        % grid on;


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

%%
tEnd = toc(tStart)/60