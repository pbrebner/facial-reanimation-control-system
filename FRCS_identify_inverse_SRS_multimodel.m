%% Control system: Identify the Inverse SRS

%Identify the Inverse SRS Model:
%   1. Simulate the response of the SRS you identify to a white input to generate a simulated output.
%   2. Then identify a model between the simulated output, as an input, and the simulated input, to estimate the inverse model.
%   3. By using the simulated signals you can avoid problems with noise.

tStart = tic;

%%
%Set Initial Parameters
time = 180;
set_output_noise_power = 0;
figNum = 300;
Fs = 1000;

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

num_signals_ident = 2;
num_signals_val = 100;

% white_input_type = false;  %white or PRBS


%%
%Generate a simulated white input for the identified SRS Model

% if white_input_type == true
%     
%     num_signals = 2;
%     
%     %mu = 0;
%     variable_amplitude = false;
%     N = 18;
%     M = 10000;
%     sig1_amplitude = 10;
%     
%     for signal = 1:num_signals
% 
%         sig1 = [0];                    %Intialize amplitude
%         if variable_amplitude == true   
%             for k = 1:N
%                 if k == 1
%                     R = sig1_amplitude;
%                 else
%                     R = rand(1,1)*sig1_amplitude;   %Randomly generate a number between 0 and 4
%                 end
% 
%                 for j = 1:M
%                     sig1 = [sig1 R];
%                 end
%             end
%             sig1 = sig1';
%         else
%             sig1 = sig1_amplitude;        %Constant Amplitude
%         end
% 
%         % sig1 = 4;
%         sig2 = 0.5;
% 
%         t_total = 0:0.001:time;
% 
%         simulated_input_uniform(:,signal) = ((rand(length(t_total),1)).*sig1)-(sig1/2);
%         simulated_input_gauss(:,signal) = ((randn(length(t_total),1))*sig2);
% 
%         %simulated_input_uniform_rect(:,signal) = max(simulated_input_uniform(:,signal),0);
%         simulated_input_uniform_rect(:,signal) = abs(simulated_input_uniform(:,signal));
%         simulated_input_gauss_rect(:,signal) = max(simulated_input_gauss(:,signal),0);
% 
%         % figure(figNum)
%         % figNum = figNum+1;
%         % plot(t_total,simulated_input)
%         % title('Simulated White Noise Input (Uniform)')
% 
%         %FFT of Desired Displacement
%         % L = length(simulated_input);
%         % 
%         % Y = fft(simulated_input);
%         % P1 = abs(Y/L);
%         % Pxx1 = P1(1:L/2+1);
%         % Pxx1(2:end-1) = 2*Pxx1(2:end-1);
%         % Pxx1 = Pxx1';
%         % 
%         % f1 = (Fs*(0:(L/2))/L)';
%         % figure(figNum)
%         % figNum = figNum+1;
%         % plot(f1,Pxx1) 
%         % title('FFT of Simulated Input (Uniform)')
%         % xlabel('f (Hz)')
%         % ylabel('|Pxx1(f)|')
% 
%         % figure(figNum)
%         % figNum = figNum+1;
%         % plot(t_total,simulated_input_gauss)
%         % title('Simulated White Noise Input (Guassian)')
% 
%         %FFT of Desired Displacement
%         % L = length(simulated_input);
%         % 
%         % Y = fft(simulated_input);
%         % P2 = abs(Y/L);
%         % Pxx2 = P2(1:L/2+1);
%         % Pxx2(2:end-1) = 2*Pxx2(2:end-1);
%         % Pxx2 = Pxx2';
%         % 
%         % f2 = (Fs*(0:(L/2))/L)';
%         % figure(figNum)
%         % figNum = figNum+1;
%         % plot(f2,Pxx2) 
%         % title('FFT of Simulated Input (Gaussian)')
%         % xlabel('f (Hz)')
%         % ylabel('|Pxx1(f)|')
% 
%         % figure(figNum)
%         % figNum = figNum+1;
%         % plot(t_total,simulated_input_rect)
%         % title('Simulated White Noise Input (Uniform-Rectified)')
% 
%         % figure(figNum)
%         % figNum = figNum+1;
%         % plot(t_total,simulated_input_gauss_rect)
%         % title('Simulated White Noise Input (Guassian-Rectified)')
%     end
% 
% end

%%
%Simulated Input Test (Using Simulink)

% %Run Simulink;
% out = sim('Simulated_Input_Simulink',time);
% 
% %Get Output from Simulink
% simulated_input_test = out.simulated_input;
% t_simulink = out.tout;
% 
% figure(figNum)
% figNum = figNum+1;
% plot(t_total,simulated_input_test)
% title('Simulated White Noise Input (Simulink)')
% 
% %FFT of Desired Displacement
% L = length(simulated_input_test);
% 
% Y = fft(simulated_input_test);
% P3 = abs(Y/L);
% Pxx3 = P3(1:L/2+1);
% Pxx3(2:end-1) = 2*Pxx3(2:end-1);
% Pxx3 = Pxx3';
% 
% f3 = (Fs*(0:(L/2))/L)';
% figure(figNum)
% figNum = figNum+1;
% plot(f3,Pxx3) 
% title('FFT of Simulated Input (Simulink)')
% xlabel('f (Hz)')
% ylabel('|Pxx1(f)|')

%%
% Simulated Input Test 2 (PRBS Input)
%PRBS Stimulus
% simulated_input_PRBS = [];
% 
% if white_input_type == false
% 
%     PRBS_stimulus_time = 300;
%     variable_amplitude = true;
%     N = PRBS_stimulus_time/10;
%     M = 10000;
%     PRBS_amplitude = 20; %mm
%     
%     simulated_input_PRBS_ident = [];
% 
%     for signal_ident = 1:num_signals_ident
% 
%         t_total = 0:0.001:PRBS_stimulus_time;
%         time = PRBS_stimulus_time;
% 
%         A = [0];                    %Intialize amplitude
%         if variable_amplitude == true   
%             for k = 1:N
%                 if k == 1
%                     R = PRBS_amplitude;
%                 else
%                     R = rand(1,1)*PRBS_amplitude;   %Randomly generate a number between 0 and 10
%                 end
% 
%                 for j = 1:M
%                     A = [A R];
%                 end
%             end
%         else
%             A = PRBS_amplitude;              %Constant Amplitude
%         end
% 
%         Range = [0,0.001]; %Specify that the single-channel PRBS value switches between -2 and 2
% 
%         %Specify the clock period of the signal as 1 sample. 
%         %That is, the signal value can change at each time step. 
%         %For PRBS signals, the clock period is specified in Band = [0 B], 
%         %where B is the inverse of the required clock period
%         %(Must be less than 1)
%         Band = [0 0.01];
% 
%         %Generate a nonperiodic PRBS of length 100 samples.
%         u = idinput(time*1000+1,'prbs',Band,Range);
% 
%         %Create an iddata object from the generated signal. 
%         %For this example, specify the sample time as 1 second.
%         u = iddata([],u,0.001);
% 
%         U = (u.InputData)';
%         desired_displacement = A.*U;
% 
%         stim_frequency = 50;
%         stim_amplitude = desired_displacement*170;
% 
%         %input_stimulus = max(stim_amplitude.*square(2*pi*stim_frequency.*t_total),0);
%         input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);
% 
%         %simulated_input_PRBS(:,signal) = input_stimulus';
%         simulated_input_PRBS_ident(:,signal_ident) = stim_amplitude;
% 
%     end
%     
%     PRBS_stimulus_time = 180;
%     variable_amplitude = true;
%     N = PRBS_stimulus_time/10;
%     M = 10000;
%     PRBS_amplitude = 20; %mm
%     
%     simulated_input_PRBS_val = [];
%     
%     for signal_val = 1:num_signals_val
% 
%         t_total = 0:0.001:PRBS_stimulus_time;
%         time = PRBS_stimulus_time;
% 
%         A = [0];                    %Intialize amplitude
%         if variable_amplitude == true   
%             for k = 1:N
%                 if k == 1
%                     R = PRBS_amplitude;
%                 else
%                     R = rand(1,1)*PRBS_amplitude;   %Randomly generate a number between 0 and 10
%                 end
% 
%                 for j = 1:M
%                     A = [A R];
%                 end
%             end
%         else
%             A = PRBS_amplitude;              %Constant Amplitude
%         end
% 
%         Range = [0,0.001]; %Specify that the single-channel PRBS value switches between -2 and 2
% 
%         %Specify the clock period of the signal as 1 sample. 
%         %That is, the signal value can change at each time step. 
%         %For PRBS signals, the clock period is specified in Band = [0 B], 
%         %where B is the inverse of the required clock period
%         %(Must be less than 1)
%         Band = [0 0.01];
% 
%         %Generate a nonperiodic PRBS of length 100 samples.
%         u = idinput(time*1000+1,'prbs',Band,Range);
% 
%         %Create an iddata object from the generated signal. 
%         %For this example, specify the sample time as 1 second.
%         u = iddata([],u,0.001);
% 
%         U = (u.InputData)';
%         desired_displacement = A.*U;
% 
%         stim_frequency = 50;
%         stim_amplitude = desired_displacement*170;
% 
%         %input_stimulus = max(stim_amplitude.*square(2*pi*stim_frequency.*t_total),0);
%         input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);
% 
%         %simulated_input_PRBS(:,signal) = input_stimulus';
%         simulated_input_PRBS_val(:,signal_val) = stim_amplitude;
% 
%     end
%     
% end


%%
%Simulate the reponse of the SRS model to the white input

% if white_input_type == true
%     simulated_input = simulated_input_uniform_rect;
% else
%     simulated_input_ident = simulated_input_PRBS_ident;
%     simulated_input_val = simulated_input_PRBS_val;
% end

% simulated_output_ident = [];
% simulated_outputs_ident = [];
% simulated_output_val = [];
% simulated_outputs_val = [];

% for signal_ident = 1:num_signals_ident
% 
%     simulated_output_ident = nlsim(SRS_model,simulated_input_ident(:,signal_ident));
%     set(simulated_output_ident, 'domainIncr',0.001)
%     
%     simulated_output_ident = max(simulated_output_ident,0);
% 
%     simulated_outputs_ident(:,signal_ident) = simulated_output_ident.dataSet;
% 
% end

% for signal_val = 1:num_signals_val
% 
%     simulated_output_val = nlsim(SRS_model,simulated_input_val(:,signal_val));
%     set(simulated_output_val, 'domainIncr',0.001)
%     
%     simulated_output_val = max(simulated_output_val,0);
% 
%     simulated_outputs_val(:,signal_val) = simulated_output_val.dataSet;
% 
% end


%%
% Identify a model between the simulated output (as an input) and the
% simulated input (as an output)
% Zcur_simulated_ident = [];
% 
% for signal_ident = 1:num_signals_ident
%     
%     Zcur_simulated = [simulated_outputs_ident(1000:end,signal_ident) simulated_input_ident(1000:end,signal_ident)];
%     Zcur_simulated = nldat(Zcur_simulated,'domainIncr',0.001,'comment','Simulated Output (as Input), Simulated Input (as Output)','chanNames', {'Simulated Output (as Input)' 'Simulated Input (as Output)'});
%     
%     Zcur_simulated_ident = [Zcur_simulated_ident Zcur_simulated];
%     
% end

% test_output_val = [];
% 
% for signal_val = 1:num_signals_val
% 
%     Test_Output = nldat(simulated_outputs_val(1000:end,signal_val));
%     set(Test_Output,'domainIncr',0.001);
%     
%     test_output_val(:,signal_val) = Test_Output.dataSet;
% 
% end


%%
% Test the signals on all the models
simulated_output_ident = [];
simulated_outputs_ident = [];

% Identification Signal Parameters
PRBS_stimulus_time = 480;
variable_amplitude = true;
N = PRBS_stimulus_time/10;
M = 10000;
PRBS_amplitude = 20; %mm

simulated_input_PRBS_ident = [];

%Number of Lags used by the models
nLags = 400;

accuracy_identification_LNL = [];
accuracy_validation_LNL = [];
accuracy_identification_Hamm = [];
accuracy_validation_Hamm = [];
accuracy_identification_Weiner = [];
accuracy_validation_Weiner = [];
accuracy_identification_IRF = [];
accuracy_validation_IRF = [];

iter_ident = 0;

for signal_ident = 1:num_signals_ident
    
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
    simulated_input_PRBS_ident = stim_amplitude;
    simulated_input_ident = simulated_input_PRBS_ident';
    
    simulated_output_ident = nlsim(SRS_model,simulated_input_ident);
    set(simulated_output_ident, 'domainIncr',0.001);
    
    simulated_output_ident = max(simulated_output_ident,0);
    simulated_outputs_ident = simulated_output_ident.dataSet;

    Zcur_simulated_ident = [simulated_outputs_ident(1000:end,:) simulated_input_ident(1000:end,:)];
    Zcur_simulated_ident = nldat(Zcur_simulated_ident,'domainIncr',0.001,'comment','Simulated Output (as Input), Simulated Input (as Output)','chanNames', {'Simulated Output (as Input)' 'Simulated Input (as Output)'});
    
%     SRS_inverse_LNL = [];
%     SRS_inverse_Hamm = [];
%     SRS_inverse_Weiner = [];
%     SRS_inverse_IRF = [];

    %LNL Model
    SRS_inverse_LNL = lnlbl;  %LNL Model
    set(SRS_inverse_LNL,'idMethod','hk','hkTolerance', 0.5,...
        'nhkMaxIts', 3, 'nhkMaxInner', 3);
    
    I1 = irf;
    set(I1,'nLags',1,'nSides',2,'domainIncr',0.001);
    P = polynom;
    % set(P,'polyType','tcheb');
    I3 = irf;
    set(I3,'nLags',nLags,'nSides', 2,'domainIncr',0.001);
    SRS_inverse_LNL.elements = {I1 P I3};
    
    SRS_inverse_LNL = nlident(SRS_inverse_LNL,Zcur_simulated_ident);
    
    figure(1)
    [R, V, yp] = nlid_resid(SRS_inverse_LNL,Zcur_simulated_ident);
    
    accuracy_identification_LNL = [accuracy_identification_LNL V];
    
    %Hammerstein
    SRS_inverse_Hamm = nlbl;  
    set(SRS_inverse_Hamm,'idMethod','hk','displayFlag',true,'threshNSE',.001);
    I2 = irf;
    set(I2,'nLags',nLags, 'nSides', 2); % Set number of lags and Sides in IRF
    SRS_inverse_Hamm{1,2} = I2;
    
    SRS_inverse_Hamm = nlident(SRS_inverse_Hamm,Zcur_simulated_ident);
    
    figure(2)
    [R, V, yp] = nlid_resid(SRS_inverse_Hamm,Zcur_simulated_ident);
    
    accuracy_identification_Hamm = [accuracy_identification_Hamm V];
    
    %Wiener
    SRS_inverse_Weiner = lnbl;
    set(SRS_inverse_Weiner,'idMethod','hk');
    I1 = irf;
    set(I1,'nLags',nLags,'nSides',2); % Set Number of lags and Sides in IRF
    SRS_inverse_Weiner{1,1} = I1;
    
    SRS_inverse_Weiner = nlident(SRS_inverse_Weiner,Zcur_simulated_ident);
    
    figure(3)
    [R, V, yp] = nlid_resid(SRS_inverse_Weiner,Zcur_simulated_ident);
    
    accuracy_identification_Weiner = [accuracy_identification_Weiner V];

    %Two-Sided IRF
    SRS_inverse_IRF = irf(Zcur_simulated_ident,'nLags',nLags,'nSides',2);

    figure(4)
    [R, V, yp] = nlid_resid(SRS_inverse_IRF,Zcur_simulated_ident);
    
    accuracy_identification_IRF = [accuracy_identification_IRF V];
    
    iter_ident = iter_ident+1
end


%% Inverse SRS Validation
simulated_output_val = [];
simulated_outputs_val = [];

% Set the Parametets of the Validation Signals
PRBS_stimulus_time = 180;
variable_amplitude = true;
N = PRBS_stimulus_time/10;
M = 10000;
PRBS_amplitude = 20; %mm

simulated_input_PRBS_val = [];
control_system_test_output_double = [];
iter_val = 0;

for signal_val = 1:num_signals_val
    
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
    simulated_input_PRBS_val = stim_amplitude;
    simulated_input_val = simulated_input_PRBS_val';
    
    simulated_output_val = nlsim(SRS_model,simulated_input_val);
    set(simulated_output_val, 'domainIncr',0.001);
    
    simulated_output_val = max(simulated_output_val,0);
    simulated_outputs_val = simulated_output_val.dataSet;
    
    Test_Output = nldat(simulated_outputs_val(1000:end,:));
    set(Test_Output,'domainIncr',0.001);
    
    test_output_val = Test_Output.dataSet;
    
    %LNL Model
    control_system_test_output = nlsim(SRS_inverse_LNL,test_output_val);
    control_system_test_output_double(:,1) = control_system_test_output.dataSet;
    inverse_SRS_validation_accuracy = vaf(simulated_input_val(1999:end-1000,:),control_system_test_output_double(1000:end-1000,1));
    
    accuracy_validation_LNL = [accuracy_validation_LNL inverse_SRS_validation_accuracy];
    
    %Hammerstein Model
    control_system_test_output = nlsim(SRS_inverse_Hamm,test_output_val);
    control_system_test_output_double(:,1) = control_system_test_output.dataSet;
    inverse_SRS_validation_accuracy = vaf(simulated_input_val(1999:end-1000,:),control_system_test_output_double(1000:end-1000,1));
    
    accuracy_validation_Hamm = [accuracy_validation_Hamm inverse_SRS_validation_accuracy];
    
    %Weiner Model
    control_system_test_output = nlsim(SRS_inverse_Weiner,test_output_val);
    control_system_test_output_double(:,1) = control_system_test_output.dataSet;
    inverse_SRS_validation_accuracy = vaf(simulated_input_val(1999:end-1000,:),control_system_test_output_double(1000:end-1000,1));
    
    accuracy_validation_Weiner = [accuracy_validation_Weiner inverse_SRS_validation_accuracy];
    
    %IRF Model
    control_system_test_output = nlsim(SRS_inverse_IRF,test_output_val);
    control_system_test_output_double(:,1) = control_system_test_output.dataSet;
    inverse_SRS_validation_accuracy = vaf(simulated_input_val(1999:end-1000,:),control_system_test_output_double(1000:end-1000,1));
    
    accuracy_validation_IRF = [accuracy_validation_IRF inverse_SRS_validation_accuracy];
    
    iter_val = iter_val+1
end

%% Mean and Std of Model Accuracies
%Identification Accuracy
accuracy_identification_LNL_mean = mean(accuracy_identification_LNL)
accuracy_identification_LNL_std = std(accuracy_identification_LNL)

accuracy_identification_Hamm_mean = mean(accuracy_identification_Hamm)
accuracy_identification_Hamm_std = std(accuracy_identification_Hamm)

accuracy_identification_Weiner_mean = mean(accuracy_identification_Weiner)
accuracy_identification_Weiner_std = std(accuracy_identification_Weiner)

accuracy_identification_IRF_mean = mean(accuracy_identification_IRF)
accuracy_identification_IRF_std = std(accuracy_identification_IRF)


%Validation Accuracy
accuracy_validation_LNL_mean = mean(accuracy_validation_LNL)
accuracy_validation_LNL_std = std(accuracy_validation_LNL)

accuracy_validation_Hamm_mean = mean(accuracy_validation_Hamm)
accuracy_validation_Hamm_std = std(accuracy_validation_Hamm)

accuracy_validation_Weiner_mean = mean(accuracy_validation_Weiner)
accuracy_validation_Weiner_std = std(accuracy_validation_Weiner)

accuracy_validation_IRF_mean = mean(accuracy_validation_IRF)
accuracy_validation_IRF_std = std(accuracy_validation_IRF)


    %% Plot all in One Figure (IRF)
%     figure(figNum)
%     figNum = figNum+1;
%     
%     subplot(3,2,1)
%     Zcur_simulated_double = Zcur_simulated.dataSet;
%     Zcur_input = Zcur_simulated_double(:,1);
%     plot(t_total(1,1:end-999),Zcur_input)
%     title('Simulated Displacement (Used as Input)','Fontsize', 10)
%     xlabel('Time (s)', 'Fontsize', 10)
%     ylabel('Displacement (m)','Fontsize', 10)
%     grid on
%     
%     subplot(3,2,2)
%     Zcur_output = Zcur_simulated_double(:,2);
%     plot(t_total(1,1:end-999),Zcur_output)
%     title('Simulated Amplitude Modulation (Used as Output)','Fontsize', 10)
%     xlabel('Time (s)', 'Fontsize', 10)
%     ylabel('Amplitude (V)','Fontsize', 10)
%     grid on
%     
%     subplot(3,2,[3 4])
%     plot(SRS_inverse)
%     title('Two-Sided Linear IRF','Fontsize', 10)
%     xlabel('Lags (s)', 'Fontsize', 10)
%     ylabel('X1','Fontsize', 10)
%     grid on
%     
%     subplot(3,2,5)
%     pred = double(yp);
%     hold on
%     plot(t_total(1,1:end-999),pred);
%     plot(t_total(1,1:end-999), Zcur_output)
%     hold off
%     V = vaf(Zcur_output(15000:end-14999,:),pred(15000:end-14999,:));
%     title(['Superimposed, VAF = ' num2str(V) '%'], 'Fontsize', 10)
%     xlabel('Time (s)', 'Fontsize', 10)
%     ylabel('Amplitude (V)', 'Fontsize', 10)
%     legend('Predicted', 'Observed', 'Fontsize', 8)
%     grid on
%     
%     subplot(3,2,6)
%     plot(R)
%     title('Residuals','Fontsize', 10)
%     xlabel('Time (s)','Fontsize', 10)
%     ylabel('Amplitude (V)','Fontsize', 10)
%     grid on

    %% Plot all in One Figure (Hammerstein or Wiener)
%     figure(figNum)
%     figNum = figNum+1;
% 
%     subplot(3,2,1)
%     Zcur_simulated_double = Zcur_simulated.dataSet;
%     Zcur_input = Zcur_simulated_double(:,1);
%     plot(t_total(1,1:end-999),Zcur_input)
%     title('(a) Simulated Displacement (Used as Input)','Fontsize', 10)
%     xlabel('Time (s)', 'Fontsize', 10)
%     ylabel('Displacement (m)','Fontsize', 10)
%     grid on
% 
%     subplot(3,2,2)
%     Zcur_output = Zcur_simulated_double(:,2);
%     plot(t_total(1,1:end-999),Zcur_output)
%     title('(b) Simulated Amplitude Modulation (Used as Output)','Fontsize', 10)
%     xlabel('Time (s)', 'Fontsize', 10)
%     ylabel('Amplitude (V)','Fontsize', 10)
%     grid on
% 
%     subplot(3,2,3)
%     plot(SRS_inverse{1,1})
%     title('(c) Static Nonlinear Element','Fontsize', 10)
%     xlabel('Displacement Input (m)', 'Fontsize', 10)
%     ylabel('Amplitude Modulation Output (V)','Fontsize', 10)
%     grid on
% 
%     subplot(3,2,4)
%     plot(SRS_inverse{1,2})
%     title('(d) Linear Element','Fontsize', 10)
%     xlabel('Lags (s)', 'Fontsize', 10)
%     ylabel('X1','Fontsize', 10)
%     grid on
% 
%     subplot(3,2,5)
%     pred = double(yp);
%     hold on
%     plot(t_total(1,1:end-2998),pred(1000:end-1000,:));
%     plot(t_total(1,1:end-2998), Zcur_output(1000:end-1000,:))
%     hold off
%     V = vaf(Zcur_output(15000:end-14999,:),pred(15000:end-14999,:));
%     title(['(e) Superimposed, VAF = ' num2str(V) '%'], 'Fontsize', 10)
%     xlabel('Time (s)', 'Fontsize', 10)
%     ylabel('Amplitude (V)', 'Fontsize', 10)
%     legend('Predicted', 'Observed', 'Fontsize', 8)
%     grid on
% 
%     subplot(3,2,6)
%     plot(R(1000:end-1000,:))
%     title('(f) Residuals','Fontsize', 10)
%     xlabel('Time (s)','Fontsize', 10)
%     ylabel('Amplitude (V)','Fontsize', 10)
%     grid on
    
%     accuracy_identification = [accuracy_identification V];

    % figure(figNum)
    % figNum = figNum+1;
    % subplot(1,2,1)
    % plot(SRS_inverse{1,1})
    % ax = gca;
    % ax.FontSize = 14;
    % title('(a) Static Nonlinear Element','Fontsize', 18)
    % xlabel('Displacement Input (m)', 'Fontsize', 18)
    % ylabel('Amplitude Modulation Output (V)','Fontsize', 18)
    % grid on
    % 
    % subplot(1,2,2)
    % plot(SRS_inverse{1,2})
    % ax = gca;
    % ax.FontSize = 14;
    % title('(b) Linear Element','Fontsize', 18)
    % xlabel('Lags (s)', 'Fontsize', 18)
    % ylabel('X1','Fontsize', 18)
    % grid on

    %% Analysis of Inverse SRS Identification Residuals
    % figure(figNum)
    % figNum = figNum+1;
    % sgtitle('Analysis of Inverse SRS Identification Residuals','Fontsize',14)
    % 
    % subplot(2,2,[1 2])
    % R_temp = double(R);
    % plot(R(1000:end-1000,:))
    % ax = gca;
    % ax.FontSize = 13;
    % title('(a) Residuals','Fontsize',18)
    % ylabel('Amplitude (V)','Fontsize',18)
    % xlabel('Time (s)','Fontsize',18)
    % grid on
    % 
    % subplot(2,2,3)
    % histogram(R_temp(1000:end-1000,:))
    % ax = gca;
    % ax.FontSize = 13;
    % title('(b) Residual Distribution','Fontsize',18)
    % xlabel('Amplitude (V)','Fontsize',18)
    % ylabel('Density','Fontsize',18)
    % grid on
    % 
    % S = spect(R);
    % subplot(2,2,4)
    % S_frequency = 0.0556:0.0556:0.0556*length(S);
    % subplot(2,2,4)
    % semilogy(S_frequency(:,1:500),S(1:500,:),'LineWidth',1.5);
    % ax = gca;
    % ax.FontSize = 13;
    % title('(c) Residual Power Spectrum','Fontsize',18);
    % ylabel('PSD (log scale)','Fontsize',18); 
    % xlabel('Frequency (Hz)','Fontsize',18);
    % grid on


    %% Inverse SRS Validation

%     control_system_test_output_double = [];
%     inverse_SRS_validation_accuracy = [];
% 
%     for signal = 1:num_signals-1
% 
%         control_system_test_output = nlsim(SRS_inverse,test_output(:,signal));
% 
%         control_system_test_output_double(:,signal) = control_system_test_output.dataSet;
% 
%         inverse_SRS_validation_accuracy(signal,:) = vaf(simulated_input(1999:end-1000,signal+1),control_system_test_output_double(1000:end-1000,signal));
% 
%     end
%     
%     accuracy_validation = [accuracy_validation inverse_SRS_validation_accuracy(1,:)];

    % Mean and Standard Deviation of Inverse SRS Validation Accuracy
%     validation_accuracy_mean = mean(inverse_SRS_validation_accuracy)
%     validation_accuracy_std = std(inverse_SRS_validation_accuracy)

    %Plot one validation example
%     figure(figNum)
%     figNum = figNum+1;
%     hold on
%     plot(t_total(1,1:end-2998),control_system_test_output_double(1000:end-1000,1))
%     plot(t_total(1,1:end-2998),simulated_input(1999:end-1000,2))
%     hold off
%     title(['Superimposed, VAF = ' num2str(inverse_SRS_validation_accuracy(1,1)) '%'], 'Fontsize', 20)
%     xlabel('Time (s)', 'Fontsize', 20)
%     ylabel('Amplitude (V)', 'Fontsize', 20)
%     legend('Predicted', 'Observed', 'Fontsize', 16)
%     grid on


%     nLags_plot = [nLags_plot nLags];
%     nLags = nLags - 100;
% 
% 
% 
% figure(figNum)
% figNum = figNum+1;
% hold on
% plot(nLags_plot,accuracy_identification,'Linewidth',2)
% plot(nLags_plot,accuracy_validation,'Linewidth',2)
% hold off
% title('Accuracy vs NLags','Fontsize',20)
% xlabel('NLags', 'Fontsize',18)
% ylabel('Accuracy (%VAF)','Fontsize',18)
% legend('Identification','Validation','Fontsize',16)
% grid on

%% Plot all in One Figure (LNL)
% figure(figNum)
% figNum = figNum+1;
% 
% Zcur_simulated_double = Zcur_simulated.dataSet;
% 
% subplot(3,3,1)
% Zcur_simulated_double = Zcur_simulated.dataSet;
% Zcur_input = Zcur_simulated_double(:,1);
% plot(t_total(1,1000:end),Zcur_input)
% title('Simulated Displacement (Used as Input)')
% xlabel('Time (s)', 'Fontsize', 10)
% ylabel('Displacement (m)','Fontsize', 10)
% 
% subplot(3,3,3)
% Zcur_output = Zcur_simulated_double(:,2);
% plot(t_total(1,1000:end),Zcur_output)
% title('Simulated Amplitude Modulation (Used as Output)')
% xlabel('Time (s)', 'Fontsize', 10)
% ylabel('Amplitude (V)','Fontsize', 10)
% 
% subplot(3,3,4)
% plot(SRS_inverse{1,1})
% 
% subplot(3,3,5)
% plot(SRS_inverse{1,2})
% 
% subplot(3,3,6)
% plot(SRS_inverse{1,3})
% 
% subplot(3,3,7)
% pred = double(yp);
% hold on
% plot(t_total(1,1000:end), Zcur_output)
% plot(t_total(1,1000:end),pred);
% hold off
% title(['Superimposed, VAF = ' num2str(V) '%'], 'Fontsize', 10)
% xlabel('Time (s)', 'Fontsize', 10)
% ylabel('Amplitude (V)', 'Fontsize', 10)
% legend('Observed', 'Predicted', 'Fontsize', 6)
% 
% subplot(3,3,9)
% plot(R)
% title('Residuals','Fontsize', 10)
% xlabel('Time (s)','Fontsize', 10)
% ylabel('Amplitude (V)','Fontsize', 10)
% grid on


%%
tEnd = toc(tStart)/60