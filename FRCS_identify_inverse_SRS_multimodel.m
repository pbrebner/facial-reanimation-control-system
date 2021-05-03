%% Control system: Identify the Inverse SRS

%Identify the Inverse SRS Model:
%   1. Simulate the response of the SRS you identify to a white input to generate a simulated output.
%   2. Then identify a model between the simulated output, as an input, and the simulated input, to estimate the inverse model.
%   3. By using the simulated signals you can avoid problems with noise.

tStart = tic;

%% Set Initial Parameters

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

    %LNL Model
    SRS_inverse_LNL = lnlbl;  %LNL Model
    set(SRS_inverse_LNL,'idMethod','hk','hkTolerance', 0.5,...
        'nhkMaxIts', 3, 'nhkMaxInner', 3);
    
    I1 = irf;
    set(I1,'nLags',1,'nSides',2,'domainIncr',0.001);
    P = polynom;
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

    input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);

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

tEnd = toc(tStart)/60