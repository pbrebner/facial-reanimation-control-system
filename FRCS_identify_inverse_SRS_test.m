%% Facial Reanimation Control system (FRCS): 

%INSTRUCTIONS: Must Run FRCS_identify_models before running this script

%

%Identifying the Inverse SRS Model:
%   1. Simulate the response of the SRS you identified to a PRBS input to generate a simulated output.
%   2. Then identify an inverse SRS model between the simulated output, as
%   input, and the simulated PRBS input, as output.

%% User Input Prompts

if compare_two_models == true
    prompt1 = 'Which SRS Model do you want to Inverse? PRBS/Phys [Phys]: ';
    str1 = input(prompt1,'s');
    if ~strcmp(str1,'PRBS') & ~strcmp(str1,'Phys') & ~isempty(str1)
        disp('Invalid Input')
        return
    elseif isempty(str1)
        str1 = 'Phys';
    end
    
    if strcmp(str1,'PRBS')
        if strcmp(model_type,'IRF')
            SRS = SRS_all{1};
        else
            SRS = SRS_all(1);
        end
        disp(['Inversing ' model_type ' SRS Model Identified from PRBS Input'])
    elseif strcmp(str1,'Phys')
        if strcmp(model_type,'IRF')
            SRS = SRS_all{2};
        else
            SRS = SRS_all(2);
        end
        disp(['Inversing ' model_type ' SRS Model Identified from Physiological Input'])
    end
    SRS_inverse_input_type = str1;
else
    if strcmp(model_type,'IRF')
        SRS = SRS_all{1};
    else
        SRS = SRS_all(1);
    end
    disp(['Inversing ' model_type ' SRS Model Identified from ' input_type ' Input'])
    SRS_inverse_input_type = input_type;
end

prompt2 = 'Inverse SRS Model Structure? LNL/Hammerstein(Hamm)/Wiener/IRF [Hamm]: ';
str2 = input(prompt2,'s');
if ~strcmp(str2,'LNL') & ~strcmp(str2,'Hamm') & ~strcmp(str2,'Wiener') & ~strcmp(str2,'IRF') & ~isempty(str2)
    disp('Invalid Input')
    return
elseif isempty(str2)
    str2 = 'Hamm';
end
SRS_inverse_model_structure = str2;

prompt3 = 'Number of Validation Trials? 1-30 [1]: ';
str3 = input(prompt3);
if str3<1 | str3>30
    disp('Invalid Input')
    return
elseif isempty(str3)
    str3 = 1;
end

tStart = tic;

%% Set Initial Parameters

figNum = 100;
Fs = 1000;

%Simulated Input (PRBS) Parameters
PRBS_stimulus_time = 480;
variable_amplitude = true;      %Can either be constant or variable amplitude
N = PRBS_stimulus_time/10;
M = 10000;
PRBS_amplitude = 20;            %PRBS Signal Amplitude (mm)

%Number of Signals (Identification and Validation)
num_signals = str3+1;

%Inverse Model Structure
SRS_inverse_LNL = false;
SRS_inverse_Hamm = false;
SRS_inverse_Weiner = false;
SRS_inverse_IRF = false;

if strcmp(str2,'LNL')
    SRS_inverse_LNL = true;
elseif strcmp(str2,'Hamm')
    SRS_inverse_Hamm = true;
elseif strcmp(str2,'Wiener')
    SRS_inverse_Weiner = true;
elseif strcmp(str2,'IRF')
    SRS_inverse_IRF = true;
end

num_signals_ident = 2;
num_signals_val = 100;

multi_lag = true;
multi_model = false;

%% Multi-Lag or Multi-Model

if multi_lag == true
    
    %Simulated Input (PRBS Input)

    simulated_input = [];

    for signal = 1:num_signals

        t_total = 0:0.001:PRBS_stimulus_time;
        time = PRBS_stimulus_time;

        A = [0];                                    %Intialize amplitude
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

        %Generate a nonperiodic PRBS of length time.
        u = idinput(time*1000+1,'prbs',Band,Range);

        %Create an iddata object from the generated signal. 
        %For this example, specify the sample time as 0.001 second.
        u = iddata([],u,0.001);

        U = (u.InputData)';
        desired_displacement = A.*U;

        stim_frequency = 50;
        stim_amplitude = desired_displacement*170;

        input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);

        simulated_input(:,signal) = stim_amplitude;

    end
    
    %% Simulate the reponse of the SRS model to the Input Realizations
    %Simulated Outputs of the SRS Model
    simulated_output = [];
    simulated_outputs = [];

    %Runs through all input realizations
    for signal = 1:num_signals

        simulated_output = nlsim(SRS,simulated_input(:,signal));
        set(simulated_output, 'domainIncr',0.001);

        simulated_output = max(simulated_output,0);
        simulated_outputs(:,signal) = simulated_output.dataSet;

    end
    
    %% Identify a model between a simulated output (as input) and a simulated input (as output)

    %Input/Output for Inverse SRS Identification
    Zcur_simulated = [simulated_outputs(1000:end,1) simulated_input(1000:end,1)];
    Zcur_simulated = nldat(Zcur_simulated,'domainIncr',0.001,'comment','Simulated Output (as Input), Simulated Input (as Output)','chanNames', {'Simulated Output (as Input)' 'Simulated Input (as Output)'});

    %plot the Input/Output
    figure(figNum)
    figNum = figNum+1;
    subplot(2,1,1)
    plot(t_total(1,1:end-999),simulated_input(1000:end,1))
    ax = gca;
    ax.FontSize = 15;
    title('Simulated Amplitude Modulation, A(t)','Fontsize',20)
    ylabel('Amplitude (V)','Fontsize',18)
    grid on

    subplot(2,1,2)
    plot(t_total(1,1:end-999),simulated_outputs(1000:end,1))
    ax = gca;
    ax.FontSize = 15;
    title('Simulated Paralyzed Displacement, Pos_P(t)','Fontsize',20)
    ylabel('Displacement (m)','Fontsize',18)
    xlabel('Time (s)','Fontsize',18)
    grid on

    %Outputs to other input realizations 
    validation_output = [];

    for signal = 2:num_signals

        Validation_Output = nldat(simulated_outputs(1000:end,signal));
        set(Validation_Output,'domainIncr',0.001);

        validation_output(:,signal-1) = Validation_Output.dataSet;

    end
    
    %% Identify the Inverse SRS
    
    nLags = 400;
    nLags_plot = [];

    accuracy_identification = [];
    accuracy_validation = [];
    
    while nLags > 50
        
        if SRS_inverse_LNL == true

            SRS_inverse = lnlbl;  %LNL Model
            set(SRS_inverse,'idMethod','hk','hkTolerance', 0.1,...
                'nhkMaxIts', 4, 'nhkMaxInner', 4);

            I1 = irf;
            set(I1,'nLags',1,'nSides',2,'domainIncr',0.001);
            P = polynom;
            I3 = irf;
            set(I3,'nLags',nLags,'nSides', 2,'domainIncr',0.001);
            SRS_inverse.elements = {I1 P I3};

            SRS_inverse = nlident(SRS_inverse,Zcur_simulated);

        elseif SRS_inverse_Hamm == true

            SRS_inverse = nlbl;  %Hammerstein
            set(SRS_inverse,'idMethod','hk','displayFlag',true,'threshNSE',.001);
            I2 = irf;
            set(I2,'nLags',nLags, 'nSides', 2,'domainIncr',0.001); % Set number of lags and Sides in IRF
            SRS_inverse{1,2} = I2;

            SRS_inverse = nlident(SRS_inverse,Zcur_simulated);

        elseif SRS_inverse_Weiner == true

            SRS_inverse = lnbl; %Wiener
            set(SRS_inverse,'idMethod','hk');
            I1 = irf;
            set(I1,'nLags',nLags,'nSides',2,'domainIncr',0.001); % Set Number of lags and Sides in IRF
            SRS_inverse{1,1} = I1;

            SRS_inverse = nlident(SRS_inverse,Zcur_simulated);

        elseif SRS_inverse_IRF == true

            %Two-Sided IRF
            SRS_inverse = irf(Zcur_simulated,'nLags',nLags,'nSides',2,'domainIncr',0.001);

        end

        figure(figNum)
        figNum = figNum+1;
        [R, V, yp] = nlid_resid(SRS_inverse,Zcur_simulated);
        
        accuracy_identification = [accuracy_identification V];
        
        %% Inverse SRS Validation
        
        control_system_test_output_double = [];
        inverse_SRS_validation_accuracy = [];

        for signal = 1:num_signals-1

            control_system_test_output = nlsim(SRS_inverse,validation_output(:,signal));
            set(control_system_test_output,'domainIncr',0.001);
            control_system_test_output_double(:,signal) = control_system_test_output.dataSet;

            inverse_SRS_validation_accuracy(signal,:) = vaf(simulated_input(1999:end-1000,signal+1),control_system_test_output_double(1000:end-1000,signal));

        end
        
        accuracy_validation = [accuracy_validation inverse_SRS_validation_accuracy(1,:)];

        nLags_plot = [nLags_plot nLags];
        nLags = nLags - 100;
        
    end
    
    %% Plot the Results
    figure(figNum)
    figNum = figNum+1;
    hold on
    plot(nLags_plot,accuracy_identification,'Linewidth',2)
    plot(nLags_plot,accuracy_validation,'Linewidth',2)
    hold off
    title('Accuracy vs NLags','Fontsize',20)
    xlabel('NLags', 'Fontsize',18)
    ylabel('Accuracy (%VAF)','Fontsize',18)
    legend('Identification','Validation','Fontsize',16)
    grid on
    
elseif multi_model == true
    
    % Test the signals on all the models
    simulated_output_ident = [];
    simulated_outputs_ident = [];

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

        simulated_output_ident = nlsim(SRS,simulated_input_ident);
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
    
    %% Display Mean and Std of Model Accuracies
    
    %Identification Accuracy
    accuracy_identification_LNL_mean = mean(accuracy_identification_LNL);
    accuracy_identification_LNL_std = std(accuracy_identification_LNL);

    accuracy_identification_Hamm_mean = mean(accuracy_identification_Hamm);
    accuracy_identification_Hamm_std = std(accuracy_identification_Hamm);

    accuracy_identification_Weiner_mean = mean(accuracy_identification_Weiner);
    accuracy_identification_Weiner_std = std(accuracy_identification_Weiner);

    accuracy_identification_IRF_mean = mean(accuracy_identification_IRF);
    accuracy_identification_IRF_std = std(accuracy_identification_IRF);
    
    disp(['Identification Accuracy Mean for LNL Inverse SRS: ' num2str(round(accuracy_identification_LNL_mean,1)) '%'])
    disp(['Identification Accuracy Std for LNL Inverse SRS: ' num2str(round(accuracy_identification_LNL_std,1)) '%'])
    
    disp(['Identification Accuracy Mean for Hammerstein Inverse SRS: ' num2str(round(accuracy_identification_Hamm_mean,1)) '%'])
    disp(['Identification Accuracy Std for Hammerstein Inverse SRS: ' num2str(round(accuracy_identification_Hamm_std,1)) '%'])
    
    disp(['Identification Accuracy Mean for Wiener Inverse SRS: ' num2str(round(accuracy_identification_Weiner_mean,1)) '%'])
    disp(['Identification Accuracy Std for Wiener Inverse SRS: ' num2str(round(accuracy_identification_Weiner_std,1)) '%'])
    
    disp(['Identification Accuracy Mean for IRF Inverse SRS: ' num2str(round(accuracy_identification_IRF_mean,1)) '%'])
    disp(['Identification Accuracy Std for IRF Inverse SRS: ' num2str(round(accuracy_identification_IRF_std,1)) '%'])

    %Validation Accuracy
    accuracy_validation_LNL_mean = mean(accuracy_validation_LNL);
    accuracy_validation_LNL_std = std(accuracy_validation_LNL);

    accuracy_validation_Hamm_mean = mean(accuracy_validation_Hamm);
    accuracy_validation_Hamm_std = std(accuracy_validation_Hamm);

    accuracy_validation_Weiner_mean = mean(accuracy_validation_Weiner);
    accuracy_validation_Weiner_std = std(accuracy_validation_Weiner);

    accuracy_validation_IRF_mean = mean(accuracy_validation_IRF);
    accuracy_validation_IRF_std = std(accuracy_validation_IRF);
    
    disp(['Validation Accuracy Mean for LNL Inverse SRS: ' num2str(round(accuracy_validation_LNL_mean,1)) '%'])
    disp(['Validation Accuracy Std for LNL Inverse SRS: ' num2str(round(accuracy_validation_LNL_std,1)) '%'])
    
    disp(['Validation Accuracy Mean for Hammerstein Inverse SRS: ' num2str(round(accuracy_validation_Hamm_mean,1)) '%'])
    disp(['Validation Accuracy Std for Hammerstein Inverse SRS: ' num2str(round(accuracy_validation_Hamm_std,1)) '%'])
    
    disp(['Validation Accuracy Mean for Wiener Inverse SRS: ' num2str(round(accuracy_validation_Weiner_mean,1)) '%'])
    disp(['Validation Accuracy Std for Wiener Inverse SRS: ' num2str(round(accuracy_validation_Weiner_std,1)) '%'])
    
    disp(['Validation Accuracy Mean for IRF Inverse SRS: ' num2str(round(accuracy_validation_IRF_mean,1)) '%'])
    disp(['Validation Accuracy Std for IRF Inverse SRS: ' num2str(round(accuracy_validation_IRF_std,1)) '%'])
    
end



