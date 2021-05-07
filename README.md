# Facial Reanimation Control System
Repository for facial reanimation control system, a thesis project for Masters of Biological and Biomedical Engineering at McGill University.

In order to run the project code, you need to download all files of the REKLAB NLID repository from: https://github.com/reklab/reklab_public and add it to your MATLAB path.

## Abstract
Facial palsy can be a devastating condition that often impairs blinking and can lead to other functional, esthetic, communication, and psychological issues. The current reanimation technologies aimed at restoring function often only achieve limited results, usually only enabling functional movement through non-spontaneous means. Due to recent advances, techniques like facial pacing have become feasible. Facial pacing is a technique proposed to restore the symmetry of facial movement that has been lost due to unilateral facial paralysis. The main idea of facial pacing is to measure the muscle activity from the healthy side of the face with electromyography (EMG) and use this to activate the muscles on the paralyzed side through functional electrical stimulation (FES). This thesis uses simulation studies to investigate the feasibility of building a facial reanimation control system (FRCS) that uses EMG to determine the stimulus required for synchronous movements on the paralyzed side. Based on our knowledge of muscle physiology, simulation models were developed that represent the healthy and paralyzed sides of the face. These models simulated the signals needed to identify the models required for the FRCS. The identified models were validated and tested for various conditions to determine if they were robust enough for use in a facial reanimation device. The results indicate that we can successfully identify these models, using the simulated data, that are able to achieve high accuracies for reasonable signal record lengths and output noise levels. The results are a proof of concept that shows that the FRCS is a feasible approach to facial pacing. The next step will be to apply the methods developed in this thesis to experimental data.

## Project Code
The project is split into three sections:
1. Identifying and evaluating the EMG Response System (ERS)
2. Identifying and evaluating the Stimulus Response System (SRS)
3. Combining the ERS and the inverse of the SRS to create the Facial Reanimation Control System (FRCS)

### ERS Model
The ERS represents the healthy side of the face and models healthy side EMG input to a healthy displacement output.

Matlab/Simulink Files:
1. ERS_simulation: The Simulink model of the ERS simulation used to generate the simulated EMG and healthy displacement data
2. ERS_model: Identifies the ERS using the simulated data from the ERS_simulation
3. ERS_model_test: Validates the models identified in the ERS_model
4. ERS_model_test_freq: Evaluates identified ERS models based on the number of movement pulses in the identification signal
5. ERS_model_test_amp: Evaluates the identified ERS models based on either the amplitudem of the identification signal or the record length of the identification signal
6. ERS_model_test_acc_vs_noise: Evaluates the identified ERS models based on the amount of output noise added to the ERS_simulation

### SRS Model
The SRS represents the paralyzed side of the face and models stimulus amplitude modulation input to a paralyzed displacement output.

Matlab/Simulink Files:
    Same as the ERS Model section but with the SRS versions

### FRCS
The Facial Reanimation Control System (FRCS) combines the ERS with the inverse of the SRS to model the healthy side displacement to a stimulus amplitude modulation output required to duplicate the healthy movement on the paralyzed side of the face

Matlab/Simulink Files:
1. FRCS_identify_models: Identifies an ERS and SRS Model using simulated data
2. FRCS_identify_inverse_SRS: takes the SRS model and estimates the inverse SRS model
3. FRCS_identify_inverse_SRS_test: Evaluates the inverse SRS performance depending on the model structure and number of lags of inverse SRS IRF(s)
4. FRCS_execute: Combines the ERS and inverse SRS to create the FRCS. An EMG input is passed into the FRCS and the outputs are evaluated.
