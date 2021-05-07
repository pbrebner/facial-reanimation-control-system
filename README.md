# Facial Reanimation Control System
Repository for facial reanimation control system project, a thesis project for Masters of Biological and Biomedical Engineering at McGill University.

In order to run the project code, you need to download all files of the REKLAB NLID repository from: https://github.com/reklab/reklab_public and add it to your MATLAB path.

The project is split into three sections:
1. Identifying and evaluating the EMG Response System (ERS)
2. Identifying and evaluating the Stimulus Response System (SRS)
3. Combining the ERS and the inverse of the SRS to create the Facial Reanimation Control System (FRCS)

## ERS Model
The ERS represents the healthy side of the face and models healthy side EMG input to a healthy displacement output.

Matlab/Simulink Files:
1. ERS_simulation: The Simulink model of the ERS simulation used to generate the simulated data
2. ERS_model: Identifies the ERS using the simulated data from the ERS_simulation
3. ERS_model_test: Validates the models identified in the ERS_model
4. ERS_model_test_freq: Evaluates identified ERS models based on the number of movement pulses in the identification signal
5. ERS_model_test_amp: Evaluates the identified ERS models based on either the amplitudem of the identification signal or the record length of the identification signal
6. ERS_model_test_acc_vs_noise: Evaluates the identified ERS models based on the amount of output noise added to the ERS_simulation

## SRS Model
The SRS represents the paralyzed side of the face and models stimulus amplitude modulation input to a paralyzed displacement output.

Matlab/Simulink Files:
    Same as the ERS Model section but with the SRS versions

## FRCS
The Facial Reanimation Control System (FRCS) combines the ERS with the inverse of the SRS to model the healthy side displacement to a stimulus amplitude modulation output required to duplicate the healthy movement on the paralyzed side of the face

Matlab/Simulink Files:
1. FRCS_identify_models: Identifies an ERS and SRS Model using simulated data
2. FRCS_identify_inverse_SRS: takes the SRS model and estimates the inverse SRS model
3. FRCS_identify_inverse_SRS_test: Evaluates the inverse SRS performance depending on the model structure and number of lags of inverse SRS IRF(s)
4. FRCS_execute: Combines the ERS and inverse SRS to create the FRCS. An EMG input is passed into the FRCS and the outputs are evaluated.
