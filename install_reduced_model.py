import os

installationDirectory = os.path.abspath(os.path.dirname(__file__))

header = '#!/bin/bash\n'
header += '\n'
header += 'cd '
header += os.path.join(installationDirectory, 'evoked_activity/reduced_model')
header += '\n'

evokedSynapseNumbers = range(25, 325, 25)
evokedSynapseOffset = 8.0
whiskers = ['PW', 'SW', 'E2']
evokedSynapseTiming = [(1.0, float(i)) for i in range(2,22)]

scriptFolder = os.path.join(installationDirectory, 'reduced_model_scripts')
if not os.path.exists(scriptFolder):
    os.makedirs(scriptFolder)

for synNumber in evokedSynapseNumbers:
    for synTiming in evokedSynapseTiming:
        for whisker in whiskers:
            suffix = 'reduced_model_' + whisker + '_nrSyn_' + str(synNumber) + '_synPeakTime_' + str(synTiming[0]) + '_synMedTime_' + str(synTiming[1]) + '.sh'
            scriptName = os.path.join(scriptFolder, suffix)
            with open(scriptName, 'w') as scriptFile:
                scriptFile.write(header)
                line = 'python L5PT_reduced_model.py '
                
                # output name
                line += 'reduced_model_'
                line += whisker
                line += '/nrSyn_'
                line += str(synNumber)
                line += '_synPeakTime_'
                line += str(synTiming[0])
                line += '_synMedTime_'
                line += str(synTiming[1])
                
                # neuron model parameters
                line += ' ../../connectome/control/C2center/86_CDK_20041214_BAC_run5_soma_Hay2013_C2center_apic_rec_scaled_diameters.param '
                
                # network model parameters
                line += 'C2_evoked_UpState_EXC_generic_INH_generic_C2center.param '
                
                # generic EXC/INH synapse parameters
                line += whisker
                line += '_evoked_UpState_EXC_INH_Generic_C2center.param '
                
                # subcellular distribution evoked synapses
                line += '../../connectome/reduced_model/'
                line += whisker
                line += '_subcellular_distribution.csv '
                
                # subcellular distribution evoked INH synapses
                line += '../../connectome/reduced_model/'
                line += whisker
                line += '_subcellular_distribution_INH.csv '
                
                # parameters
                line += str(synNumber)
                line += ' '
                line += str(evokedSynapseOffset)
                line += ' '
                line += str(synTiming[0])
                line += ' '
                line += str(synTiming[1])
                line += '\n'
                
                scriptFile.write(line)
            
            # make sure file is also executable
            os.chmod(scriptName, 0777)

