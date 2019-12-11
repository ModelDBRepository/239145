import os

installationDirectory = os.path.abspath(os.path.dirname(__file__))

header = '#!/bin/bash\n'
header += '\n'
header += 'cd '
header += os.path.join(installationDirectory, 'evoked_activity/control')
header += '\n'

whiskers = ['B1','B2','B3','C1','C2','C3','D1','D2','D3','E2']
cellLocations = ['B1border','B2border','B3border','C1border','C2center',\
                'C3border','D1border','D2border','D3border']
                
scriptFolder = os.path.join(installationDirectory, 'control_scripts')
if not os.path.exists(scriptFolder):
    os.makedirs(scriptFolder)

for whisker in whiskers:
    for cellLocation in cellLocations:
        suffix = whisker + '_deflection_control_cell_location_' + cellLocation + '.sh'
        scriptName = os.path.join(scriptFolder, suffix)
        with open(scriptName, 'w') as scriptFile:
            scriptFile.write(header)
            line = 'python L5PT_control.py '

            # output name
            line += whisker
            line += '_deflection_control_cell_location_'
            line += cellLocation
            
            # neuron model parameters
            line += ' ../../connectome/control/'
            line += cellLocation
            line += '/86_CDK_20041214_BAC_run5_soma_Hay2013_'
            line += cellLocation
            line += '_apic_rec_scaled_diameters.param '
            
            # network model parameters
            line += whisker
            line += '_deflection_control_cell_location_'
            line += cellLocation
            line += '.param\n'
            
            scriptFile.write(line)
        
        # make sure file is also executable
        os.chmod(scriptName, 0777)
