import os, glob
from sumatra.parameters import build_parameters, NTParameterSet

installationDirectory = os.path.abspath(os.path.dirname(__file__))

controlFolder = os.path.join(installationDirectory, 'evoked_activity/control')
manipulation2Folder = os.path.join(installationDirectory, 'evoked_activity/manipulation2')
if not os.path.exists(manipulation2Folder):
    os.makedirs(manipulation2Folder)

for fname in glob.glob(os.path.join(controlFolder, '*')):
    if fname.endswith('.param'):
        outName = fname.replace('control', 'manipulation2')
        # load control file as ParameterSet
        controlParameters = build_parameters(fname)
        # pop population (manipulation 2: other EXC cell types in E2)
        newParameters = controlParameters.tree_copy()
        otherExcTypesE2 = ['L2_E2','L34_E2','L4py_E2','L4sp_E2','L4ss_E2','L5st_E2','L5tt_E2','L6ccinv_E2','L6ct_E2']
        for population in newParameters.network.keys():
            if population in otherExcTypesE2:
                newParameters.network.pop(population)
        # save ParameterSet as new file
        newParameters.save(outName)
    else:
        continue

header = '#!/bin/bash\n'
header += '\n'
header += 'cd '
header += manipulation2Folder
header += '\n'

whiskers = ['B1','B2','B3','C1','C2','C3','D1','D2','D3','E2']
cellLocations = ['B1border','B2border','B3border','C1border','C2center',\
                'C3border','D1border','D2border','D3border']

scriptFolder = os.path.join(installationDirectory, 'manipulation2_scripts')
if not os.path.exists(scriptFolder):
    os.makedirs(scriptFolder)

for whisker in whiskers:
    for cellLocation in cellLocations:
        suffix = whisker + '_deflection_manipulation2_cell_location_' + cellLocation + '.sh'
        scriptName = os.path.join(scriptFolder, suffix)
        with open(scriptName, 'w') as scriptFile:
            scriptFile.write(header)
            line = 'python '
            line += os.path.join(controlFolder, 'L5PT_control.py')
            line += ' '

            # output name
            line += whisker
            line += '_deflection_manipulation2_cell_location_'
            line += cellLocation
            
            # neuron model parameters
            line += ' ../../connectome/control/'
            line += cellLocation
            line += '/86_CDK_20041214_BAC_run5_soma_Hay2013_'
            line += cellLocation
            line += '_apic_rec_scaled_diameters.param '
            
            # network model parameters
            line += whisker
            line += '_deflection_manipulation2_cell_location_'
            line += cellLocation
            line += '.param\n'
            
            scriptFile.write(line)
        
        # make sure file is also executable
        os.chmod(scriptName, 0777)
