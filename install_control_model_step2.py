import os, os.path

installationDirectory = os.path.abspath(os.path.dirname(__file__))

whiskers = ['B1','B2','B3','C1','C2','C3','D1','D2','D3','E2']

cellLocations = ['B1border','B2border','B3border','C1border','C2center',\
                'C3border','D1border','D2border','D3border']

pythonScriptName = os.path.join(installationDirectory, 'evoked_network_param_from_template.py')
locationBaseName = os.path.join(installationDirectory, 'connectome', 'control')
locationFolderNames = {'B1border': '86_L5_CDK20041214_nr3L5B_dend_PC_neuron_transform_registered_C2_B1border_synapses_20150504-1602_10393',\
                        'B2border': '86_L5_CDK20041214_nr3L5B_dend_PC_neuron_transform_registered_C2_B2border_synapses_20150504-2001_11709',\
                        'B3border': '86_L5_CDK20041214_nr3L5B_dend_PC_neuron_transform_registered_C2_B3border_synapses_20150504-1959_11710',\
                        'C1border': '86_L5_CDK20041214_nr3L5B_dend_PC_neuron_transform_registered_C2_C1border_synapses_20150504-1606_10377',\
                        'C2center': '86_L5_CDK20041214_nr3L5B_dend_PC_neuron_transform_registered_C2center_synapses_20150504-1611_10389',\
                        'C3border': '86_L5_CDK20041214_nr3L5B_dend_PC_neuron_transform_registered_C2_C3border_synapses_20150504-1602_10391',\
                        'D1border': '86_L5_CDK20041214_nr3L5B_dend_PC_neuron_transform_registered_C2_D1border_synapses_20150504-1612_10395',\
                        'D2border': '86_L5_CDK20041214_nr3L5B_dend_PC_neuron_transform_registered_C2_D2border_synapses_20150504-1629_10397',\
                        'D3border': '86_L5_CDK20041214_nr3L5B_dend_PC_neuron_transform_registered_C2_D3border_synapses_20150504-1630_10399'}

synParamName = os.path.join(installationDirectory, 'celltype_PSTH/ongoing_activity_celltype_template_exc_conductances_fitted.param')
outPath = os.path.join(installationDirectory, 'evoked_activity/control/')

fname = os.path.join(installationDirectory, 'install_control_model_step3.sh')
outName = 'deflection_control_cell_location'

with open(fname, 'w') as scriptFile:
    header = '#!/bin/bash\n\n'
    scriptFile.write(header)
    for cellLocation in cellLocations:
        for whisker in whiskers:
            line = 'python '
            line += pythonScriptName
            line += ' '
            line += synParamName
            line += ' '
            line += os.path.join(locationBaseName, cellLocation, locationFolderNames[cellLocation], 'NumberOfConnectedCells.csv')
            line += ' '
            line += os.path.join(locationBaseName, cellLocation, locationFolderNames[cellLocation], locationFolderNames[cellLocation])
            line += '.syn'
            line += ' '
            line += whisker
            # output name
            line += ' '
            line += outPath
            line += whisker
            line += '_'
            line += outName
            line += '_'
            line += cellLocation
            line += '.param\n'
            scriptFile.write(line)

# make sure file is also executable
os.chmod(fname, 0777)
