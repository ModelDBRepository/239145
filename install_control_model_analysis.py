import os

installationDirectory = os.path.abspath(os.path.dirname(__file__))

resultsFolder = os.path.join(installationDirectory, 'evoked_activity/control/results')
spikeTimesScript = os.path.join(installationDirectory, 'lib/visualization/spike_raster_plots.py')
synapseHistogramScript = os.path.join(installationDirectory, 'lib/visualization/active_synapse_histograms.py')
spatiotemporalScript = os.path.join(installationDirectory, 'lib/visualization/spatiotemporal_synaptic_input_pca.py')

analysisFolder = os.path.join(installationDirectory, 'control_analysis')
if not os.path.exists(analysisFolder):
    os.makedirs(analysisFolder)
    os.makedirs(os.path.join(analysisFolder, 'spike_raster_plots'))
    os.makedirs(os.path.join(analysisFolder, 'active_synapses'))
    os.makedirs(os.path.join(analysisFolder, 'spatiotemporal_synaptic_input_PCA'))

whiskers = ('B1', 'B2', 'B3', 'C1', 'C2', 'C3', 'D1', 'D2', 'D3', 'E2')

scriptName = os.path.join(analysisFolder, 'control_analysis_script.sh')
with open(scriptName, 'w') as scriptFile:
    header = '#!/bin/bash\n'
    header += '\n'
    scriptFile.write(header)
    
    # generate spike time files
    line = 'echo \"******************************************\"\n'
    line += 'echo \"ANALYSIS STEP 1: generate spike time files\"\n'
    line += 'echo \"******************************************\"\n'
    line += 'python '
    line += spikeTimesScript
    line += ' '
    line += resultsFolder
    line += ' '
    line += os.path.join(analysisFolder, 'spike_raster_plots')
    line += '\n\n'
    scriptFile.write(line)
    
    # generate spike raster plots/PSTHs
    line = 'echo \"******************************************\"\n'
    line += 'echo \"ANALYSIS STEP 2: generate spike raster plots/PSTHs\"\n'
    line += 'echo \"******************************************\"\n'
    line += 'python '
    line += spikeTimesScript
    line += ' '
    line += resultsFolder
    line += ' vm_all_traces.csv'
    line += ' '
    line += os.path.join(analysisFolder, 'spike_raster_plots')
    line += '\n'

    for whisker in whiskers:
        line += 'python '
        line += spikeTimesScript
        line += ' '
        line += resultsFolder
        line += ' vm_all_traces.csv'
        line += ' '
        line += whisker + '_deflection '
        line += os.path.join(analysisFolder, 'spike_raster_plots', 'summary', whisker + 'deflection_location_all')
        line += '\n'
    line += '\n'
    scriptFile.write(line)
    
    # generate synapse activation histograms
    line = 'echo \"******************************************\"\n'
    line += 'echo \"ANALYSIS STEP 3: generate active synapse histograms\"\n'
    line += 'echo \"******************************************\"\n'
    for whisker in whiskers:
        line += 'python '
        line += synapseHistogramScript
        line += ' '
        line += resultsFolder
        line += ' synapses.csv '
        line += whisker
        line += '_deflection '
        suffix = whisker + '_deflection'
        line += os.path.join(analysisFolder, 'active_synapses', suffix)
        line += '\n'
    line += '\n'
    for whisker in whiskers:
        line += 'python '
        line += synapseHistogramScript
        line += ' '
        line += resultsFolder
        line += ' synapses.csv spike_times.csv '
        line += whisker
        line += '_deflection '
        suffix = whisker + '_deflection'
        line += os.path.join(analysisFolder, 'active_synapses', suffix)
        line += '\n'
    line += '\n'
    scriptFile.write(line)
    
    # perform spatiotemporal synaptic input PCA
    line = 'echo \"******************************************\"\n'
    line += 'echo \"ANALYSIS STEP 4: perform spatiotemporal synaptic input PCA\"\n'
    line += 'echo \"******************************************\"\n'
    for whisker in whiskers:
        line += 'python '
        line += spatiotemporalScript
        line += ' '
        line += resultsFolder
        line += ' '
        line += whisker
        line += '_deflection'
        line += ' spike_times.csv synapses.csv '
        suffix = whisker + '_deflection'
        line += os.path.join(analysisFolder, 'spatiotemporal_synaptic_input_PCA', suffix)
        line += '\n'
    line += '\n'
    scriptFile.write(line)
    
os.chmod(scriptName, 0777)
