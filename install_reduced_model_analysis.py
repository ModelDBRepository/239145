import os

installationDirectory = os.path.abspath(os.path.dirname(__file__))

resultsFolder = os.path.join(installationDirectory, 'evoked_activity/reduced_model/results')
spikeTimesScript = os.path.join(installationDirectory, 'lib/visualization/spike_raster_plots.py')
probabilityScript = os.path.join(installationDirectory, 'lib/visualization/reduced_model_iso-probability_contours.py')

analysisFolder = os.path.join(installationDirectory, 'reduced_model_analysis')
spikeRasterFolder = os.path.join(analysisFolder, 'spike_raster_plots')
if not os.path.exists(analysisFolder):
    os.makedirs(analysisFolder)
if not os.path.exists(spikeRasterFolder):
    os.makedirs(os.path.join(spikeRasterFolder, 'PW'))
    os.makedirs(os.path.join(spikeRasterFolder, 'SW'))
    os.makedirs(os.path.join(spikeRasterFolder, 'E2'))

scriptName = os.path.join(analysisFolder, 'reduced_model_analysis_script.sh')
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
    line += os.path.join(resultsFolder, 'reduced_model_PW')
    line += ' '
    line += os.path.join(spikeRasterFolder, 'PW')
    line += '\n'
    line += 'python '
    line += spikeTimesScript
    line += ' '
    line += os.path.join(resultsFolder, 'reduced_model_SW')
    line += ' '
    line += os.path.join(spikeRasterFolder, 'SW')
    line += '\n'
    line += 'python '
    line += spikeTimesScript
    line += ' '
    line += os.path.join(resultsFolder, 'reduced_model_E2')
    line += ' '
    line += os.path.join(spikeRasterFolder, 'E2')
    line += '\n'
    line += '\n\n'
    scriptFile.write(line)
    
    # generate iso-probability contour plots
    line = 'echo \"******************************************\"\n'
    line += 'echo \"ANALYSIS STEP 2: generate iso-probability contour plots\"\n'
    line += 'echo \"******************************************\"\n'
    line += 'python '
    line += probabilityScript
    line += ' '
    line += os.path.join(resultsFolder, 'reduced_model_PW')
    line += ' spike_times.csv PW '
    line += os.path.join(analysisFolder, 'Iso-probability_PW')
    line += '\n'
    line += 'python '
    line += probabilityScript
    line += ' '
    line += os.path.join(resultsFolder, 'reduced_model_SW')
    line += ' spike_times.csv SW '
    line += os.path.join(analysisFolder, 'Iso-probability_SW')
    line += '\n'
    line += 'python '
    line += probabilityScript
    line += ' '
    line += os.path.join(resultsFolder, 'reduced_model_E2')
    line += ' spike_times.csv E2 '
    line += os.path.join(analysisFolder, 'Iso-probability_E2')
    line += '\n\n'
    scriptFile.write(line)
    
os.chmod(scriptName, 0777)
