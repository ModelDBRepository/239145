'''
Created on Jan 28, 2013

ongoing activity L2 neuron model

@author: robert
'''

import sys
import os, os.path
import glob
import single_cell_analyzer as sca
import single_cell_parser as scp
import numpy as np
import matplotlib.pyplot as plt

def create_active_synapse_histogram(folder, suffix, contains, outName):
    '''
    load all traces, compute spike times
    and create raster plots
    '''
    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    inhTypes = ('L1','L23Trans','L45Sym','L45Peak','L56Trans',\
                'SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6')
    plotTypes = ('INH', 'L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    columns = ('A1','A2','A3','A4','Alpha','B1','B2','B3','B4','Beta',\
               'C1','C2','C3','C4','D1','D2','D3','D4','Delta','E1','E2','E3','E4','Gamma')
    
    fnames = []
#    scan_directory(folder, fnames, suffix)
    scan_directory2(folder, fnames, suffix, contains)
    nrOfFiles = len(fnames)
    print 'Creating active synapse plots from %d files' % nrOfFiles
    
    synapseTimes = {}
    synapseTimesProximal = {}
    synapseTimesDistal = {}
    for col in columns:
        synapseTimes[col] = {}
        synapseTimesProximal[col] = {}
        synapseTimesDistal[col] = {}
        for excType in excTypes:
            synapseTimes[col][excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimesProximal[col][excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimesDistal[col][excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        synapseTimes[col]['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        synapseTimes[col]['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
        synapseTimesProximal[col]['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        synapseTimesProximal[col]['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
        synapseTimesDistal[col]['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        synapseTimesDistal[col]['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
    
    for n in range(len(fnames)):
        fname = fnames[n]
        print 'Loading synapse activation times from file %s (file %d of %d)\r' % (fname,n+1,nrOfFiles) ,
        sys.stdout.flush()
        activeSyns = scp.read_complete_synapse_activation_file(fname)
        
        for col in columns:
            for excType in excTypes:
                synapseTimes[col][excType]['ApicalDendrite'].append([])
                synapseTimes[col][excType]['Dendrite'].append([])
                synapseTimes[col][excType]['Total'].append([])
                synapseTimesProximal[col][excType]['ApicalDendrite'].append([])
                synapseTimesProximal[col][excType]['Dendrite'].append([])
                synapseTimesProximal[col][excType]['Total'].append([])
                synapseTimesDistal[col][excType]['ApicalDendrite'].append([])
                synapseTimesDistal[col][excType]['Dendrite'].append([])
                synapseTimesDistal[col][excType]['Total'].append([])
            synapseTimes[col]['EXC']['ApicalDendrite'].append([])
            synapseTimes[col]['EXC']['Dendrite'].append([])
            synapseTimes[col]['EXC']['Total'].append([])
            synapseTimesProximal[col]['EXC']['ApicalDendrite'].append([])
            synapseTimesProximal[col]['EXC']['Dendrite'].append([])
            synapseTimesProximal[col]['EXC']['Total'].append([])
            synapseTimesDistal[col]['EXC']['ApicalDendrite'].append([])
            synapseTimesDistal[col]['EXC']['Dendrite'].append([])
            synapseTimesDistal[col]['EXC']['Total'].append([])
            synapseTimes[col]['INH']['ApicalDendrite'].append([])
            synapseTimes[col]['INH']['Dendrite'].append([])
            synapseTimes[col]['INH']['Soma'].append([])
            synapseTimes[col]['INH']['Total'].append([])
            synapseTimesProximal[col]['INH']['ApicalDendrite'].append([])
            synapseTimesProximal[col]['INH']['Dendrite'].append([])
            synapseTimesProximal[col]['INH']['Soma'].append([])
            synapseTimesProximal[col]['INH']['Total'].append([])
            synapseTimesDistal[col]['INH']['ApicalDendrite'].append([])
            synapseTimesDistal[col]['INH']['Dendrite'].append([])
            synapseTimesDistal[col]['INH']['Soma'].append([])
            synapseTimesDistal[col]['INH']['Total'].append([])
        
        for synType in activeSyns.keys():
            preCellType = synType.split('_')[0]
            preCol = synType.split('_')[1]
            for col in columns:
                if col == preCol:
                    for excType in excTypes:
                        if excType == preCellType:
                            for syn in activeSyns[synType]:
                                somaDist = syn[1]
                                structure = syn[4]
                                synTimes = syn[5]
                                synapseTimes[col][excType][structure][n].extend(synTimes)
                                synapseTimes[col][excType]['Total'][n].extend(synTimes)
                                synapseTimes[col]['EXC'][structure][n].extend(synTimes)
                                synapseTimes[col]['EXC']['Total'][n].extend(synTimes)
                                if somaDist < 500.0:
                                    synapseTimesProximal[col][excType][structure][n].extend(synTimes)
                                    synapseTimesProximal[col][excType]['Total'][n].extend(synTimes)
                                    synapseTimesProximal[col]['EXC'][structure][n].extend(synTimes)
                                    synapseTimesProximal[col]['EXC']['Total'][n].extend(synTimes)
                                if somaDist >= 500.0:
                                    synapseTimesDistal[col][excType][structure][n].extend(synTimes)
                                    synapseTimesDistal[col][excType]['Total'][n].extend(synTimes)
                                    synapseTimesDistal[col]['EXC'][structure][n].extend(synTimes)
                                    synapseTimesDistal[col]['EXC']['Total'][n].extend(synTimes)
                    for inhType in inhTypes:
                        if inhType == preCellType:
                            for syn in activeSyns[synType]:
                                somaDist = syn[1]
                                structure = syn[4]
                                synTimes = syn[5]
                                synapseTimes[col]['INH'][structure][n].extend(synTimes)
                                synapseTimes[col]['INH']['Total'][n].extend(synTimes)
                                if somaDist < 500.0:
                                    synapseTimesProximal[col]['INH'][structure][n].extend(synTimes)
                                    synapseTimesProximal[col]['INH']['Total'][n].extend(synTimes)
                                if somaDist >= 500.0:
                                    synapseTimesDistal[col]['INH'][structure][n].extend(synTimes)
                                    synapseTimesDistal[col]['INH']['Total'][n].extend(synTimes)
    
    print ''
    tOffset = 100.0
    tOngoing = 120.0
    tOngoingWindow = 100.0
    tStim = 245.0
    tStimWindow = 25.0
    binWidth = 1.0
    maxCount = 0
    synapseHistogramsEvoked = {}
    synapseHistogramsOngoing = {}
    synapseHistogramsEvokedProximal = {}
    synapseHistogramsOngoingProximal = {}
    synapseHistogramsEvokedDistal = {}
    synapseHistogramsOngoingDistal = {}
    for col in columns:
        synapseHistogramsEvoked[col] = {}
        synapseHistogramsOngoing[col] = {}
        synapseHistogramsEvokedProximal[col] = {}
        synapseHistogramsOngoingProximal[col] = {}
        synapseHistogramsEvokedDistal[col] = {}
        synapseHistogramsOngoingDistal[col] = {}
        for cellType in synapseTimes[col].keys():
            synapseHistogramsEvoked[col][cellType] = {}
            synapseHistogramsOngoing[col][cellType] = {}
            synapseHistogramsEvokedProximal[col][cellType] = {}
            synapseHistogramsOngoingProximal[col][cellType] = {}
            synapseHistogramsEvokedDistal[col][cellType] = {}
            synapseHistogramsOngoingDistal[col][cellType] = {}
            for structure in synapseTimes[col][cellType].keys():
                synTimes1 = synapseTimes[col][cellType][structure]
                hist1, bins1 = sca.PSTH_from_spike_times(synTimes1, binWidth, tStim, tStim+tStimWindow)
                synapseHistogramsEvoked[col][cellType][structure] = np.sum(hist1)
                hist2, bins2 = sca.PSTH_from_spike_times(synTimes1, binWidth, tOngoing, tOngoing+tOngoingWindow)
                synapseHistogramsOngoing[col][cellType][structure] = np.sum(hist2)
                
                synTimes2 = synapseTimesProximal[col][cellType][structure]
                hist3, bins3 = sca.PSTH_from_spike_times(synTimes2, binWidth, tStim, tStim+tStimWindow)
                synapseHistogramsEvokedProximal[col][cellType][structure] = np.sum(hist3)
                hist4, bins4 = sca.PSTH_from_spike_times(synTimes2, binWidth, tOngoing, tOngoing+tOngoingWindow)
                synapseHistogramsOngoingProximal[col][cellType][structure] = np.sum(hist4)
                
                synTimes3 = synapseTimesDistal[col][cellType][structure]
                hist5, bins5 = sca.PSTH_from_spike_times(synTimes3, binWidth, tStim, tStim+tStimWindow)
                synapseHistogramsEvokedDistal[col][cellType][structure] = np.sum(hist5)
                hist6, bins6 = sca.PSTH_from_spike_times(synTimes3, binWidth, tOngoing, tOngoing+tOngoingWindow)
                synapseHistogramsOngoingDistal[col][cellType][structure] = np.sum(hist6)
    
    ongoingOutName = outName + '_ongoing_syn.csv'
    with open(ongoingOutName, 'w') as outputTable:
        header = 'Column\tcell type\tnr. of active synapses (100ms)\n'
        outputTable.write(header)
        for cellType in plotTypes:
            for col in columns:
                line = col
                line += '\t'
                line += cellType
                line += '\t'
                line += str(synapseHistogramsOngoing[col][cellType]['Total'])
                line += '\n'
                outputTable.write(line)
    
    writeHeader = True
    evokedOutName = outName + '_total_evoked_syn.csv'
    with open(evokedOutName, 'w') as outputTable:
        header = 'Column\tcell type\tnr. of active synapses (25ms post-stimulus)\n'
        outputTable.write(header)
        for cellType in plotTypes:
            for col in columns:
                line = col
                line += '\t'
                line += cellType
                line += '\t'
                line += str(synapseHistogramsEvoked[col][cellType]['Total'])
                line += '\n'
                outputTable.write(line)
    
    ongoingProximalOutName = outName + '_ongoing_proximal_syn.csv'
    with open(ongoingProximalOutName, 'w') as outputTable:
        header = 'Column\tcell type\tnr. of active synapses (100ms)\n'
        outputTable.write(header)
        for cellType in plotTypes:
            for col in columns:
                line = col
                line += '\t'
                line += cellType
                line += '\t'
                line += str(synapseHistogramsOngoingProximal[col][cellType]['Total'])
                line += '\n'
                outputTable.write(line)
    
    writeHeader = True
    evokedProximalOutName = outName + '_total_evoked_proximal_syn.csv'
    with open(evokedProximalOutName, 'w') as outputTable:
        header = 'Column\tcell type\tnr. of active synapses (25ms post-stimulus)\n'
        outputTable.write(header)
        for cellType in plotTypes:
            for col in columns:
                line = col
                line += '\t'
                line += cellType
                line += '\t'
                line += str(synapseHistogramsEvokedProximal[col][cellType]['Total'])
                line += '\n'
                outputTable.write(line)
    
    ongoingDistalOutName = outName + '_ongoing_distal_syn.csv'
    with open(ongoingDistalOutName, 'w') as outputTable:
        header = 'Column\tcell type\tnr. of active synapses (100ms)\n'
        outputTable.write(header)
        for cellType in plotTypes:
            for col in columns:
                line = col
                line += '\t'
                line += cellType
                line += '\t'
                line += str(synapseHistogramsOngoingDistal[col][cellType]['Total'])
                line += '\n'
                outputTable.write(line)
    
    writeHeader = True
    evokedDistalOutName = outName + '_total_evoked_distal_syn.csv'
    with open(evokedDistalOutName, 'w') as outputTable:
        header = 'Column\tcell type\tnr. of active synapses (25ms post-stimulus)\n'
        outputTable.write(header)
        for cellType in plotTypes:
            for col in columns:
                line = col
                line += '\t'
                line += cellType
                line += '\t'
                line += str(synapseHistogramsEvokedDistal[col][cellType]['Total'])
                line += '\n'
                outputTable.write(line)


def create_active_synapse_histogram_spike_no_spike(folder, suffix, vmSuffix, contains, outName):
    '''
    load all traces, compute spike times
    and create raster plots
    '''
    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    inhTypes = ('L1','L23Trans','L45Sym','L45Peak','L56Trans',\
                'SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6')
    plotTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM', 'EXC', 'INH')
    
    vmNames = []
    scan_directory2(folder, vmNames, vmSuffix, contains)
    synNames = []
    scan_directory2(folder, synNames, suffix, contains)
    nrOfFiles = len(synNames)
    print 'Creating active synapse plots from %d files' % nrOfFiles
    synData = {}
    for fname in synNames:
        activeSyns = scp.read_complete_synapse_activation_file(fname)
        synData[fname] = activeSyns
    
    synapseTimes = {}
    spikeTrialSyns = {}
    noSpikeTrialSyns = {}
    for excType in excTypes:
        synapseTimes[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        spikeTrialSyns[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        noSpikeTrialSyns[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
    synapseTimes['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
    synapseTimes['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
    spikeTrialSyns['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
    spikeTrialSyns['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
    noSpikeTrialSyns['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
    noSpikeTrialSyns['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
    
    tStim = 253.0
    earlyWindow = 17.0
    trials = []
    trialSpikeTimes = [[] for j in range(len(vmNames))]
    trialWithSpikes = {}
    
    for n in range(len(vmNames)):
        fname = vmNames[n]
        print 'Loading spike times from file %s' % fname
        trialWithSpikes[n] = []
        trialSpikeTimes_ = scp.read_spike_times_file(fname)
        nrOfTrials = len(trialSpikeTimes_.keys())
        trials.append(nrOfTrials)
        for trial in trialSpikeTimes_.keys():
            trialSpikeTimes[n].append([])
            trialWithSpikes_ = False
            for tSpike in trialSpikeTimes_[trial]:
                if tSpike >= tStim and tSpike < tStim + earlyWindow:
                    trialSpikeTimes[n][-1].append(tSpike-tStim)
                    trialWithSpikes_ = True
            trialWithSpikes[n].append(trialWithSpikes_)
    
    for n in range(len(vmNames)):
        nrSpikeTrials = 0
        nrNoSpikeTrials = 0
        earlyProxSyns = 0
        for trialNr in range(trials[n]):
            print 'Counting active synapses in trial %d of %d\r' % (trialNr+1, trials[n]),
            sys.stdout.flush()
            synTrialStr = 'simulation_run%04d_synapses.csv' % trialNr
            synTrialFile = ''
            tmpVmName = vmNames[n]
            for name in synNames:
                if synTrialStr in name and os.path.split(tmpVmName)[0] == os.path.split(name)[0]:
                    synTrialFile = name
                    break
            if synTrialFile == '':
                errstr = 'Could not find synapse activation file for trial nr. %d' % trialNr
                raise RuntimeError(errstr)
            activeSyns = synData[synTrialFile]
            
            for excType in excTypes:
                synapseTimes[excType]['ApicalDendrite'].append([])
                synapseTimes[excType]['Dendrite'].append([])
                synapseTimes[excType]['Total'].append([])
                if trialWithSpikes[n][trialNr]:
                    spikeTrialSyns[excType]['ApicalDendrite'].append([])
                    spikeTrialSyns[excType]['Dendrite'].append([])
                    spikeTrialSyns[excType]['Total'].append([])
                else:
                    noSpikeTrialSyns[excType]['ApicalDendrite'].append([])
                    noSpikeTrialSyns[excType]['Dendrite'].append([])
                    noSpikeTrialSyns[excType]['Total'].append([])
            synapseTimes['EXC']['ApicalDendrite'].append([])
            synapseTimes['EXC']['Dendrite'].append([])
            synapseTimes['EXC']['Total'].append([])
            synapseTimes['INH']['ApicalDendrite'].append([])
            synapseTimes['INH']['Dendrite'].append([])
            synapseTimes['INH']['Soma'].append([])
            synapseTimes['INH']['Total'].append([])
            if trialWithSpikes[n][trialNr]:
                spikeTrialSyns['EXC']['ApicalDendrite'].append([])
                spikeTrialSyns['EXC']['Dendrite'].append([])
                spikeTrialSyns['EXC']['Total'].append([])
                spikeTrialSyns['INH']['ApicalDendrite'].append([])
                spikeTrialSyns['INH']['Dendrite'].append([])
                spikeTrialSyns['INH']['Soma'].append([])
                spikeTrialSyns['INH']['Total'].append([])
            else:
                noSpikeTrialSyns['EXC']['ApicalDendrite'].append([])
                noSpikeTrialSyns['EXC']['Dendrite'].append([])
                noSpikeTrialSyns['EXC']['Total'].append([])
                noSpikeTrialSyns['INH']['ApicalDendrite'].append([])
                noSpikeTrialSyns['INH']['Dendrite'].append([])
                noSpikeTrialSyns['INH']['Soma'].append([])
                noSpikeTrialSyns['INH']['Total'].append([])
            
            for synType in activeSyns.keys():
                preCellType = synType.split('_')[0]
                for excType in excTypes:
                    if excType == preCellType:
                        for syn in activeSyns[synType]:
                            somaDist = syn[1]
                            structure = syn[4]
                            synTimes = syn[5]
                            if somaDist < 500.0:
                                synapseTimes[excType][structure][-1].extend(synTimes)
                                synapseTimes[excType]['Total'][-1].extend(synTimes)
                                synapseTimes['EXC'][structure][-1].extend(synTimes)
                                synapseTimes['EXC']['Total'][-1].extend(synTimes)
                                if trialWithSpikes[n][trialNr]:
                                    spikeTrialSyns[excType][structure][-1].extend(synTimes)
                                    spikeTrialSyns[excType]['Total'][-1].extend(synTimes)
                                    spikeTrialSyns['EXC'][structure][-1].extend(synTimes)
                                    spikeTrialSyns['EXC']['Total'][-1].extend(synTimes)
                                else:
                                    noSpikeTrialSyns[excType][structure][-1].extend(synTimes)
                                    noSpikeTrialSyns[excType]['Total'][-1].extend(synTimes)
                                    noSpikeTrialSyns['EXC'][structure][-1].extend(synTimes)
                                    noSpikeTrialSyns['EXC']['Total'][-1].extend(synTimes)
                for inhType in inhTypes:
                    if inhType == preCellType:
                        for syn in activeSyns[synType]:
                            somaDist = syn[1]
                            structure = syn[4]
                            synTimes = syn[5]
                            if somaDist < 500.0:
                                synapseTimes['INH'][structure][-1].extend(synTimes)
                                synapseTimes['INH']['Total'][-1].extend(synTimes)
                                if trialWithSpikes[n][trialNr]:
                                    spikeTrialSyns['INH'][structure][-1].extend(synTimes)
                                    spikeTrialSyns['INH']['Total'][-1].extend(synTimes)
                                else:
                                    noSpikeTrialSyns['INH'][structure][-1].extend(synTimes)
                                    noSpikeTrialSyns['INH']['Total'][-1].extend(synTimes)
        print ''
    
    tOffset = 100.0
    tPlotBegin = 220.0
    tPlotBeginWindow = 50.0
    binWidth = 1.0
    maxCount = 0
    synapseHistograms = {}
    spikeTrialHistograms = {}
    noSpikeTrialHistograms = {}
    for cellType in synapseTimes.keys():
        synapseHistograms[cellType] = {}
        spikeTrialHistograms[cellType] = {}
        noSpikeTrialHistograms[cellType] = {}
        for structure in synapseTimes[cellType].keys():
            synTimes = synapseTimes[cellType][structure]
            hist, bins = sca.PSTH_from_spike_times(synTimes, binWidth, tOffset, tPlotBegin+tPlotBeginWindow)
            synapseHistograms[cellType][structure] = hist, bins
            spikeTrialSynTimes = spikeTrialSyns[cellType][structure]
            hist2, bins2 = sca.PSTH_from_spike_times(spikeTrialSynTimes, binWidth, tOffset, tPlotBegin+tPlotBeginWindow)
            spikeTrialHistograms[cellType][structure] = hist2, bins2
            noSpikeTrialSynTimes = noSpikeTrialSyns[cellType][structure]
            hist3, bins3 = sca.PSTH_from_spike_times(noSpikeTrialSynTimes, binWidth, tOffset, tPlotBegin+tPlotBeginWindow)
            noSpikeTrialHistograms[cellType][structure] = hist3, bins3
    
    tableOutName = outName + '_all_trials.csv'
    with open(tableOutName, 'w') as outputTable:
        hist, bins = synapseHistograms['EXC']['Total']
        header = 'Bin start\tbin end\tbin center'
        for cellType in plotTypes:
            header += '\t'
            header += cellType
        header += '\n'
        outputTable.write(header)
        for i in range(len(bins)-1):
            line = str(bins[i])
            line += '\t'
            line += str(bins[i+1])
            line += '\t'
            line += str(0.5*(bins[i]+bins[i+1]))
            for cellType in plotTypes:
                line += '\t'
                line += str(synapseHistograms[cellType]['Total'][0][i])
            line += '\n'
            outputTable.write(line)
    
    tableOutName2 = outName + '_spike_trials.csv'
    with open(tableOutName2, 'w') as outputTable:
        hist, bins = spikeTrialHistograms['EXC']['Total']
        header = 'Bin start\tbin end\tbin center'
        for cellType in plotTypes:
            header += '\t'
            header += cellType
        header += '\n'
        outputTable.write(header)
        for i in range(len(bins)-1):
            line = str(bins[i])
            line += '\t'
            line += str(bins[i+1])
            line += '\t'
            line += str(0.5*(bins[i]+bins[i+1]))
            for cellType in plotTypes:
                line += '\t'
                line += str(spikeTrialHistograms[cellType]['Total'][0][i])
            line += '\n'
            outputTable.write(line)
    
    tableOutName3 = outName + '_no_spike_trials.csv'
    with open(tableOutName3, 'w') as outputTable:
        hist, bins = noSpikeTrialHistograms['EXC']['Total']
        header = 'Bin start\tbin end\tbin center'
        for cellType in plotTypes:
            header += '\t'
            header += cellType
        header += '\n'
        outputTable.write(header)
        for i in range(len(bins)-1):
            line = str(bins[i])
            line += '\t'
            line += str(bins[i+1])
            line += '\t'
            line += str(0.5*(bins[i]+bins[i+1]))
            for cellType in plotTypes:
                line += '\t'
                line += str(noSpikeTrialHistograms[cellType]['Total'][0][i])
            line += '\n'
            outputTable.write(line)


def scan_directory(path, fnames, suffix):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory(fname, fnames, suffix)
        elif fname.endswith(suffix):
            fnames.append(fname)
        else:
            continue

def scan_directory2(path, fnames, suffix, contains):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory2(fname, fnames, suffix, contains)
        elif fname.endswith(suffix) and fname.find(contains) != -1:
            fnames.append(fname)
        else:
            continue

if __name__ == '__main__':
    if len(sys.argv) == 5:
        folder = sys.argv[1]
        suffix = sys.argv[2]
        contains = sys.argv[3]
        outName = sys.argv[4]
        create_active_synapse_histogram(folder, suffix, contains, outName)
    elif len(sys.argv) == 6:
        folder = sys.argv[1]
        suffix = sys.argv[2]
        vmSuffix = sys.argv[3]
        contains = sys.argv[4]
        outName = sys.argv[5]
        create_active_synapse_histogram_spike_no_spike(folder, suffix, vmSuffix, contains, outName)
    else:
        print 'Error! Number of arguments is %d; should be 4 or 5' % (len(sys.argv)-1)
    
    
    
    
