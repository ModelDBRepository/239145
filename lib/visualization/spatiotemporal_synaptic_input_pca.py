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

def automated_analysis(folder, contains, vmSuffix, synapseSuffix, outName):
    vmNames = []
#    scan_directory(folder, vmNames, vmSuffix)
    scan_directory2(folder, vmNames, vmSuffix, contains)
    synNames = []
#    scan_directory(folder, synNames, synapseSuffix)
    scan_directory2(folder, synNames, synapseSuffix, contains)
    
    print 'Loading %d synapse activation files...' % len(synNames)
    synapseData = {}
    for synTrialFile in synNames:
        synapseData[synTrialFile] = scp.read_complete_synapse_activation_file(synTrialFile)
    
    synaptic_input_PCA(vmNames, synapseData, outName)
    

def synaptic_input_PCA(vmNames, synData, outName):
    '''
    calculate PCA of temporal synaptic input patterns for all trials
    with/without spikes to identify "synaptic pattern dimension" with highest variance
    Parameterization: Per trial, E/I syn. proximal/distal in 1ms bins (0-25ms)
    --> 4*25-dimensional space
    '''
    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    inhTypes = ('L1','L23Trans','L45Sym','L45Peak','L56Trans',\
                'SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6')
    # Generic excitatory synapse analysis
    plotTypes = ('EXC', 'INH')
    
    spikeTrialSynsEarlyProx = {}
    spikeTrialSynsEarlyDistal = {}
    noSpikeTrialSynsEarlyProx = {}
    noSpikeTrialSynsEarlyDistal = {}
    for cellType in plotTypes:
        spikeTrialSynsEarlyProx[cellType] = []
        spikeTrialSynsEarlyDistal[cellType] = []
        noSpikeTrialSynsEarlyProx[cellType] = []
        noSpikeTrialSynsEarlyDistal[cellType] = []
    
    tOffset = 245.0
    tStim = 253.0 # after VPM activation only
    tStimWindow = 50.0
    earlySynWindow = 25.0
    earlyWindow = 17.0 # after VPM activation only: 25-8
    binWidth = 1.0
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
    
    synNames = synData.keys()
    for n in range(len(vmNames)):
        nrSpikeTrials = 0
        nrNoSpikeTrials = 0
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
            
            synapseTimes = {}
            synapseTimesProx = {}
            synapseTimesDistal = {}
            for excType in excTypes:
                synapseTimes[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
                synapseTimesProx[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
                synapseTimesDistal[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimes['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimesProx['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimesDistal['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimes['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
            synapseTimesProx['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
            synapseTimesDistal['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
            
            for synType in activeSyns.keys():
                preCellType = synType.split('_')[0]
                for excType in excTypes:
                    if excType == preCellType:
                        for syn in activeSyns[synType]:
                            somaDist = syn[1]
                            structure = syn[4]
                            synTimes = syn[5]
                            if somaDist < 500.0:
                                synapseTimesProx['EXC']['Total'].extend(synTimes)
                            else:
                                synapseTimesDistal['EXC']['Total'].extend(synTimes)
                for inhType in inhTypes:
                    if inhType == preCellType:
                        for syn in activeSyns[synType]:
                            somaDist = syn[1]
                            structure = syn[4]
                            synTimes = syn[5]
                            if somaDist < 500.0:
                                synapseTimesProx['INH']['Total'].extend(synTimes)
                            else:
                                synapseTimesDistal['INH']['Total'].extend(synTimes)
            
            if not trialWithSpikes[n][trialNr]:
                nrNoSpikeTrials += 1
                for cellType in plotTypes:
                    synTimesProx = synapseTimesProx[cellType]['Total']
                    synTimesDistal = synapseTimesDistal[cellType]['Total']
                    tSynProxTmpList = []
                    tSynDistalTmpList = []
                    for tSyn in synTimesProx:
                        if tOffset <= tSyn < tOffset + earlySynWindow:
                            tSynProxTmpList.append(tSyn-tOffset)
                    for tSyn in synTimesDistal:
                        if tOffset <= tSyn < tOffset + earlySynWindow:
                            tSynDistalTmpList.append(tSyn-tOffset)
                    proxHist, tmpBins = np.histogram(tSynProxTmpList, bins=range(26))
                    distalHist, tmpBins = np.histogram(tSynDistalTmpList, bins=range(26))
                    noSpikeTrialSynsEarlyProx[cellType].append(proxHist)
                    noSpikeTrialSynsEarlyDistal[cellType].append(distalHist)
            
            elif trialWithSpikes[n][trialNr]:
                nrSpikeTrials += 1
                tSpikeReference = np.min(trialSpikeTimes[n][trialNr])
                for cellType in plotTypes:
                    synTimesProx = synapseTimesProx[cellType]['Total']
                    synTimesDistal = synapseTimesDistal[cellType]['Total']
                    tSynProxTmpList = []
                    tSynDistalTmpList = []
                    for tSyn in synTimesProx:
                        if tOffset <= tSyn < tOffset + earlySynWindow:
                            tSynProxTmpList.append(tSyn-tOffset)
                    for tSyn in synTimesDistal:
                        if tOffset <= tSyn < tOffset + earlySynWindow:
                            tSynDistalTmpList.append(tSyn-tOffset)
                    proxHist, tmpBins = np.histogram(tSynProxTmpList, bins=range(26))
                    distalHist, tmpBins = np.histogram(tSynDistalTmpList, bins=range(26))
                    spikeTrialSynsEarlyProx[cellType].append(proxHist)
                    spikeTrialSynsEarlyDistal[cellType].append(distalHist)
            
        print ''
        print 'Nr of trials with spike: %d' % nrSpikeTrials
        print 'Nr of trials without spike: %d' % nrNoSpikeTrials
        print 'mean spike time = %.1fms' % np.mean([tSpike for trace in trialSpikeTimes[n] for tSpike in trace])
    
    trialSpikeList = []
    totalHistList = []
    totalSpikeTrials = len(spikeTrialSynsEarlyProx['EXC'])
    totalNoSpikeTrials = len(noSpikeTrialSynsEarlyProx['EXC'])
    # concatenate E/I/prox/distal histograms in the order:
    # prox. E/prox. I/distal E/distal I
    for i in range(totalSpikeTrials):
        trialSpikeList.append(1)
        trialHist = []
        trialHist.extend(spikeTrialSynsEarlyProx['EXC'][i])
        trialHist.extend(spikeTrialSynsEarlyProx['INH'][i])
        trialHist.extend(spikeTrialSynsEarlyDistal['EXC'][i])
        trialHist.extend(spikeTrialSynsEarlyDistal['INH'][i])
        totalHistList.append(trialHist)
    for i in range(totalNoSpikeTrials):
        trialSpikeList.append(0)
        trialHist = []
        trialHist.extend(noSpikeTrialSynsEarlyProx['EXC'][i])
        trialHist.extend(noSpikeTrialSynsEarlyProx['INH'][i])
        trialHist.extend(noSpikeTrialSynsEarlyDistal['EXC'][i])
        trialHist.extend(noSpikeTrialSynsEarlyDistal['INH'][i])
        totalHistList.append(trialHist)
    
    allData = np.array(totalHistList)
    dataMean = np.mean(allData)
    allData = allData - dataMean
    eigenvectors, eigenvals, V = np.linalg.svd(allData.T, full_matrices=False)
    projectedData = np.dot(allData, eigenvectors).transpose()
    PC1LoadVec = eigenvectors.transpose()[0]
    proxELoad = PC1LoadVec[:25]
    proxILoad = PC1LoadVec[25:50]
    distalELoad = PC1LoadVec[50:75]
    distalILoad = PC1LoadVec[75:]
    PC2LoadVec = eigenvectors.transpose()[1]
    proxELoad2 = PC2LoadVec[:25]
    proxILoad2 = PC2LoadVec[25:50]
    distalELoad2 = PC2LoadVec[50:75]
    distalILoad2 = PC2LoadVec[75:]
    
    with open(outName+'_PC1_PC2_Load.csv', 'w') as outFile1:
        header = 'time (ms)\tPC1 E prox load\tPC1 I prox load\tPC1 E distal load\tPC1 I distal load\t'
        header += 'PC2 E prox load\tPC2 I prox load\tPC2 E distal load\tPC2 I distal load\n'
        outFile1.write(header)
        for i in range(25):
            line = str(i+0.5)
            line += '\t'
            line += str(proxELoad[i])
            line += '\t'
            line += str(proxILoad[i])
            line += '\t'
            line += str(distalELoad[i])
            line += '\t'
            line += str(distalILoad[i])
            line += '\t'
            line += str(proxELoad2[i])
            line += '\t'
            line += str(proxILoad2[i])
            line += '\t'
            line += str(distalELoad2[i])
            line += '\t'
            line += str(distalILoad2[i])
            line += '\n'
            outFile1.write(line)
    with open(outName+'_PC1_PC2.csv', 'w') as outFile2:
        header = 'spike trial 1/0\tPC1\tPC2\n'
        outFile2.write(header)
        for i in range(len(projectedData[0])):
            line = str(trialSpikeList[i])
            line += '\t'
            line += str(projectedData[0][i])
            line += '\t'
            line += str(projectedData[1][i])
            line += '\n'
            outFile2.write(line)


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
    if len(sys.argv) == 6:
        folder = sys.argv[1]
        contains = sys.argv[2]
        vmSuffix = sys.argv[3]
        synapseSuffix = sys.argv[4]
        outName = sys.argv[5]
        automated_analysis(folder, contains, vmSuffix, synapseSuffix, outName)
    else:
        print 'Error! Number of arguments is %d; should be 5' % (len(sys.argv)-1)
    
    
    
    
