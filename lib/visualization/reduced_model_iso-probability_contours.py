import sys
import os, os.path
import glob
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import single_cell_analyzer as sca
import single_cell_parser as scp


def iso_probability_plot(folder, suffix, whisker, outName):
    fnames = []
    scan_directory(folder, fnames, suffix)
    
    synNrTimingFilenames = {}
    for fname in fnames:
        basePath = fname.split('/')[-3]
        nSyn = int(basePath.split('_')[1])
        synTiming = float(basePath.split('_')[5])
        if not synNrTimingFilenames.has_key(nSyn):
            synNrTimingFilenames[nSyn] = {}
        synNrTimingFilenames[nSyn][synTiming] = fname
    
    tBegin = 8.0
    if whisker == 'PW':
        window = 7.0
    elif whisker == 'SW':
        window = 9.0
    elif whisker == 'E2':
        window = 17.0
    
    synNrTimingProbs = {}
    offset = 245.0
    for nSyn in synNrTimingFilenames.keys():
        synNrTimingProbs[nSyn] = {}
        for synTiming in synNrTimingFilenames[nSyn].keys():
            fname = synNrTimingFilenames[nSyn][synTiming]
            print 'Analyzing spike times in file %s' % fname
            trialSpikeTimes = scp.read_spike_times_file(fname)
            nrOfTrials = len(trialSpikeTimes.keys())
            spikeTrials = 0.0
            for trial in trialSpikeTimes.keys():
                for tSpike in trialSpikeTimes[trial]:
                    if tBegin <= tSpike - offset < tBegin + window:
                        spikeTrials += 1
                        break
            spikeProb = spikeTrials/nrOfTrials
            synNrTimingProbs[nSyn][synTiming] = spikeProb
            print 'Window: %.1f - %.1f ms: spike prob = %.2f' % (tBegin, tBegin+window, spikeProb)
    
    numbers = np.array(range(25,350,25), dtype=np.float64)
    timings = np.array(range(2,25,1), dtype=np.float64)
    #numberTimingMesh = np.meshgrid(numbers,timings)
    #spikeProbMesh = np.zeros_like(numberTimingMesh[0])
    #for i in range(len(spikeProbMesh)-1):
        #for j in range(len(spikeProbMesh[i])-1):
            #number = int(numberTimingMesh[0][i][j])
            #timing = numberTimingMesh[1][i][j]
            #try:
                #spikeProbMesh[i][j] = synNrTimingProbs[number][timing]
            #except KeyError:
                #spikeProbMesh[i][j] = 0.0
    numberTimingMesh = np.meshgrid(timings, numbers)
    spikeProbMesh = np.zeros_like(numberTimingMesh[0])
    for i in range(len(spikeProbMesh)-1):
        for j in range(len(spikeProbMesh[i])-1):
            timing = numberTimingMesh[0][i][j]
            number = int(numberTimingMesh[1][i][j])
            try:
                spikeProbMesh[i][j] = synNrTimingProbs[number][timing]
            except KeyError:
                spikeProbMesh[i][j] = 0.0
    
    plt.figure(1)
    if whisker == 'PW':
        isocontours = [0.02, 0.04, 0.08, 0.15, 0.18, 0.3, 0.6, 0.9]
    elif whisker == 'SW' or whisker == 'E2':
        isocontours = [0.01, 0.02, 0.04, 0.08, 0.15, 0.3, 0.6, 0.9]
    CS = plt.contour(numberTimingMesh[0], numberTimingMesh[1], spikeProbMesh, isocontours)
    plt.clabel(CS, inline=1, fontsize=10, inline_spacing=-15)
    titleStr = 'Reduced SPD model - iso-probability contours for %s touch' % whisker
    plt.title(titleStr)
    plt.xlabel('Synapse timing (median; ms)')
    plt.ylabel('Synapse number')
    plt.xlim(2.0,24.0)
    plt.ylim(0,325)
    plt.xticks(np.arange(2,24,1.0), [str(i) for i in range(2,24,1)])
    if whisker == 'PW':
        plt.yticks(np.arange(25,325,25.0), [str(0.76*i) for i in range(25,325,25)]) # PW: 76% syn <500 microns
    elif whisker == 'SW':
        plt.yticks(np.arange(25,325,25.0), [str(0.69*i) for i in range(25,325,25)]) # SW: 69% syn <500 microns
    elif whisker == 'E2':
        plt.yticks(np.arange(25,325,25.0), [str(0.68*i) for i in range(25,325,25)]) # E2: 68% syn <500 microns
    cbar = plt.colorbar()
    cbar.set_label('Spike prob')
    plt.savefig(outName+'_isocontours.pdf')
    

def scan_directory(path, fnames, suffix):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory(fname, fnames, suffix)
        elif fname.endswith(suffix):
            fnames.append(fname)
        else:
            continue

if __name__ == '__main__':
    if len(sys.argv) == 5:
        folder = sys.argv[1]
        suffix = sys.argv[2]
        whisker = sys.argv[3]
        outName = sys.argv[4]
        iso_probability_plot(folder, suffix, whisker, outName)
    else:
        print 'Error! Number of arguments is %d; should be 4' % (len(sys.argv)-1)
