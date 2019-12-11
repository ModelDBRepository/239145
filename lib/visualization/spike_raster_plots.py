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
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def create_spike_times_files(folder, outName):
    '''
    automatically creates spike time files
    from all somatic recording files found below folder
    '''
    suffix = 'vm_all_traces.csv'
    fnames = []
    scan_directory(folder, fnames, suffix)
    
    for fname in fnames:
        outName1 = fname[:-4] + '_spike_times.csv'
        #if os.path.exists(outName1):
        #    print 'Spike times file %s already exists; skipping...' % outName1
        #    continue
        print 'Loading traces from file %s' % fname
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        t = data[0]
        allSpikeTimes = {}
        for i in range(1, len(data)):
            v = data[i]
            spikeTimes = sca.simple_spike_detection(t, v)
            allSpikeTimes[i-1] = spikeTimes
        
        scp.write_spike_times_file(outName1, allSpikeTimes)
        splitFName = fname.split('/')
        spikeTimesName = outName
        if not spikeTimesName.endswith('/'):
            spikeTimesName += '/'
        whisker = splitFName[-3].split('_')[0]
        location = splitFName[-3].split('_')[-1]
        spikeTimesName += whisker
        spikeTimesName += '_deflection_location_'
        spikeTimesName += location
        spikeTimesName += '_spike_times.csv'
        scp.write_spike_times_file(spikeTimesName, allSpikeTimes)

def create_spike_raster_plots(folder, suffix, outName):
    '''
    load all traces, compute spike times
    and create raster plots
    '''
    fnames = []
    scan_directory(folder, fnames, suffix)
    
    print 'Creating spike raster plots from %d files' % len(fnames)
    
    tOffset = 100.0
    tStop = 270.0
    for n in range(len(fnames)):
        fname = fnames[n]
        print 'Loading traces from file %s' % fname
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        t = data[0]
        allSpikeTimes = []
        ax = []
        fig = plt.figure(2*n+1)
        for i in range(1, len(data)):
            fig.add_subplot(2,1,1)
            ax.append(plt.plot(t, data[i]))
            fig.add_subplot(2,1,2)
            v = data[i]
            spikeTimes = sca.simple_spike_detection(t, v)
            allSpikeTimes.append(spikeTimes)
            spikes = [i for time in spikeTimes]
            ax.append(plt.plot(spikeTimes, spikes, 'k|'))
        fig.add_subplot(2,1,1)
        plt.xlabel('t [ms]')
        plt.ylabel('Vm [mV]')
        plt.xlim([tOffset,tStop])
        fig.add_subplot(2,1,2)
        plt.xlabel('t [ms]')
        plt.ylabel('trial nr.')
        plt.xlim([tOffset,tStop])
        
        splitFName = fname.split('/')
        whisker = splitFName[-3].split('_')[0]
        location = splitFName[-3].split('_')[-1]
        rasterName = outName
        if not rasterName.endswith('/'):
            rasterName += '/'
        rasterName += whisker
        rasterName += '_deflection_location_'
        rasterName += location
        rasterName += '_spike_raster_plot.pdf'
        plt.savefig(rasterName)
        
        hist, bins = sca.PSTH_from_spike_times(allSpikeTimes, 1.0, tOffset, tStop)
        fig = plt.figure(2*n+2)
        plt.bar(bins[:-1], hist, color='b', width=1)
        plt.xlabel('t [ms]')
        plt.ylabel('AP/stim')
        plt.xlim([tOffset,tStop])
        psthName = outName
        if not psthName.endswith('/'):
            psthName += '/'
        psthName += whisker
        psthName += '_deflection_location_'
        psthName += location
        psthName += '_PSTH_total_1ms.pdf'
        plt.savefig(psthName)
        scp.write_PSTH(psthName[:-4]+'.csv', hist, bins)

def create_PSTH(folder, suffix, contains, outName):
    '''
    load all traces, compute spike times
    and create raster plots
    '''
    fnames = []
    scan_directory2(folder, fnames, suffix, contains)
    
    print 'Creating PSTH from %d files' % len(fnames)
    
    if not os.path.exists(os.path.dirname(outName)):
        os.makedirs(os.path.dirname(outName))
    
    tOffset = 100.0
    tStop = 270.0
    allSpikeTimes = []
    for n in range(len(fnames)):
        fname = fnames[n]
        print 'Loading traces from file %s' % fname
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        t = data[0]
        for i in range(1, len(data)):
            v = data[i]
            spikeTimes = sca.simple_spike_detection(t, v)
            allSpikeTimes.append(spikeTimes)
    
    binWidth = 1.0
    hist, bins = sca.PSTH_from_spike_times(allSpikeTimes, binWidth, tOffset, tStop)
    offset = 0.5*(bins[1] - bins[0])
    fig = plt.figure(2*n+2)
    plt.bar(bins[:-1], hist, color='b', width=binWidth)
    plt.xlabel('t [ms]')
    plt.ylabel('AP/stim')
    plt.xlim([tOffset,tStop])
    outName += '_PSTH_total_%.1fms' % binWidth
    plt.savefig(outName+'.pdf')
    scp.write_PSTH(outName+'.csv', hist, bins)

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
    if len(sys.argv) == 3:
        folder = sys.argv[1]
        outName = sys.argv[2]
        create_spike_times_files(folder, outName)
    elif len(sys.argv) == 4:
        folder = sys.argv[1]
        suffix = sys.argv[2]
        outName = sys.argv[3]
        create_spike_raster_plots(folder, suffix, outName)
    elif len(sys.argv) == 5:
        folder = sys.argv[1]
        suffix = sys.argv[2]
        contains = sys.argv[3]
        outName = sys.argv[4]
        create_PSTH(folder, suffix, contains, outName)
    else:
        print 'Error! Number of arguments is %d; should be 2, 3 or 4' % (len(sys.argv)-1)
    
    
    
   
