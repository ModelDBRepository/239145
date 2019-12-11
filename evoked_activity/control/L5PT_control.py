'''
L5PT neuron model

@author: robert
'''

import sys
import time
import os, os.path
import neuron
import single_cell_parser as scp
import single_cell_analyzer as sca
import numpy as np
h = neuron.h

def evoked_activity(simName, cellName, evokedUpParamName):
    '''
    pre-stimulus ongoing activity
    and evoked activity
    '''
    neuronParameters = scp.build_parameters(cellName)
    evokedUpNWParameters = scp.build_parameters(evokedUpParamName)
    scp.load_NMODL_parameters(neuronParameters)
    scp.load_NMODL_parameters(evokedUpNWParameters)
    cellParam = neuronParameters.neuron
    paramEvokedUp = evokedUpNWParameters.network
    
    cell = scp.create_cell(cellParam)
            
    uniqueID = str(os.getpid())
    dirName = 'results/'
    dirName += simName
    if not simName.endswith('/'):
        dirName += '/'
    dirName += time.strftime('%Y%m%d-%H%M')
    dirName += '_'
    dirName += uniqueID
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    
    vTraces = []
    tTraces = []
    recordingSiteFiles = neuronParameters.sim.recordingSites
    recSiteManagers = []
    for recFile in recordingSiteFiles:
        recSiteManagers.append(sca.RecordingSiteManager(recFile, cell))
    
    nSweeps = 200
    tOffset = 100.0 # avoid numerical transients
    tStim = 245.0
    tStop = 270.0
    neuronParameters.sim.tStop = tStop
    dt = neuronParameters.sim.dt
    offsetBin = int(tOffset/dt + 0.5)
    
    nRun = 0
    while nRun < nSweeps:
        synParametersEvoked = paramEvokedUp
        
        startTime = time.time()
        evokedNW = scp.NetworkMapper(cell, synParametersEvoked, neuronParameters.sim)
        evokedNW.create_saved_network2()
        stopTime = time.time()
        setupdt = stopTime - startTime
        print 'Network setup time: %.2f s' % setupdt
        
        synTypes = cell.synapses.keys()
        synTypes.sort()
        
        print 'Testing evoked response properties run %d of %d' % (nRun+1, nSweeps)
        tVec = h.Vector()
        tVec.record(h._ref_t)
        startTime = time.time()
        scp.init_neuron_run(neuronParameters.sim, vardt=False)
        stopTime = time.time()
        simdt = stopTime - startTime
        print 'NEURON runtime: %.2f s' % simdt
        
        vmSoma = np.array(cell.soma.recVList[0])
        t = np.array(tVec)
        vTraces.append(np.array(vmSoma[offsetBin:])), tTraces.append(np.array(t[offsetBin:]))
        for RSManager in recSiteManagers:
            RSManager.update_recordings()
        
        print 'writing simulation results'
        fname = 'simulation'
        fname += '_run%04d' % nRun
        
        synName = dirName + '/' + fname + '_synapses.csv'
        print 'computing active synapse properties'
        sca.compute_synapse_distances_times(synName, cell, t, synTypes)
        preSynCellsName = dirName + '/' + fname + '_presynaptic_cells.csv'
        scp.write_presynaptic_spike_times(preSynCellsName, evokedNW.cells)
        
        nRun += 1
        
        cell.re_init_cell()
        evokedNW.re_init_network()

        print '-------------------------------'
    
    vTraces = np.array(vTraces)
    scp.write_all_traces(dirName+'/'+uniqueID+'_vm_all_traces.csv', t[offsetBin:], vTraces)
    for RSManager in recSiteManagers:
        for recSite in RSManager.recordingSites:
            tmpTraces = []
            for vTrace in recSite.vRecordings:
                tmpTraces.append(vTrace[offsetBin:])
            recSiteName = dirName +'/' + uniqueID + '_' + recSite.label + '_vm_dend_traces.csv'
            scp.write_all_traces(recSiteName, t[offsetBin:], tmpTraces)
    
    print 'writing simulation parameter files'
    neuronParameters.save(dirName+'/'+uniqueID+'_neuron_model.param')
    evokedUpNWParameters.save(dirName+'/'+uniqueID+'_network_model.param')
    

if __name__ == '__main__':
    if len(sys.argv) == 4:
        name = sys.argv[1]
        cellName = sys.argv[2]
        networkName = sys.argv[3]
        evoked_activity(name, cellName, networkName)
    else:
        print 'Error! Number of arguments is %d; should be 3' % (len(sys.argv)-1)
