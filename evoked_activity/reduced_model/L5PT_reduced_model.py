'''
L5PT neuron model (reduced 2-parameter model)

@author: robert
'''

import sys
import time
import os, os.path
import glob
import neuron
import single_cell_parser as scp
import single_cell_analyzer as sca
import numpy as np
h = neuron.h

def evoked_activity(simName, cellName, evokedINHParamName, evokedGenericParamName, evokedSynDistributionFile, evokedInhSynDistributionFile, nrOfSynapses, synapseTiming):
    '''
    pre-stimulus ongoing activity
    and evoked activity
    '''
    neuronParameters = scp.build_parameters(cellName)
    evokedUpNWParameters = scp.build_parameters(evokedINHParamName)
    evokedUpGenericParameters = scp.build_parameters(evokedGenericParamName)
    scp.load_NMODL_parameters(neuronParameters)
    scp.load_NMODL_parameters(evokedUpNWParameters)
    scp.load_NMODL_parameters(evokedUpGenericParameters)
    cellParam = neuronParameters.neuron
    paramEvokedUp = evokedUpNWParameters.network
    paramEvokedUpGeneric = evokedUpGenericParameters.network
    
    cell = scp.create_cell(cellParam)
    
    evokedSynDistribution = np.loadtxt(evokedSynDistributionFile, skiprows=1, unpack=True, delimiter=',')
    cellPointsLUT = cell_points_by_distance_LUT(cell, evokedSynDistribution, 'EXC')
    evokedInhSynDistribution = np.loadtxt(evokedInhSynDistributionFile, skiprows=1, unpack=True, delimiter=',')
    cellPointsLUT_INH = cell_points_by_distance_LUT(cell, evokedInhSynDistribution, 'INH')
    offset = synapseTiming[0]
    mode = synapseTiming[1]
    median = synapseTiming[2]
    mu = np.log(median)
    sigma = np.sqrt(mu - np.log(mode))
    synTimingLognormal = (offset, mu, sigma)
            
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
    
    if 'PW' in simName:
        nrOfInhSynapses = 375 # Avg PW deflection
    elif 'SW' in simName:
        nrOfInhSynapses = 209 # Avg SW deflection
    elif 'E2' in simName:
        nrOfInhSynapses = 255 # Avg E2 deflection
    
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
        place_synapses_by_distance(evokedSynDistribution, cell, cellPointsLUT, nrOfSynapses, 'Generic')
        activate_distance_dependent_synapses(cell, tStim, synTimingLognormal, 'Generic', paramEvokedUpGeneric)
        place_synapses_by_distance(evokedInhSynDistribution, cell, cellPointsLUT_INH, nrOfInhSynapses, 'GenericINH')
        activate_distance_dependent_synapses(cell, tStim, None, 'GenericINH', paramEvokedUpGeneric)
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
        
        nRun += 1
        
        cell.re_init_cell()
        cell.synapses.pop('Generic')
        cell.synapses.pop('GenericINH')
        evokedNW.re_init_network()

        print '-------------------------------'
    
    vTraces = np.array(vTraces)
    dendTraces = []

    scp.write_all_traces(dirName+'/'+uniqueID+'_vm_all_traces.csv', t[offsetBin:], vTraces)
    for RSManager in recSiteManagers:
        for recSite in RSManager.recordingSites:
            tmpTraces = []
            for vTrace in recSite.vRecordings:
                tmpTraces.append(vTrace[offsetBin:])
            recSiteName = dirName +'/' + uniqueID + '_' + recSite.label + '_vm_dend_traces.csv'
            scp.write_all_traces(recSiteName, t[offsetBin:], tmpTraces)
            dendTraces.append(tmpTraces)
    dendTraces = np.array(dendTraces)
    
    print 'writing simulation parameter files'
    neuronParameters.save(dirName+'/'+uniqueID+'_neuron_model.param')
    evokedUpNWParameters.save(dirName+'/'+uniqueID+'_network_model.param')
    

def cell_points_by_distance_LUT(cell, distribution, preType):
    '''
    Generate LUT that contains all cell points (i.e., (section, pt) indices)
    sorted into the distance bins of the probability distribution
    Parameter distribution contains three arrays:
        distribution[0]: beginning of bins
        distribution[1]: end of bins
        distribution[2]: relative frequency of bin
    '''
    excludedStructures = []
    if preType == 'EXC':
        excludedStructures.append('Soma')
        excludedStructures.append('AIS')
        excludedStructures.append('Myelin')
    elif preType == 'INH':
        excludedStructures.append('AIS')
        excludedStructures.append('Myelin')
    else:
        errstr = 'Unknown preType: %s' % preType
        raise RuntimeError(errstr)
    LUT = [[] for i in range(len(distribution[0]))]
    for i in range(len(cell.sections)):
        sec = cell.sections[i]
	if sec.label in excludedStructures:
	    continue
        for j in range(cell.sections[i].nrOfPts):
            x = sec.relPts[j]
            somaDist = cell.distance_to_soma(sec, x)
            distIndex = -1
            for k in range(len(distribution[0])):
                if distribution[0][k] <= somaDist < distribution[1][k]:
                    distIndex = k
                    break
            
            if distIndex < 0:
                errstr = 'Could not assign point %d of section %d with distance %f to LUT!' % (j, i, somaDist)
                raise RuntimeError(errstr)
            LUT[distIndex].append((i,j))
            # make INH synapses on soma more likely:
            if preType == 'INH' and sec.label == 'Soma':
                for n in range(9):
                    LUT[distIndex].append((i,j))
    
    return LUT

def place_synapses_by_distance(distribution, cell, cellPointsByDistance, nrOfSynapses, synType):
    '''
    Place synapses onto cell according to normalized distance-dependent
    probability distribution.
    Parameter distribution contains three arrays:
        distribution[0]: beginning of bins
        distribution[1]: end of bins
        distribution[2]: relative frequency of bin
    Parameter cellPointsByDistance is a LUT that contains all 
    cell points (i.e., (section, pt) indices) sorted into the distance
    bins of the probability distribution
    '''
    #distBins = np.random.choice(len(distribution[0]), nrOfSynapses, p=distribution[2])
    samples = np.random.random_sample(nrOfSynapses)
    cumDistribution = [(np.sum(distribution[2][:i]), np.sum(distribution[2][:i+1])) for i in range(len(distribution[2]))]
    distBins = []
    for sample in samples:
        for i in range(len(cumDistribution)):
            interval = cumDistribution[i]
            if interval[0] <= sample < interval[1]:
                distBins.append(i)
                break
    
    totalDist = 0.0
    for i in range(nrOfSynapses):
        totalDist += 0.5*(distribution[0][distBins[i]] + distribution[1][distBins[i]])
        ptBin = np.random.randint(len(cellPointsByDistance[distBins[i]]))
        secID, ptID = cellPointsByDistance[distBins[i]][ptBin]
        x = cell.sections[secID].relPts[ptID]
        cell.add_synapse(secID, ptID, x, synType)
    
    print '------------------------------------'
    print 'Average synapse distance: %.1f microns' % (totalDist/nrOfSynapses)
    print '------------------------------------'

def activate_distance_dependent_synapses(cell, tStim, synapseTiming, synType, networkParam):
    '''
    Activate distance-dependent synapses according to timing distribution
    Parameter synapseTiming contains three parameters:
        synapseTiming[0]: offset
        synapseTiming[1]: mu (log-normal dist.)
        synapseTiming[2]: sigma (log-normal dist.)
    '''
    synParameters = networkParam[synType].synapses
    nrOfSynapses = len(cell.synapses[synType])
    print 'Activating %d distance-dependent synapses of type %s...' % (nrOfSynapses, synType)
    
    if synType == 'Generic':
        offset = synapseTiming[0]
        mu = synapseTiming[1]
        sigma = synapseTiming[2]
    elif synType == 'GenericINH':
        spikeProb = networkParam[synType].celltype.pointcell.probabilities
        timeBins = networkParam[synType].celltype.pointcell.intervals
        cumSpikeProb = [(np.sum(spikeProb[:i]), np.sum(spikeProb[:i+1])) for i in range(len(spikeProb))]
        samples = np.random.random_sample(nrOfSynapses)
        synTimeBins = []
        for sample in samples:
            for i in range(len(cumSpikeProb)):
                interval = cumSpikeProb[i]
                if interval[0] <= sample < interval[1]:
                    synTimeBins.append(i)
                    break
    else:
        errstr = 'Unknown synType: %s' % synType
        raise RuntimeError(errstr)
    
    #for syn in cell.synapses[synType]:
    for n in range(len(cell.synapses[synType])):
        syn = cell.synapses[synType][n]
        if synType == 'Generic':
            releaseTime = tStim + offset + np.random.lognormal(mu, sigma)
        elif synType == 'GenericINH':
            timeBin = timeBins[synTimeBins[n]]
            releaseTime = tStim + timeBin[0] + np.random.rand()*(timeBin[1] - timeBin[0])
        releaseSite = scp.PointCell([releaseTime])
        preSynCell = scp.PointCell([releaseTime])
        releaseSite.play()
        receptors = synParameters.receptors
        for recepStr in receptors.keys():
            receptor = receptors[recepStr]
            if syn.weight is None:
                syn.weight = {}
            syn.weight[recepStr] = []
            if synType == 'Generic':
                for i in range(len(receptor.weight)):
                    syn.weight[recepStr].append(receptor.weight[i])
            elif synType == 'GenericINH':
                syn.weight[recepStr].append(receptor.weight)
        syn.activate_hoc_syn(releaseSite, preSynCell, cell, receptors)
        #set properties for all receptors here
        for recepStr in receptors.keys():
            recep = receptors[recepStr]
            for param in recep.parameter.keys():
                #try treating parameters as hoc range variables,
                #then as hoc global variables
                try:
                    paramStr = 'syn.receptors[\'' + recepStr + '\'].'
                    paramStr += param + '=' + str(recep.parameter[param])
                    exec(paramStr)
                except LookupError:
                    paramStr = param + '_' + recepStr + '='
                    paramStr += str(recep.parameter[param])
                    h(paramStr)

if __name__ == '__main__':
    if len(sys.argv) == 11:
        name = sys.argv[1]
        cellName = sys.argv[2]
        networkName = sys.argv[3]
        genericNetworkName = sys.argv[4]
        genericSynFile = sys.argv[5]
        genericInhSynFile = sys.argv[6]
        nrSynapses = int(sys.argv[7])
        offset = float(sys.argv[8])
        mean = float(sys.argv[9])
        median = float(sys.argv[10])
        evoked_activity(name, cellName, networkName, genericNetworkName, genericSynFile, genericInhSynFile, nrSynapses, (offset, mean, median))
    else:
        print 'Error! Number of arguments is %d; should be 10' % (len(sys.argv)-1)
    
    
    
    
