'''
    Perform aignal amplitude to energy conversion in an existing TTree. 
    Energy to charge conversion will also be automatically perfomed.

    Author: Giorgio Alberto Lucia
    contact: giorgio.lucia@edu.unito.it
'''

import numpy as np
import uproot
from ROOT import TTree, TFile
from array import array

def EnergyConversion(rootFilePath: str, treeName: str, signalBranchName: str, a, b):
    '''
        A linear relation amplitude = a + b * energy is assumed and will be used for the conversion.
        It is also assumed that in this relation the amplitude is expressed in mV and the energy in keV.

        Parameters
        ----------
            - rootFilePath (str): path to root file where the input tree is stored
            - treeName (str): name of the TTree
            - signalBranchName (str): name of the branch where the signal values are stored
            - a (float): calibration offset
            - b (float): calibration slope
    '''
    
    with uproot.open(rootFilePath) as inFile:   tree = inFile[treeName].arrays(library='pd')
    signal_amplitude = tree[signalBranchName]

    energy = np.asarray([(x - a)/b for x in signal_amplitude], dtype=float)
    electrons = np.asarray([int(x/0.0036) for x in energy], dtype=int)

    outFile = uproot.recreate(rootFilePath)
    #del outFile[treeName]
    tree['energy'] = np.asarray(energy)
    tree['electrons'] = np.asarray(electrons)
    outFile[treeName] = tree
    
    #outFile[treeName].extend({"energy": np.asarray(energy), "electrons": np.asarray(electrons)})

if __name__ == '__main__':
        rootFilePaths = ['/home/curved/apts_new/apts-dpts-ce65-daq-software/Data/2023_06_28_SourceScan/0.0/apts_DAQ-000901010542292A_20230628_1145341M0v1.2v.root', 
                         '/home/curved/apts_new/apts-dpts-ce65-daq-software/Data/2023_06_28_SourceScan/1.2/apts_DAQ-000901010542292A_20230628_1146311M0v1.2v.root']
        treeNames = ['dataFe',
                     'dataFe' ]
        offsets = [0.1, -0.1]
        slopes = [7.67, 13.5]

        for rootFilePath, treeName, a, b in zip(rootFilePaths, treeNames, offsets, slopes):
            EnergyConversion(rootFilePath, treeName, 'seed_signal_mV', a, b)