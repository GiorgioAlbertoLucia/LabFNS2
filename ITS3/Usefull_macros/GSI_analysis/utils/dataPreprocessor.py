'''
    Script to load, read and display data from .npy and .npz files.
    Independent classes to be included in stand-alone scripts

    Author: Giorgio Alberto Lucia
    contact: giorgio.lucia@edu.unito.it
'''

import os
import numpy as np
from scipy.interpolate import interp1d
from array import array
from alive_progress import alive_bar

from ROOT import TTree, TFile
#from utilsGG.patternIdentifier import PatternIdentifier

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

class NpyFileManager:
    '''
    Class to upload a .npy file, display its contents and convert it in a TTree storing it inside an output .root file.
    The np.array contained in the .npy file should have the shape (n_events, n_cols, n_rows, n_timeframes), where n_events is the number of events collected,
    n_cols and n_rows are the number of columns and rows in the pixel matrix and n_timeframes the number of timeframes observed for each event.

    To begin, load the data with loadData()

    Attributes
    ----------
        - _inputFilePath (str): path to the input file
        - _threshold (int): threshold value (in ADC chn) used during data acquisition
        - _thresholds ([int]): threshold value for each pixel (different if a calibration is done)
        - _data (np.array): data stored in the input file
        - _data_cal (np.array): data stored in the input file calibrated in mV
        - _outputFile (TFile): file to save the TTree inside
        - _outputTree (TTree): TTree to store data in
        - _nPixel (int): number of pixel of the sensor
        - _branchDict (dict): dictionary with branches to load on the TTree
        - _calibrated (bool): whether the data has been calibrated or not
        - _maxPixelComb (int): maximum number of pixel whose combinations will be considered as pattern
    '''

    def __init__(self, inputFilePath: str, threshold: int, maxPixelComb: int = 3):
        
        self._inputFilePath = inputFilePath
        self._data = None
        self._data_cal = None
        self._threshold = threshold
        self._thresholds = None
        
        self._outputFile = None
        self._outputTree = None

        self._nPixel = 0
        self._branchDict = None
        self._calibrated = False

        self._maxPixelComb = maxPixelComb
        

    def loadData(self):
        '''
        Load data from file and initialize number of pixels and the branch dictionary.
        The threshold list is also initialized.
        '''
        
        try:                           
            with alive_bar(title=(color.PURPLE +'Loading file...'+ color.DARKCYAN + self._inputFilePath + color.END)) as bar: 
                self._data = np.load(self._inputFilePath)
        except FileNotFoundError:       raise FileNotFoundError(f"File '{+ color.DARKCYAN + self._inputFilePath + color.END}' not found.")
        except Exception as e:          raise Exception(f"An error occurred while loading the file: {str(e)}")
        
        self._nPixel = self._data.shape[1] * self._data.shape[2]
        
        self._branchDict = {'seed_signal_chn':      array('d', [0.]),
                            'matrix_signal_chn':    array('d', [0.]),
                            'signal_over_th_chn':   array('d', [0.]),
                            'seed_signal_mV':       array('d', [0.]),
                            'matrix_signal_mV':     array('d', [0.]),
                            'signal_over_th_mV':    array('d', [0.]),
                            'seed':                 array('i', [0]),
                            'cluster_size':         array('i', [0]),
                            'cluster_shape':        array('i', [0])         
                            }
        self._thresholds = np.empty(self._nPixel)
        self._thresholds.fill(self._threshold)

        
    def display(self):
        '''
        Display the input data inside the .npy file. 
        Information displayed are
        * the entire content of the file
        * the shape of the file
        * the type of the data
        '''

        print("Loaded data:", flush=True)
        print(self._data, flush=True)
        print("Data shape:", self._data.shape, flush=True)
        print("Data type:", self._data.dtype, flush=True)

    def calibrate(self, npzFilePath: str, save: bool = False):
        '''
        Calibrate ADC data to voltage values using a calibration .npz file.

        Parameters
        ----------
            - npzFilePath (str): path to the .npz file needed for calbration
            - save (bool): whether a copy of the calibrated data should be saved as an .npy
        '''
        
        assert self._data is not None, 'Error: need to load data first'
        
        calibration = np.load(npzFilePath)
        ys = calibration['vres_range']
        xs = calibration['baseline_all']
        data_calibrated = np.zeros((self._data.shape[0], self._data.shape[1], self._data.shape[2], self._data.shape[3]))

        #with alive_bar(self._nPixel, title='Calibration...') as bar:
        print('\n' + color.PURPLE + 'Calibration...' + color.END, flush=True)
        pixel_idx=0
        for row in range(self._data.shape[1]):
            for column in range(self._data.shape[2]): 
                interpolation = interp1d(xs[pixel_idx, :], ys, kind='cubic', fill_value='extrapolate')
                data_calibrated[:, row, column, :] = interpolation(self._data[:, row, column, :])
                #self._thresholds[pixel_idx] = interpolation(self._threshold)
                pixel_idx += 1
                #bar()
        self._data_cal = data_calibrated
        
        if save:    np.save(os.path.splitext(self._inputFilePath)[0]+'_calibrated.npy', self._data)
        self._calibrated = True


    def plantTree(self, treeName: str):
        '''
        Function to initialize the TTree. 

        Parameters
        ----------
            - treeName (str): name of the TTree created
        '''

        assert self._data is not None, 'Error: need to load data first'
        
        self._outputTree = TTree(treeName, 'Output Tree')

        self._outputTree.Branch('seed_signal_chn', self._branchDict['seed_signal_chn'], 'seed_signal_chn/D')  
        self._outputTree.Branch('matrix_signal_chn', self._branchDict['matrix_signal_chn'], 'matrix_signal_chn/D')  
        self._outputTree.Branch('signal_over_th_chn', self._branchDict['signal_over_th_chn'], 'signal_over_th_chn/D') 
   
        self._outputTree.Branch('seed_signal_mV', self._branchDict['seed_signal_mV'], 'seed_signal_mV/D')  
        self._outputTree.Branch('matrix_signal_mV', self._branchDict['matrix_signal_mV'], 'matrix_signal_mV/D')  
        self._outputTree.Branch('signal_over_th_mV', self._branchDict['signal_over_th_mV'], 'signal_over_th_mV/D') 

        self._outputTree.Branch('seed', self._branchDict['seed'], 'seed/I')
        self._outputTree.Branch('cluster_size', self._branchDict['cluster_size'], 'cluster_size/I')
        self._outputTree.Branch('cluster_shape', self._branchDict['cluster_shape'], 'cluster_shape/I')

    def fillTree(self, **kwargs):
        '''
        Fill the TTree with data stored in data.

        Parameters
        ----------
            Possible kwargs:
            - option (str): defines the type of data and how the seed signal should be defined. 
                            Accepted values are 'data' (default), 'pedestal', 'zero_point'
            - v_reset (float): vreset value in mV. Mandatory for 'pedestal' option.
        '''

        assert self._data is not None, 'Error: need to load data first'
        assert self._outputTree is not None, 'Error: need to plant tree first'

        with alive_bar(int(self._outputTree.GetEntries()/10000), title=(color.PURPLE + 'Filling the TTree...' + color.END)) as bar:

            #pattern_identifier = PatternIdentifier(self._data.shape[0], self._data.shape[1], self._maxPixelComb)
            for ev_idx, event in enumerate(self._data):

                if ev_idx%10000==0: print(f'Processing event {ev_idx}...', flush=True)

                signals_chn = np.zeros(self._nPixel)
                signals_mV = np.zeros(self._nPixel)

                self._branchDict['matrix_signal_chn'][0] = 0.
                self._branchDict['signal_over_th_chn'][0] = 0.
                self._branchDict['matrix_signal_mV'][0] = 0.
                self._branchDict['signal_over_th_mV'][0] = 0.

                hit_pixels = []

                pixel_idx = 0
                for row_idx, row in enumerate(event):
                    for col_idx, column in enumerate(row): 

                        baseline = np.mean(column[0:100])
                        if 'option' in kwargs:
                            if kwargs['option'] == 'zero_point':    signal = column[100] - baseline
                            if kwargs['option'] == 'pedestal':
                                assert 'v_reset' in kwargs, 'For the "pedestal" option, a "v_reset" value should be provided'
                                print('\n' + color.BOLD + color.YELLOW + 'Warning: ' + color.END + 'with the pedestal option only the valtage results will be coherent, the ADC results will match the data option.')
                                signal = abs(np.min(column) - baseline)
                        if 'option' not in kwargs or kwargs['option'] == 'data':    signal = abs(np.min(column) - baseline)

                        signals_chn[pixel_idx] = signal
                        if signal > self._thresholds[pixel_idx]:    hit_pixels.append(pixel_idx)

                        if self._calibrated:
                            column_mV = self._data_cal[ev_idx][row_idx][col_idx]
                            baseline_mV = np.mean(column_mV[0:100])
                            if 'option' in kwargs:
                                if kwargs['option'] == 'zero_point':    signal_mV = column_mV[100] - baseline_mV
                                if kwargs['option'] == 'pedestal':
                                    assert 'v_reset' in kwargs, 'For the "pedestal" option, a "v_reset" value should be provided'
                                    assert kwargs['v_reset'] is not None, '"v_reset" value is None.'
                                    signal_mV = abs(baseline_mV - kwargs['v_reset'])
                            if 'option' not in kwargs or kwargs['option'] == 'data':    signal_mV = abs(np.min(column_mV) - baseline_mV)
                            signals_mV[pixel_idx] = signal_mV

                        pixel_idx += 1

                if 'option' in kwargs:
                    if kwargs['option'] == 'pedestal':      self._branchDict['seed_signal_chn'][0] = np.max(signals_chn)
                    elif kwargs['option'] == 'zero_point':  self._branchDict['seed_signal_chn'][0] = np.mean(signals_chn)
                if 'option' not in kwargs or kwargs['option'] == 'data':    
                    self._branchDict['seed_signal_chn'][0] = np.max(signals_chn)
                
                for signal in signals_chn:  self._branchDict['matrix_signal_chn'][0] += signal
                for hit_idx in hit_pixels:  self._branchDict['signal_over_th_chn'][0] += signals_chn[hit_idx]

                if self._calibrated:  
                    if 'option' in kwargs:
                        if kwargs['option'] == 'pedestal' or kwargs['option'] == 'zero_point':
                            self._branchDict['seed_signal_mV'][0] = np.mean(signals_mV)
                    if 'option' not in kwargs or kwargs['option'] == 'data':    
                        self._branchDict['seed_signal_mV'][0] = np.max(signals_mV)  
                    
                    for signal in signals_mV:  self._branchDict['matrix_signal_mV'][0] += signal
                    for hit_idx in hit_pixels:  self._branchDict['signal_over_th_mV'][0] += signals_mV[hit_idx]

                self._branchDict['seed'][0] = np.argmax(signals_chn)
                self._branchDict['cluster_size'][0] = len(hit_pixels)

                #cluster_matrix = [[0 for _ in range(self._data.shape[1])] for _ in range(self._data.shape[0])]
                #for hit in hit_idx: cluster_matrix[hit//self._data.shape[1]][hit%self._data.shape[1]] = 1
                #self._branchDict['cluster_shape'][0] = pattern_identifier.identifyPattern(cluster_matrix)

                self._outputTree.Fill()
                if ev_idx%10000==0: bar()
            
    def saveTree(self, outFilePath: str):
        '''
        Save the TTree in a TFile.
        Delete the TTree and data afterwards.

        Parameters
        ----------
            - outFilePath (str): path to the output file
        '''

        #if self._calibrated:    outFilePath = os.path.splitext(outFilePath)[0] + '_calibrated.root'
        self._outputFile = TFile(outFilePath, 'recreate')
        self._outputTree.Write()
        self._outputFile.Close()
        print(f'TTree stored in {color.DARKCYAN + outFilePath + color.END}', flush=True)
        del self._outputTree
        


# Example usage
if __name__ == "__main__":
    file_path = "/home/curved/apts_new/apts-dpts-ce65-daq-software/Data/dataAcqFe/SecondAcq_10082023/0.0/source/apts_DAQ-000901010542292A_20230810_182606.npy"
    
    reader = NpyFileManager(file_path, 88)
    reader.loadData()
    #reader.display()
    reader.calibrate("/home/curved/apts_new/apts-dpts-ce65-daq-software/Data/dataAcqFe/SecondAcq_10082023/0.0/gain/apts_gain_DAQ-000901010542292A_20230810_182606_analysed.npz", True)
    reader.plantTree('tree')
    reader.fillTree()
    reader.saveTree("/home/curved/apts_new/apts-dpts-ce65-daq-software/Data/dataAcqFe/SecondAcq_10082023/0.0/source/apts_DAQ-000901010542292A_20230810_182606.root")
