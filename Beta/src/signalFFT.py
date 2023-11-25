'''
    Code to clean signal from noise by doing a FFT of it
'''

import numpy as np
import pandas as pd
from tqdm import tqdm

from scipy.fftpack import fft,  ifft, fftfreq

import uproot


class SignalFFT:

    def __init__(self, noise_freq:float = 0.1e10, freq_window:float = 0.1e9) -> None:
        '''
            Parameters
            ----------
            noise_freq (float) :    Frequency of the noise (GHz)
        '''
        
        self.noise_freq = noise_freq
        self.freq_window = freq_window

    def clean_signal(self, original_wfm, time) -> np.ndarray:
        '''
            Clean signal from noise by doing a FFT of it.
            Filter the signal by removing the frequencies outside a window around the signal frequency and inside that window below
            a certain amplitude.

            Parameters
            ----------
            original_wfm (np.ndarray) :     Original waveform
            time (np.ndarray) :             Time array of the waveform

            Returns
            -------
            np.ndarray :                    Filtered waveform
        '''
        
        # compute the FFT of the signal
        fft_amplitude = fft(original_wfm)
        xf = fftfreq(len(original_wfm), time[1] - time[0])
        df_fft = pd.DataFrame({'Frequency': xf, 'Amplitude': fft_amplitude})

        # filter the signal
        idmax_ampl = df_fft.query(f'Frequency > {self.noise_freq}', inplace=False)['Amplitude'].idxmax()
        freq_max = df_fft.iloc[idmax_ampl]['Frequency']
        amplitude_th =df_fft.iloc[idmax_ampl]['Amplitude']*0.9

        fft_amplitude_filtered = np.copy(fft_amplitude)
        freq_low_th = np.real(freq_max - self.freq_window)
        freq_high_th = np.real(freq_max + self.freq_window)
        fft_amplitude_filtered[(np.abs(df_fft['Frequency']) < freq_low_th) & (np.abs(df_fft['Frequency']) > freq_high_th)] = 0 
        fft_amplitude_filtered[np.abs(fft_amplitude_filtered) < amplitude_th] = 0

        return np.real(ifft(fft_amplitude_filtered))



def clean_signal(inData, wfmBranch, timeBranch, noise_freq, freq_window):
    '''
        Clean signal from noise by doing a FFT of it

        Parameters
        ----------
        inData (uproot.tree.TTreeMethods) :   Input data
        wfmBranch (str) :                     Branch of the waveform
        timeBranch (str) :                    Branch of the time array
        noise_freq (float) :                  Frequency of the noise (GHz)

        Returns
        -------
        signal_clean (array) :  Cleaned signal
    '''
    
    waveform_array = inData[f'{wfmBranch}'].array(library="np")
    time_array = inData[f'{timeBranch}'].array(library="np")
    clean_wfm_array = np.empty_like(waveform_array)

    signalFFT = SignalFFT(noise_freq=noise_freq, freq_window=freq_window)
    for idx, (wfm, time) in tqdm(enumerate(zip(waveform_array, time_array)), total=len(waveform_array), desc="Cleaning signal"):
        clean_wfm_array[idx] = signalFFT.clean_signal(wfm, time)  
    
    return clean_wfm_array



def main():
    inFile = '/Users/giogi/Documents/LabFNS2/Beta/data/input/BetaMeas_Lab2.root'
    treeName = 'wfm'
    inData = uproot.open(f'{inFile}:{treeName}')
    
    cw2 = clean_signal(inData, 'w2', 't2', 0.25e10, 0.1e9)
    cw3 = clean_signal(inData, 'w3', 't3', 0.25e10, 0.1e9)

    # add the cleaned waveform to a new tree in a new file with the same structure of the original one
    outFileName = '/Users/giogi/Documents/LabFNS2/Beta/data/output/BetaMeas_Lab2_clean.root'
    outFile = uproot.recreate(outFileName)

    t2 = inData['t2'].array(library="np")
    t3 = inData['t3'].array(library="np")
    event = inData['event'].array(library="np")

    # specify that they will be saved as vectors (required in analysis.C)
    # event_vector = event.tolist()
    # t2_vector = t2.tolist()
    # t3_vector = t3.tolist()
    # w2_vector = cw2.tolist()
    # w3_vector = cw3.tolist()
    
    outFile[treeName] = {"event": event, "t2": t2, "t3": t3, "w2": cw2, "w3": cw3}
    outFile.close()

if __name__ == "__main__":
    main()