'''
    Code to clean signal from noise by doing a FFT of it
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq, fftshift

import uproot

def clean_signal(inFile, treeName, branchName, dt):
    '''
        Clean signal from noise by doing a FFT of it

        Parameters
        ----------
        inFile (str) :          Input file name
        treeName (str) :        Name of the tree
        branchName (str) :      Name of the branch
        dt (float) :            Time step of the signal

        Returns
        -------
        signal_clean (array) :  Cleaned signal
    '''

    # Number of points
    

    inData = uproot.open(f'{inFile}:{treeName}')
    waveform_array = inData[f'{branchName}'].array(library="np")
    
    nsample = waveform_array[0].size                  # Number of points sampled for each event 
    print(waveform_array[0].size)

    waveform = np.concatenate(waveform_array)         # Concatenate all the events in a single array
    #print(np.mean(waveform))  
    #print(np.std(waveform[:200]))  

    waveform_fft = np.fft.fft(waveform)               # FFT of the signal
    plt.plot(waveform_fft)
    plt.show()

    threshold = 800
    mask = np.abs(waveform_fft) < threshold
    filtered_fft = waveform_fft * mask

    filtered_signal = np.fft.ifft(filtered_fft)

    # Plot the waveform and filtered signal
    #fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    #axs[0].plot(waveform[:10000])
    #axs[0].set_title('Waveform')
    #axs[1].plot(filtered_signal[:10000])
    #axs[1].set_title('Filtered Signal')
    #plt.show()
    
    # Plot the first 3000 points of waveform
    #plt.plot(waveform[:10000])
    #plt.show()

    



def main():
    inFile = '/Users/giogi/Documents/LabFNS2/Beta/data/input/BetaMeas_Lab2.root'
    treeName = 'wfm'
    branchName = 'w3'
    dt = 0.1
    clean_signal(inFile, treeName, branchName, dt)

if __name__ == "__main__":
    main()