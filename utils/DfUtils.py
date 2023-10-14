import pandas as pd
import numpy as np

from ROOT import TGraphErrors

def GetGraphErrorsFromPandas(df:pd.DataFrame):
    '''
    Function that returns a TGraphError based on a pandas dataframe, containing the following columns:
    x_values, x_errors, y_values, y_errors
    '''
    return TGraphErrors(len(df), np.asarray(df.iloc[:,0],'d'), np.asarray(df.iloc[:,2],'d'), np.asarray(df.iloc[:,1],'d'), np.asarray(df.iloc[:,3],'d'))

def GetGraphErrorsFromCSV(infile:str, **kwargs):
    '''
    Function that returns a TGraphError based on a input csv file, containing the following columns:
    x_values, x_errors, y_values, y_errors
    '''
    df = pd.read_csv(infile, **kwargs)
    return GetGraphErrorsFromPandas(df)