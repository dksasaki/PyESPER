# %% EXAMPLE
import glodap

from PyESPER.lir import lir
from PyESPER.nn import nn
 
data = glodap.atlantic() 

L = slice(0, -1)
PredictorMeasurements = {
    k: data[k][L].values.tolist()
    for k in [
        "salinity",
        "temperature",
        "phosphate",
        "nitrate",
        "silicate",
        "oxygen",
    ]
}
OutputCoordinates = {
    k: data[k][L].values.tolist()
    for k in [
        "longitude",
        "latitude",
        "depth",
    ]
}
MeasUncerts = {
    'sal_u': [0.001], 
    'temp_u': [0.3], 
    'phosphate_u': [0.14], 
    'nitrate_u':[0.5], 
    'silicate_u': [0.03], 
    'oxygen_u': [0.025]
}
EstDates = data.year[L].values.tolist()
Path = "" # path works relative to the location of this script
             
EstimatesLIR, CoefficientsLIR, UncertaintiesLIR = lir(
    ['TA'], 
    Path, 
    OutputCoordinates, 
    PredictorMeasurements, 
    EstDates=EstDates
    )

EstimatesNN, UncertaintiesNN = nn(
    ['TA'], 
    Path, 
    OutputCoordinates, 
    PredictorMeasurements, 
    EstDates=EstDates 
    )

# DEBUG
print(EstimatesLIR['TA1'][25:30]) 
print(CoefficientsLIR['DIC2']['Coef A'][35:40])
print(EstimatesNN['pH11'][100:105])
print(UncertaintiesNN['phosphate16'][0:5])
print(UncertaintiesLIR['TA2'][15:20])

# Optional format to pandas and save
import pandas as pd
(pd.DataFrame.from_dict(data=EstimatesLIR, orient='index')
    .to_csv('TALIR.csv'))
(pd.DataFrame.from_dict(data=EstimatesNN, orient='index')
    .to_csv('TANN.csv'))
