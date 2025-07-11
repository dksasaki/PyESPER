# %% EXAMPLE
import glodap

from PyESPER.lir import lir
from PyESPER.nn import nn
from PyESPER.mixed import mixed
 
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
#MeasUncerts = {
#    'sal_u': [0.001], 
#    'temp_u': [0.3], 
#    'phosphate_u': [0.14], 
#    'nitrate_u':[0.5], 
#    'silicate_u': [0.03], 
#    'oxygen_u': [0.025]
#}
EstDates = data.year[L].values.tolist()
Path = "" # path works relative to the location of this script
             
EstimatesLIR, CoefficientsLIR, UncertaintiesLIR = lir(
    ['TA'], 
    Path, 
    OutputCoordinates, 
    PredictorMeasurements, 
    EstDates=EstDates,
    Equations=[1]
    )

#EstimatesNN, UncertaintiesNN = nn(
#    ['TA', 'DIC', 'pH', 'phosphate', 'nitrate', 'silicate', 'oxygen'], 
#    Path, 
#    OutputCoordinates, 
#    PredictorMeasurements, 
#    EstDates=EstDates
#    )

#EstimatesMixed, UncertaintiesMixed = mixed(
#    ['TA'], 
#    Path,
#    OutputCoordinates, 
#    PredictorMeasurements,
#    EstDates=EstDates,
#    Equations=[15]
#    )

# DEBUG
#print(EstimatesLIR['TA1']) 
#print(CoefficientsLIR['DIC2']['Coef A'][5:10])
#print(EstimatesNN['TA1'])
#print(UncertaintiesNN['TA1'])
#print(UncertaintiesLIR['TA1'])

# Optional format to pandas and save
import pandas as pd
(pd.DataFrame.from_dict(data=EstimatesLIR, orient='index')
    .to_csv('TA1LIR.csv'))
#(pd.DataFrame.from_dict(data=EstimatesNN, orient='index')
#    .to_csv('TANN.csv'))
