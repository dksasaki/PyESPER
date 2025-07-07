def adjust_pH_DIC(DesiredVariables, VerboseTF, Dates, Est_pre={}, PredictorMeasurements={}, OutputCoordinates={}, **kwargs):

    """
    Adjusting pH and DIC for anthropogenic carbon.
    """

    import numpy as np
    import math    
    from PyESPER.simplecantestimatelr import simplecantestimatelr
    
    YouHaveBeenWarnedCanth = False

    Cant_adjusted={}
    combos2 = list(Est_pre.keys())
    values2 = list(Est_pre.values())
    values2 = np.array(values2)#.T
    Cant = []
    Cant2002 = []
    if "EstDates" in kwargs and ("DIC" in DesiredVariables or "pH" in DesiredVariables):      
        if not YouHaveBeenWarnedCanth:
            if VerboseTF:
                print("Estimating anthropogenic carbon for PyESPER_NN.")
            longitude = np.mod(OutputCoordinates["longitude"], 360)
            latitude = np.array(OutputCoordinates["latitude"])
            depth = np.array(OutputCoordinates["depth"])
            Canta, Cant2002a = simplecantestimatelr(Dates, longitude, latitude, depth)
            Canta = np.array(Canta)
            Cant2002a = np.array(Cant2002a)  
            YouHaveBeenWarnedCanth = True
            for element in Canta:
                Cant.append(element)
            for element2002 in Cant2002a:
                Cant2002.append(element2002)
                
        for combo, a in zip(combos2, values2):
            a2 = np.array(a)
            a2 = np.transpose(a2)
            dic = []
            if combo.startswith("DIC"):
                for vala, Canta, Cant2002a in zip(a, Cant, Cant2002):
                    if math.isnan(vala): 
                        dic.append("nan")
                    else:
                        dic.append(vala + Canta - Cant2002a)
                        
            else:
                dic = a2
        
            Cant_adjusted[combo] = dic
    else:
        Cant = [0] * len(Dates)
        Cant2002 = [0] * len(Dates)
   
    return Cant_adjusted, Cant, Cant2002
