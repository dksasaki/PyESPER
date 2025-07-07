def temperature_define(
    DesiredVariables, 
    PredictorMeasurements, 
    InputAll, 
    **kwargs
):

    """ 
    A small function to define temperature as needed and adjust 
        InputAll and PredictorMeasurements acccordingly.
    """

    # Printing a custom warning if temperature is absent but needed 
    if "EstDates" in kwargs and "pH" in DesiredVariables:
        if "temperature" not in PredictorMeasurements:
            print(
                "Warning: Carbonate system calculations will be used to adjust the pH, but no temperature is "
                "specified so 10 C will be assumed. If this is a poor estimate for your region, consider supplying "
                "your own value in the PredictorMeasurements input."
            )
            Temperature = [10] * n
        else:
            Temperature = InputAll["Temperature"]
    
        PredictorMeasurements["temperature"] = Temperature
        InputAll["temperature"] = Temperature 

    return PredictorMeasurements, InputAll

    
