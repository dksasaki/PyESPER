def defaults (DesiredVariables, OutputCoordinates={}, **kwargs):

    """
    Set default values and bookkeep inputs.
    """

    import numpy as np

    # Check and define Equations based on user-defined kwargs, or use default values
    Equations = kwargs.get("Equations", list(range(1, 17)))
    
    # Reading dimensions of user input
    n = max(len(v) for v in OutputCoordinates.values()) # number of rows out
    e = len(Equations) # number of Equations
    p = len(DesiredVariables) # number of Variables
                
    # Checking kwargs for presence of VerboseTF and EstDates, and Equations, and defining defaults, as needed
    VerboseTF = kwargs.get("VerboseTF", True)
        
    # Set EstDates based on kwargs, defaulting to 2002.0 if not provided
    if "EstDates" in kwargs:
        d = np.array(kwargs["EstDates"])
        EstDates = (
            [item for sublist in [kwargs["EstDates"]] * (n + 1) for item in sublist]
            if len(d) != n else list(d)
        )
    else:
        EstDates = [2002.0] * n
        
    # Bookkeeping coordinates
    C = {}
    longitude = np.array(OutputCoordinates["longitude"])
    longitude[longitude > 360] = np.remainder(longitude[longitude > 360], 360)
    longitude[longitude < 0] = longitude[longitude<0] + 360
    C["longitude"] = longitude
    C["latitude"] = OutputCoordinates["latitude"]
    C["depth"] = OutputCoordinates["depth"]   
    
    # Defining or reading in PerKgSwTF
    PerKgSwTF = kwargs.get("PerKgSwTF", True)

    # Defining Measurement Uncertainties
    MeasUncerts = kwargs.get("MeasUncerts", {})

    # Validate MeasUncerts dimensions
    if MeasUncerts:
        if max(len(v) for v in MeasUncerts.values()) != n:
            if min(len(v) for v in MeasUncerts.values()) != 1:
                raise CustomError(
                    "MeasUncerts must be undefined, a vector with the same number of elements as "
                    "PredictorMeasurements has columns, or a matrix of identical dimension to PredictorMeasurements."
                )
        if len(MeasUncerts) != len(PredictorMeasurements):
            print("Warning: Different numbers of input uncertainties and input measurements.")

    return Equations, n, e, p, VerboseTF, EstDates, C, PerKgSwTF, MeasUncerts
