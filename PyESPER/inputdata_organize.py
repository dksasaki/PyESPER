def inputdata_organize(
    EstDates, 
    C={}, 
    PredictorMeasurements={}, 
    Uncertainties={} 
):

    """
    This function preprocesses data into a pandas DataFrame
    """
  
    import pandas as pd

    n = max(len(v) for v in C.values()) # number of rows out

    # Redefining and organizing all data thus far
    order = list(range(n))
    input_data = {
        "Order": order,
        "Longitude": C["longitude"],
        "Latitude": C["latitude"],
        "Depth": C["depth"],
        "Salinity": PredictorMeasurements["salinity"],
        "Dates": EstDates,
        "Salinity_u": Uncertainties["sal_u"],
        "Temperature_u": Uncertainties["temp_u"],
        "Phosphate_u": Uncertainties["phosphate_u"],
        "Nitrate_u": Uncertainties["nitrate_u"],
        "Silicate_u": Uncertainties["silicate_u"],
        "Oxygen_u": Uncertainties["oxygen_u"]
    }

    # Map PredictorMeasurements keys to input_data keys
    for key, label in {
        "temperature": "Temperature",
        "phosphate": "Phosphate",
        "nitrate": "Nitrate",
        "silicate": "Silicate",
        "oxygen": "Oxygen"
    }.items():
        if key in PredictorMeasurements:
            input_data[label] = PredictorMeasurements[key]
   
    # Create a DataFrame with order stamp and drop all NaNs from a replicate dataframe
    InputAll = pd.DataFrame(input_data)

    return InputAll



