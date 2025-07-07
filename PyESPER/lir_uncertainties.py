def process_uncertainties(param, default_factor, MeasUncerts, PredictorMeasurements, n):
   
    """ 
    Helps proces uncertainties as needed.
    """
    
    import numpy as np

    if param in MeasUncerts:
        result = np.array(MeasUncerts.get(param))
        if len(result) < n:
            result = np.tile(result, n)
        if param.replace('_u', '') in PredictorMeasurements:
            dresult = np.array([i * default_factor for i in PredictorMeasurements[param.replace('_u', '')]])
        else:
            dresult = result
    else:
        if param.replace('_u', '') in PredictorMeasurements:
            result = np.array([i * default_factor for i in PredictorMeasurements[param.replace('_u', '')]])
            dresult = result
        else:
            result = np.tile('nan', n)
            dresult = np.tile(0, n)
    return result, dresult

def measurement_uncertainty_defaults(n, PredictorMeasurements={}, MeasUncerts={}):
    
    """
    Defines the default measurement uncertainties, as needed
    """
    
    import numpy as np
    import pandas as pd

    MeasUncerts_processed, DefaultUAll = {}, {}

    # Default salinity uncertainties
    sal_u = np.array(MeasUncerts.get("sal_u", [0.003]))
    sal_u = np.tile(sal_u, n) if len(sal_u) < n else sal_u
    sal_defu = np.tile(0.003, n) 

    # Temperature uncertainties
    temp_u = np.tile(np.array(MeasUncerts.get("temp_u", [0.003])), n) if "temp_u" in MeasUncerts or "temperature" in PredictorMeasurements else np.tile("nan", n)
    temp_defu = np.tile(0.003 if "temp_u" in MeasUncerts or "temperature" in PredictorMeasurements else 0, n)

    # Process other parameters
    parameters = {
        "phosphate_u": 0.02,
        "nitrate_u": 0.02,
        "silicate_u": 0.02,
        "oxygen_u": 0.01
    }

    for param, factor in parameters.items():
        MeasUncerts_processed[param], DefaultUAll[f"{param.replace('_u', '_defu')}"] = process_uncertainties(
            param, factor, MeasUncerts, PredictorMeasurements, n
        )

    # Update MeasUncerts and DefaultUAll dictionaries
    meas_uncerts_keys = ["sal_u", "temp_u", *parameters.keys()]
    default_uall_keys = ["sal_defu", "temp_defu", *[k.replace('_u', '_defu') for k in parameters.keys()]]

    MeasUncerts.update(dict(zip(meas_uncerts_keys, [sal_u, temp_u, *MeasUncerts_processed.values()])))
    DefaultUAll.update(dict(zip(default_uall_keys, [sal_defu, temp_defu, *DefaultUAll.values()])))

    # Create DataFrames
    keys = meas_uncerts_keys
 
    Uncerts = np.column_stack([MeasUncerts[k] for k in keys])
    Uncertainties_pre = pd.DataFrame(Uncerts, columns=keys)

    DUncerts = np.column_stack([DefaultUAll[k] for k in default_uall_keys])
    DUncertainties_pre = pd.DataFrame(DUncerts, columns=keys)

    return Uncertainties_pre, DUncertainties_pre
