def emlr_nn(Path, DesiredVariables, Equations, OutputCoordinates={}, PredictorMeasurements={}, **kwargs):

    """
    Estimating EMLR for nn's
    """

    from PyESPER.fetch_polys_NN import fetch_polys_NN
    import pandas as pd
    import numpy as np
    from scipy.interpolate import griddata

    EMLR = []
    for dv in DesiredVariables:
        DV = f"{dv}"
        print(DV)
        NN_data = fetch_polys_NN(Path, [DV])
        data_arrays = [
            np.nan_to_num(np.array([
                NN_data[1][i][c][b][a]
                for a in range(16)
                for b in range(11)
                for c in range(8)
            ]))
            for i in range(4)
        ]
        
        # Create DataFrame with meaningful column names
        UGridArray = pd.DataFrame({
            'UDepth': data_arrays[0],
            'USal': data_arrays[1],
            'Eqn': data_arrays[2],
            'RMSE': data_arrays[3],
        })
        
        UGridPoints = (UGridArray['UDepth'], UGridArray['USal'], UGridArray['Eqn'])
        UGridValues = UGridArray['RMSE']
        
        no_equations = len(Equations)
        # Perform estimation for each equation
        EM = [
            griddata(
                UGridPoints, UGridValues,
                (OutputCoordinates['depth'], PredictorMeasurements['salinity'], [Equations[eq]] * len(PredictorMeasurements['salinity'])),
                method='linear'
            )
            for eq in range(no_equations)
        ]
      
        EMLR.append(EM)
    
    return EMLR

