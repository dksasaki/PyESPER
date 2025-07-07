def coefs_AAinds(Equations, LIR_data):

    """
    Separates coefficients from MATLAB ESPERv1 into Atlantic and Arctic or other regions.
    """

    import pandas as pd
    import numpy as np
    
    # Use boolean for AA or Else to separate coefficients into Atlantic or not
    GridCoords, Cs, AAInds = LIR_data[:3]    
    DVs, CsVs = list(Cs.keys()), list(Cs.values())   
    ListVars, NumVars = list(range(len(AAInds))), len(AAInds)
    GridValues, AAIndValues = list(GridCoords.values())[0], list(AAInds.values())[0]
    
    lon_grid, lat_grid, d2d_grid, aainds = np.array((GridValues[0])), np.array((GridValues[1])), \
       np.array(GridValues[2])/25, np.array(AAIndValues[0])
    names = ['lon', 'lat', 'd2d', "C_alpha", "C_S", "C_T", "C_A", "C_B", "C_C", 'AAInds']

    Gdf, CsDesired = {}, {}   
    for l, name in zip(ListVars, DVs):
        Cs2 = CsVs[:][l][:]
        for e in Equations:
            CsName = f'Cs{name}{e}'
            CsDesired[CsName] = Cs2[e-1][:]
            Cs3 = Cs2[e-1][:]
            C_alpha, C_S, C_T, C_A, C_B, C_C = np.array(Cs3[0]), np.array(Cs3[1]), np.array(Cs3[2]), np.array(Cs3[3]), \
                np.array(Cs3[4]), np.array(Cs3[5])
            grid_indices = np.column_stack((lon_grid, lat_grid, d2d_grid, C_alpha, C_S, C_T, C_A, C_B, C_C, aainds))
            Gdf[f"{name}{e}"] = pd.DataFrame(data=grid_indices, columns=names)

    return Gdf, CsDesired

