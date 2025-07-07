def fetch_data (DesiredVariables, Path):
    
    """
    Gathers the necessary LIR files similar to ESPERs
    """

    from scipy.io import loadmat
    import os
    import numpy as np
    import pandas as pd

    AAIndsCs, GridCoords, Cs = {}, {}, {}

    for v in DesiredVariables:               
        fname1 = os.path.join(Path, f"Mat_fullgrid/LIR_files_{v}_fullCs1.mat")
        fname2 = os.path.join(Path, f"Mat_fullgrid/LIR_files_{v}_fullCs2.mat")
        fname3 = os.path.join(Path, f"Mat_fullgrid/LIR_files_{v}_fullCs3.mat")
        fname4 = os.path.join(Path, f"Mat_fullgrid/LIR_files_{v}_fullGrids.mat")

        Cs1 = loadmat(fname1)
        Cs2 = loadmat(fname2)
        Cs3 = loadmat(fname3)
        Grid = loadmat(fname4)

        UncGrid = Grid["UncGrid"][0][0] 
        GridCoodata, AAInds = np.array(Grid["GridCoords"]), np.array(Grid["AAIndsM"]) 
        Csdata1, Csdata2, Csdata3 = np.array(Cs1["Cs1"]), np.array(Cs2["Cs2"]), np.array(Cs3["Cs3"])        
        AAIndsCs[v] = pd.DataFrame(data=AAInds)
        GridCoords[v] = pd.DataFrame(data=GridCoodata[:, :])
        Csdata = np.concatenate((Csdata1, Csdata2, Csdata3), axis=1)
        Cs[v] = [pd.DataFrame(data=Csdata[:, :, i]) for i in range(16)]
      
    LIR_data = [GridCoords, Cs, AAIndsCs, UncGrid]
    return LIR_data
