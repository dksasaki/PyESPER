def emlr_estimate(Equations, DesiredVariables, Path, OutputCoordinates={}, PredictorMeasurements={}, UDict={}, DUDict={}, Coefficients={}, **kwargs):
    
    """
    Uncertainty estimation step 1
    """
    
    import numpy as np
    import pandas as pd        
    from scipy.interpolate import griddata
    from PyESPER.fetch_data import fetch_data
    
    EMLR, varnames, EqM = {}, [], []
            
    for dv in DesiredVariables:
        # Fetch LIR data and process into grid arrays
        LIR_data = fetch_data([dv], Path)
        grid_names = ['UDepth', 'USal', 'Eqn', 'RMSE']
        UGridArray = pd.DataFrame(
            [np.nan_to_num([LIR_data[3][i][c][b][a] for a in range(16) for b in range(11) for c in range(8)]) for i in range(4)]
        ).T
        UGridArray.columns = grid_names
        UGridPoints, UGridValues = (UGridArray['UDepth'], UGridArray['USal'], UGridArray['Eqn']), UGridArray['RMSE']
        
        for eq in range(len(Equations)):
            varnames.append(dv + str(Equations[eq]))
            EM = []
            eq_str = str(Equations[eq])
            eq_repeated = [Equations[eq]] * len(PredictorMeasurements['salinity'])
            UGridPointsOut = (OutputCoordinates['depth'], PredictorMeasurements['salinity'], eq_repeated)
            emlr = griddata(UGridPoints, UGridValues, UGridPointsOut, method='linear')
            combo = f"{dv}{eq_str}"
            Coefs = {k: np.nan_to_num(np.array(Coefficients[combo][k])) for k in ["Intercept", "Coef S", "Coef T", "Coef A", "Coef B", "Coef C"]}
            uncdfs, duncdfs = UDict[combo], DUDict[combo]
            keys = uncdfs.columns.to_numpy()
            USu, UTu, UAu, UBu, UCu = [np.nan_to_num(uncdfs[key].fillna(0).astype(float)) for key in keys]
            DUSu, DUTu, DUAu, DUBu, DUCu = [np.nan_to_num(duncdfs[key].fillna(0).astype(float)) for key in keys]
            USu2, UTu2, UAu2, UBu2, UCu2 = [np.nan_to_num(uncdfs[key].fillna(-9999).astype(float)) for key in keys]
            DUSu2, DUTu2, DUAu2, DUBu2, DUCu2 = [np.nan_to_num(duncdfs[key].fillna(-9999).astype(float)) for key in keys]
            
            C0u2 = Coefs["Intercept"] * 0
            Csum, DCsum = [], []
                    
            for cucombo in range(len(Coefs["Coef S"])):
    
                s1 = (Coefs["Coef S"][cucombo]*USu[cucombo])**2
                t1 = (Coefs["Coef T"][cucombo]*UTu[cucombo])**2
                a1 = (Coefs["Coef A"][cucombo]*UAu[cucombo])**2
                b1 = (Coefs["Coef B"][cucombo]*UBu[cucombo])**2
                c1 = (Coefs["Coef C"][cucombo]*UCu[cucombo])**2
                sum2 = s1+t1+a1+b1+c1
                ds1 = (Coefs["Coef S"][cucombo]*DUSu[cucombo])**2
                dt1 = (Coefs["Coef T"][cucombo]*DUTu[cucombo])**2
                da1 = (Coefs["Coef A"][cucombo]*DUAu[cucombo])**2
                db1 = (Coefs["Coef B"][cucombo]*DUBu[cucombo])**2
                dc1 = (Coefs["Coef C"][cucombo]*DUCu[cucombo])**2
                dsum2 = ds1+dt1+da1+db1+dc1
                        
                uncestimate = (sum2 - dsum2 + emlr[cucombo]**2)**0.5
                EM.append(uncestimate)            
            EqM.append(EM)
            
    EqM2 = []
    for i in EqM:
        UncertEst = np.array(i)
        UncertEst = UncertEst.astype('float')
        UncertEst[USu2==-9999]=['nan']
        if Equations[eq] == 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8:
            UncertEst[UTu2==-9999]=['nan']
        if Equations[eq] == 1 | 2 | 5 | 6 | 9 | 10 | 13 | 14:
            UncertEst[UAu2==-9999]=['nan']
        if Equations[eq] == 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15:
            UncertEst[UBu2==-9999]=['nan']
        if Equations[eq] == 1 | 2 | 3 | 4 | 9 | 10 | 11 | 12:
            UncertEst[UCu2==-9999]=['nan']
        EqM2.append(UncertEst)
    for key in range(0, len(varnames)):
        EMLR[varnames[key]] = EqM2[key]
        EMLR = pd.DataFrame(EMLR)
              
    return EMLR

