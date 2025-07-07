def process_netresults(Equations, code={}, df={}, EstAtl={}, EstOther={}):

    """ 
    Regional smoothing and processing net outputs
    """

    from statistics import mean
    import pandas as pd
    import numpy as np

    def process_estimates(estimates):
        keys = list(estimates.keys())
        values = list(estimates.values())
        result = {}
        for i, key in enumerate(keys):
            result[key] = [mean(values[i][v]) for v in range(len(values[0]))]
        return result

    Esta = process_estimates(EstAtl)
    Esto = process_estimates(EstOther)
    
    # Processing regionally in the Atlantic and Bering
    EstA, EstB, EB2, ESat, ESat2, ESaf, Estimate  = {}, {}, {}, {}, {}, {}, {}

    for i in code:
        code[i]["AAInds"] = df["AAInds"]
        code[i]["BeringInds"] = df["BeringInds"]
        code[i]["SAtlInds"] = df["SAtlInds"]
        code[i]["SoAfrInds"] = df["SoAfrInds"]

    for codename, codedata in code.items():
        Estatl, Estb, eb2, Estsat, esat2, esafr, esaf2 = [], [], [], [], [], [], []
        aainds, beringinds, satlinds, latitude, safrinds = (
            codedata[key] for key in ["AAInds", "BeringInds", "SAtlInds", "Latitude", "SoAfrInds"]
        )

        Estatl = [Esta[codename][i] if aa_ind else Esto[codename][i] for i, aa_ind in enumerate(aainds)]

        for l in range(0, len(Estatl)):
            repeated_values = (latitude[l]-62.5)/7.5
            B = np.tile(repeated_values, (1, len(Equations)))
            C = Esta[codename][l]
            B1 = C * B
            repeated_values2 = (70-latitude[l])/7.5
            D = np.tile(repeated_values2, (1, len(Equations)))
            E = Esto[codename][l]
            B2 = E * D
            Estb.append(B1[0][0] + B2[0][0])

        eb2 = [Estb[j] if b_ind else Estatl[j] for j, b_ind in enumerate(beringinds)]

        for n in range(0, len(satlinds)): 
            repeated_values = (latitude[n]+44)/10
            F1 = Esta[codename][n]
            F = np.tile(repeated_values, (1, len(Equations)))
            G1 = F1 * F
            repeated_values2 = (-34-latitude[n])/10
            H1 = Esto[codename][n]
            H = np.tile(repeated_values2, (1, len(Equations)))
            G2 = H1 * H
            Estsat.append(G1[0][0] + G2[0][0])
    
        EstA[codename], EstB[codename], EB2[codename], ESat[codename] = Estatl, Estb, eb2, Estsat

    # Regional processing for S. Atlantic
        ESat2[codename] = [
            ESat[codename][i] if satlinds[i] == "True" else EB2[codename][i]
            for i in range(len(satlinds))
        ]
            
    # Regional processing for S. Africa
        for s in range(0, len(safrinds)):
            repeated_values = (27-df["Lon"][s])/8
            F1 = ESat2[codename][s]
            F = np.tile(repeated_values, (1, len(Equations)))
            G1 = F1*F
            repeated_values2 = (df["Lon"][s]-19)/8
            H1 = Esto[codename][s]
            H = np.tile(repeated_values2, (1, len(Equations)))
            G2 = H1 * H
            esafr.append(G1[0][0] + G2[0][0])

        ESaf[codename] = esafr

        Estimate[codename] = [
            ESaf[codename][i] if safrinds[i] == "True" else ESat2[codename][i]
            for i in range(len(safrinds))
        ]

    # Bookkeeping blanks back to NaN as needed
    Estimate = {k: ('NaN' if v == '' else v) for k, v in Estimate.items()}
    no_equations = len(Equations)

    return Estimate, no_equations

