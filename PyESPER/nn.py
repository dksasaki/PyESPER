def nn(DesiredVariables, Path, OutputCoordinates={}, PredictorMeasurements={}, **kwargs):

    """
    Neural networks for seawater property estimation as part of PyESPERsv1.0.0
    """

    import time
    import pandas as pd
    from PyESPER.errors import errors
    from PyESPER.defaults import defaults
    from PyESPER.lir_uncertainties import measurement_uncertainty_defaults
    from PyESPER.inputdata_organize import inputdata_organize
    from PyESPER.temperature_define import temperature_define
    from PyESPER.iterations import iterations
    from PyESPER.fetch_polys_NN import fetch_polys_NN
    from PyESPER.define_polygons import define_polygons
    from PyESPER.run_nets import run_nets
    from PyESPER.process_netresults import process_netresults
    from PyESPER.emlr_nn import emlr_nn
    from PyESPER.organize_nn_output import organize_nn_output
    from PyESPER.pH_DIC_nn_adjustment import pH_DIC_nn_adjustment
    from PyESPER.final_formatting import final_formatting

    # Starting the timer 
    tic = time.perf_counter()
    
    # Function that provides custom error messages for erroneous inputs
    errors(OutputCoordinates, PredictorMeasurements)

    # Function which calculates default measurement uncertainties
    Equations, n, e, p, VerboseTF, EstDates, C, PerKgSwTF, MeasUncerts = defaults(
        DesiredVariables,
        OutputCoordinates,
        **kwargs
    )

    # Function that processes the input values and default uncertainties and makes sense of it
    Uncertainties_pre, DUncertainties_pre = measurement_uncertainty_defaults(
        n,
        PredictorMeasurements,
        MeasUncerts
    )

    # Creating a pandas DataFrame of input data
    InputAll = inputdata_organize(
        EstDates,
        C,
        PredictorMeasurements,
        Uncertainties_pre
    )

    PredictorMeasurements, InputAll = temperature_define(
        DesiredVariables, 
        PredictorMeasurements,
        InputAll,
        **kwargs
    )

    code, unc_combo_dict, dunc_combo_dict = iterations(
        DesiredVariables,
        Equations,
        PerKgSwTF,
        C,
        PredictorMeasurements,
        InputAll,
>>>>>>> 591e951 (added examples.py)
        Uncertainties_pre,
        DUncertainties_pre
    )

<<<<<<< HEAD
    for d, var in enumerate(DesiredVariables):
        Pertdv, DPertdv, Unc, DUnc = [], [], [], []
        var = [var]  # Wrap single variable in a list
        keys = ["sal_u", "temp_u", "phosphate_u", "nitrate_u", "silicate_u", "oxygen_u"]

        PredictorMeasurements2, Est, Uncertainties, DUncertainties, emlr = preprocess_applynets(
            var, 
            Equations, 
            EstDates, 
            ['False'], 
            OutputCoordinates, 
            PredictorMeasurements, 
            Uncertainties_pre, 
            DUncertainties_pre
        )
        
        names = list(PredictorMeasurements2.keys())
        PMs = list(PredictorMeasurements2.values())

        # Replace "nan" with 0 in PMs using list comprehensions
        PMs_nonan = [[0 if val == "nan" else val for val in pm] for pm in PMs]

        # Transpose PMs_nonan
        PMs = np.transpose(PMs_nonan)

        PMs3, DMs3 = {}, {}

        for pred in range(len(PredictorMeasurements2)):
            num_coords = len(OutputCoordinates['longitude'])
            num_preds = len(PredictorMeasurements2)

            # Initialize perturbation arrays
            Pert = np.zeros((num_coords, num_preds))
            DefaultPert = np.zeros((num_coords, num_preds))

            # Populate perturbation arrays
            Pert[:, pred] = Uncertainties[keys[pred]]
            DefaultPert[:, pred] = DUncertainties[keys[pred]]

            # Apply perturbations
            PMs2 = PMs + Pert
            DMs2 = PMs + DefaultPert

           # Update PMs3 and DMs3 dictionaries
            for col, name in enumerate(names):
                PMs3[name] = PMs2[:, col].tolist()
                DMs3[name] = DMs2[:, col].tolist()

            # Run preprocess_applynets for perturbed and default data
            VTF = False
            _, PertEst, _, _, _ = preprocess_applynets(
                var, Equations, EstDates, VTF, OutputCoordinates, PMs3, Uncertainties_pre, DUncertainties_pre
            )
            _, DefaultPertEst, _, _, _ = preprocess_applynets(
                var, Equations, EstDates, VTF, OutputCoordinates, DMs3, Uncertainties_pre, DUncertainties_pre
            )

            # Extract estimates and perturbation results
            combo, estimates = list(Est.keys()), list(Est.values())
            pertests, defaultpertests = list(PertEst.values()), list(DefaultPertEst.values())

            # Initialize result lists
            PertDiff, DefaultPertDiff, Unc_sub2, DUnc_sub2 = [], [], [], []
        
            for c in range(len(Equations)):
                # Compute differences and squared differences using list comprehensions
                PD = [estimates[c][e] - pertests[c][e] for e in range(len(estimates[c]))]
                DPD = [estimates[c][e] - defaultpertests[c][e] for e in range(len(estimates[c]))]
                Unc_sub1 = [(estimates[c][e] - pertests[c][e])**2 for e in range(len(estimates[c]))]
                DUnc_sub1 = [(estimates[c][e] - defaultpertests[c][e])**2 for e in range(len(estimates[c]))]

                # Append results to their respective lists
                PertDiff.append(PD)
                DefaultPertDiff.append(DPD)
                Unc_sub2.append(Unc_sub1)
                DUnc_sub2.append(DUnc_sub1)
            Pertdv.append(PertDiff)
            DPertdv.append(DefaultPertDiff)
            Unc.append(Unc_sub2)
            DUnc.append(DUnc_sub2)
        PD_final.append(Pertdv)
        DPD_final.append(DPertdv)
        Unc_final.append(Unc)
        DUnc_final.append(DUnc) # CHECK THIS WHOle shebang next

    est = list(Est_pre.values())
    Uncertainties = []
    propu = []
    for dv in range(0, len(DesiredVariables)):
        dvu = []
        for eq in range(0, len(Equations)):
            sumu = []
            for n in range(0, len(est[0])):
                u, du = [], []
                for pre in range(0, len(PredictorMeasurements)):
                    u.append(Unc_final[dv][pre][eq][n])
                    du.append(DUnc_final[dv][pre][eq][n])
                eu = EMLR[dv][eq][n]
                sumu.append((sum(u) - sum(du) + eu**2)**(1/2))
            dvu.append(sumu)
        propu.append(dvu)
    Uncertainties.append(propu)
    YouHaveBeenWarnedCanth = False


    def SimpleCantEstimateLR(EstDates, longitude, latitude, depth):
        # Load interpolation points and values
        CantIntPoints = pd.read_csv('SimpleCantEstimateLR_full.csv')
        pointsi = (
            CantIntPoints['Int_long'] * 0.25,
            CantIntPoints['Int_lat'],
            CantIntPoints['Int_depth'] * 0.025,
        )
        values = CantIntPoints['values']

        # Scale input coordinates
        pointso = (
            np.array(longitude) * 0.25,
            np.array(latitude),
            np.array(depth) * 0.025,
        )

        # Interpolate and compute Cant2002
        Cant2002 = griddata(pointsi, values, pointso, method='linear')

        # Adjust for estimation dates
        CantMeas = [
            c * math.exp(0.018989 * (date - 2002)) for c, date in zip(Cant2002, EstDates)
        ]

        return CantMeas, Cant2002

    Cant_adjusted = {}
    combos2 = list(Est_pre.keys())
    values2 = list(Est_pre.values())

    if "EstDates" in kwargs and ("DIC" in DesiredVariables or "pH" in DesiredVariables):      
        if not YouHaveBeenWarnedCanth:
            if VerboseTF:
                print("Estimating anthropogenic carbon for PyESPER_NN.")
            longitude = np.mod(OutputCoordinates["longitude"], 360)
            latitude = np.array(OutputCoordinates["latitude"])
            depth = np.array(OutputCoordinates["depth"])
            Cant, Cant2002 = SimpleCantEstimateLR(EstDates, longitude, latitude, depth)
            YouHaveBeenWarnedCanth = True

        for combo, a in zip(combos2, values2):
            dic = []
            if combo.startswith("DIC"):
                for vala, Canta, Cant2002a in zip(a, Cant, Cant2002):
                    if math.isnan(vala): 
                        dic.append("nan")
                    else:
                        dic.append(vala + Canta - Cant2002a)
            else:
                dic = list(a)
            Cant_adjusted[combo] = dic
                 
        if "pH" in DesiredVariables:
            warning = []
            for combo, values in zip(combos2, values2):
                if combo.startswith("pH"):
                    salinity = PredictorMeasurements["salinity"]
                    PM_pH = {'salinity': salinity}
                    eq = [16]
                    alkpm, alkest, _, _, _ = preprocess_applynets(
                        ["TA"], eq, EstDates, ['False'], C, PM_pH, Uncertainties_pre, DUncertainties_pre
                    )
                    EstAlk = np.array(alkest["TA16"])
                    EstSi = EstP = [0] * len(EstAlk)
                    Pressure = sw.pres(OutputCoordinates["depth"], OutputCoordinates["latitude"])
                    Est = np.array(values)

                    # CO2SYS calculations
                    Out = pyco2.sys(
                        par1=EstAlk, par2=Est, par1_type=1, par2_type=3, salinity=salinity,
                        temperature=PredictorMeasurements["temperature"], temperature_out=PredictorMeasurements["temperature"],
                        pressure=Pressure, pressure_out=Pressure, total_silicate=EstSi, total_phosphate=EstP, opt_total_borate=2
                    )
                    DICadj = Out["dic"] + Cant - Cant2002
                    OutAdj = pyco2.sys(
                        par1=EstAlk, par2=DICadj, par1_type=1, par2_type=2, salinity=salinity,
                        temperature=PredictorMeasurements["temperature"], temperature_out=PredictorMeasurements["temperature"],
                        pressure=Pressure, pressure_out=Pressure, total_silicate=EstSi, total_phosphate=EstP, opt_total_borate=2
                    )
                    pHadj = OutAdj["pH"]

                    # Check for convergence warnings
                    if any(np.isnan(pHadj)):
                        warning_message = (
                            "Warning: CO2SYS took >20 iterations to converge. The corresponding estimate(s) will be NaN. "
                            "This typically happens when ESPER_NN is poorly suited for estimating water with the given properties "
                            "(e.g., very high or low salinity or estimates in marginal seas)."
                        )
                        warning.append(warning_message)
                else:
                    pHadj = np.array(values)

                Cant_adjusted[combo] = pHadj.tolist()

            # Print warnings if any
            if warning:
                print(warning[0])

    elif "EstDates" not in kwargs and ("DIC" or "pH" in DesiredVariables) and VerboseTF and not YouHaveBeenWarnedCanth:
        print("Warning: DIC or pH is a requested output but the user did not provide dates for the desired estimates.  The estimates will be specific to 2002.0 unless the optional EstDates input is provided (recommended).")
        YouHaveBeenWarnedCanth = True

    if kwargs.get("pHCalcTF") and "pH" in DesiredVariables:
        if VerboseTF:
            print("Recalculating the pH to be appropriate for pH values calculated from TA and DIC.")
        for combo, pH_values in zip(combos2, values2):
            if combo.startswith("pH"):
                pH_adjcalc_Est = [(pH + 0.3168) / 1.0404 for pH in pH_values]
                Cant_adjusted[combo] = pH_adjcalc_Est

    # Prepare data for processing
    combos3 = Cant_adjusted.keys()
    values3 = Cant_adjusted.values()
    Us = Uncertainties[0]
    Us2 = [u2 for u in Us for u2 in u]

    # Convert combos and values to lists for iteration
    k2, v2 = list(combos2), list(values2)
    k3, v3 = list(combos3), list(values3)

    # Initialize estimates and uncertainties dictionaries
    Estimates, Uncerts = {}, {}

    for key2, value2 in zip(k2, v2):
        # Adjust values in v2 based on matches in k3
        adjusted_array = np.array(value2)
        for key3, value3 in zip(k3, v3):
            adjusted_array[key2 == key3] = value3

        # Store adjusted values and uncertainties
        Estimates[key2] = adjusted_array
        Uncerts[key2] = Us2[k2.index(key2)]

    # Convert results to DataFrame
    Estimates = pd.DataFrame(Estimates)
    Uncertainties = pd.DataFrame(Uncerts)

    toc = time.perf_counter()
    print(f"PyESPER_NN took {toc - tic:0.4f} seconds, or {(toc-tic)/60:0.4f} minutes to run")    
=======
    NN_data = fetch_polys_NN(
        Path, 
        DesiredVariables
    )

    df = define_polygons(C)

    EstAtl, EstOther = run_nets(
        DesiredVariables, 
        Equations, 
        code
    )

    Estimates, no_equations = process_netresults(
        Equations, 
        code, 
        df, 
        EstAtl, 
        EstOther
    )

    EMLR = emlr_nn(
        Path, 
        DesiredVariables,  
        Equations,
        OutputCoordinates, 
        PredictorMeasurements, 
    )

    Uncertainties = organize_nn_output(
        Path,
        DesiredVariables,
        OutputCoordinates,
        PredictorMeasurements,
        **kwargs
    )

    YouHaveBeenWarnedCanth=False
    Cant_adjusted = pH_DIC_nn_adjustment(
        Path,
        DesiredVariables, 
        Estimates,
        YouHaveBeenWarnedCanth,
        OutputCoordinates,
        PredictorMeasurements,
        **kwargs
    )
    
    Estimates = final_formatting(
        DesiredVariables,
        Cant_adjusted,
        Estimates
    )
       
    toc = time.perf_counter()
    print(f"PyESPER_NN took {toc - tic:0.4f} seconds, or {(toc-tic)/60:0.4f} minutes to run")
>>>>>>> 591e951 (added examples.py)

    return Estimates, Uncertainties
