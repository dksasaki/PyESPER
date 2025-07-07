def final_formatting(DesiredVariables, Cant_adjusted={}, Est_pre={}):

    """
    Formatting the final bits.
    """
    
    import pandas as pd
    
    if ("pH" or "DIC") in DesiredVariables:
        Estimates=pd.DataFrame(Cant_adjusted)
        print("anthropogenic C has been incorporated into some estimates")
    else:
        Estimates=Est_pre
        print("anthropogenic carbon is not considered for these estimates")        

    return Estimates
