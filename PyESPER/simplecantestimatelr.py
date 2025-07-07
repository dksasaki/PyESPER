def simplecantestimatelr(EstDates, longitude, latitude, depth):

    """
    Simple estimate of contribution of anthropogenic carbon to pH and DIC estimates.
    """

    import numpy as np
    import math
    import pandas as pd
    from scipy.interpolate import griddata

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
