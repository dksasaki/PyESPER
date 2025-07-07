def interpolate(Gdf={}, AAdata={}, Elsedata={}):

    """
    This LIR function performs the interpolation
    """

    import numpy as np
    from scipy.spatial import Delaunay
    import scipy.interpolate
  
    Gkeys, Gvalues = list(Gdf.keys()), list(Gdf.values())
    AAOkeys, AAOvalues, ElseOkeys, ElseOvalues = list(AAdata.keys()), list(AAdata.values()), list(Elsedata.keys()), \
        list(Elsedata.values())
       
    def process_grid(grid_values, data_values):
        results = []
        for i in range(len(grid_values)):
            grid = grid_values[i]
            points = np.array([list(grid['lon']), list(grid['lat']), list(grid['d2d'])]).T
            tri = Delaunay(points)
        
            values = np.array([
                list(grid['C_alpha']),
                list(grid['C_S']),
                list(grid['C_T']),
                list(grid['C_A']),
                list(grid['C_B']),
                list(grid['C_C'])
            ]).T
            interpolant = scipy.interpolate.LinearNDInterpolator(tri, values) 
                        
            data = data_values[i]
            points_to_interpolate = (list(data['Longitude']), list(data['Latitude']), list(data['d2d']))
            results.append(interpolant(points_to_interpolate))

        return results, interpolant

    # Process AA and EL grids
    aaLCs, aaInterpolants_pre = process_grid(Gvalues, AAOvalues)
    elLCs, elInterpolants_pre = process_grid(Gvalues, ElseOvalues)   
   
    return aaLCs, aaInterpolants_pre, elLCs, elInterpolants_pre
