"""
Script that constructs a structure to hold the data of a TRT. The first step is to import the data,
and then storing it in a inmutable structure.
"""

using CSV, DataFrames

function TRTData_import(filepath::String)
    """
        TRTData_import(path, file)
    
    Requires the file path and the file name to search for a .CSV file containing 4 columns ans 1
    line of header. The first time step should be at time 0 with no heating. Then, the first time
    should be with the heating power on:
        - Time steps [s]
        - Heating power [W]
        - Inlet fluid temperature [degC]
        - Outlet fluid temperature [degC]
    Then, the data are stored in the structure TRTData, which is immutable.
    """
    # Read the .csv file
    DATA = CSV.read(filepath, DataFrame)   

    # Write the structure
    data = TRTData(DATA[:, 1], DATA[:, 2], DATA[:, 3], DATA[:, 4], DATA[:, 4] .- DATA[1, 4])
    
    # Correct the time array
    data.t[1] = data.t[1] == 0 ? 1e-4 : data.t[1]
    return data
end