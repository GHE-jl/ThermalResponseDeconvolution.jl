"""
Importer that builds a `TRTData` structure from a CSV file.
"""

using CSV, DataFrames

"""
    TRTData_import(filepath)

Read a `.csv` file containing four columns and a single header line, then store
the data in a `TRTData` structure. The first time step should be at time 0 with
no heating; the following rows describe the test with the heater on:

  1. Time steps [s]
  2. Heating power [W]
  3. Inlet fluid temperature [degC]
  4. Outlet fluid temperature [degC]

`Texp` is computed as `Tout - Tout[1]`. A time of exactly 0 at the first sample
is nudged to `1e-4` to avoid a singular first time step.
"""
function TRTData_import(filepath::String)
    DATA = CSV.read(filepath, DataFrame)

    data = TRTData(
        Float64.(DATA[:, 1]),
        Float64.(DATA[:, 2]),
        Float64.(DATA[:, 3]),
        Float64.(DATA[:, 4]),
        Float64.(DATA[:, 4]) .- DATA[1, 4],
    )

    # Avoid a singular first time step
    data.t[1] = data.t[1] == 0 ? 1e-4 : data.t[1]
    return data
end
