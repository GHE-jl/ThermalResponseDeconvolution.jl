# References

The method implemented in this package is drawn from the following sources.

## Deconvolution method

- **Dion, G., Pasquier, P., & Marcotte, D.** (2022). Deconvolution of experimental thermal
  response test data to recover short-term g-function. *Geothermics*, 100, 102302.
  <https://doi.org/10.1016/j.geothermics.2021.102302>
  Source of the multi-objective deconvolution formulation, the node parameterization, and the
  exponential-integral initial guess (Eqs. 20–21).

- **Dion, G., Pasquier, P., Marcotte, D., & Beaudry, G.** (2023). Multi-deconvolution in
  non-stationary conditions applied to experimental thermal response test analysis to obtain
  short-term transfer functions. *Science and Technology for the Built Environment*, 30(3), 1–14.
  <https://doi.org/10.1080/23744731.2023.2217729>
  Source of the multi-objective weighting scheme (`SetWeigth`) balancing the temperature-misfit,
  first-, and second-derivative terms.

- **Dion, G., Pasquier, P., & Marcotte, D.** (2024). Application of deconvolution to
  interpretation of distributed thermal response test (DTRT) and to determination of thermal
  conductivity profiles. *Applied Thermal Engineering*, 236, 121680.
  <https://doi.org/10.1016/j.applthermaleng.2023.121680>
  Extension of the method to distributed thermal response tests.

## Convolution

- **Marcotte, D., & Pasquier, P.** (2008). Fast fluid and ground temperature computation for
  geothermal ground-loop heat exchanger systems. *Geothermics*, 37(6), 651–665.
  <https://doi.org/10.1016/j.geothermics.2008.08.003>
  Non-circular convolution scheme underlying [`convolution`](@ref).

- **Pasquier, P., & Marcotte, D.** (2013). Efficient computation of heat flux signals to ensure
  the reproduction of prescribed temperatures at several interacting heat sources. *Applied
  Thermal Engineering*, 59(1–2), 515–526.
  <https://doi.org/10.1016/j.applthermaleng.2013.06.018>
