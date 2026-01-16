LinkeTurbidities.h5 - Global Linke Turbidity Climatology
=========================================================

Source: pvlib-python (https://github.com/pvlib/pvlib-python)
License: BSD 3-Clause (compatible with pvflux GPL-3)
Downloaded from: https://github.com/pvlib/pvlib-python/blob/main/pvlib/data/LinkeTurbidities.h5

Data Description
----------------
This file contains monthly climatological Linke turbidity values on a global
5 arcminute resolution grid (2160 × 4320 × 12 matrix).

- Coverage: Global (90°S to 90°N, -180° to 180°)
- Resolution: 5 arcminutes (0.0833°)
- Temporal: Monthly climatology (January-December)
- Format: HDF5 dataset named "LinkeTurbidity"
- Encoding: Values stored as uint8, multiply by 20 to get actual turbidity
- Size: ~15.6 MB

The Linke turbidity coefficient quantifies atmospheric turbidity (aerosols
and water vapor) relative to a clean, dry atmosphere. Typical values:
- 2.0: Very clean, clear sky
- 3.0: Clean, clear sky (rural)
- 4.0: Moderately turbid (urban)
- 5.0: Turbid
- 6-7: Very turbid (polluted/humid)

Original Data Source
--------------------
SoDa (Solar radiation Data) service
Based on the worldwide Linke turbidity database

Citation
--------
Remund, J., Wald, L., Lefèvre, M., Ranchin, T., & Page, J. (2003).
Worldwide Linke turbidity information. Proceedings of the ISES Solar
World Congress, Göteborg, Sweden.

Additional reference (meteonorm):
Section 8 of "Aerosol optical depth and Linke turbidity climatology"
Available from: https://meteonorm.com

Usage in pvflux
---------------
Use the lookup_linke_turbidity() function to retrieve turbidity values
for any location and time:

  library(pvflux)
  time <- seq(as.POSIXct("2026-01-01", tz = "UTC"),
              by = "month", length.out = 12)
  tl <- lookup_linke_turbidity(time, lat = -30.6279, lon = 24.0054)

License Compatibility
---------------------
The LinkeTurbidities.h5 file is distributed with pvlib-python under the
BSD 3-Clause License, which is permissive and compatible with pvflux's
GPL-3 license. The file is redistributed here under the same BSD 3-Clause
License terms. See pvlib-python LICENSE at:
https://github.com/pvlib/pvlib-python/blob/main/LICENSE
