# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

pvflux is an R package for solar PV power forecasting with a modular pipeline architecture. Users can independently select:
- Transposition models (Hay-Davies, Reindl, Perez, Olmo)
- Decomposition models (Erbs, Boland-Ridley)
- Cell temperature models (Skoplaki, Faiman)
- Clear-sky models (Ineichen-Perez, Haurwitz)

Run ensemble analysis across all model combinations, or estimate clear-sky power for performance monitoring.

**Important:** The Olmo transposition model has known validation issues (RMSE of 21-52%) outside its calibration region (Granada, Spain). Default is now "haydavies" - see `?olmo_transposition` for details.

## Architecture

The package uses a **three-level abstraction** with two parallel pathways:

### Main Pathway (Measured Irradiance)
1. **Individual model functions** - Atomic components for each physical model
2. **Modular pipeline** - `pv_dc_pipeline()` allows independent model selection
3. **Convenience functions** - `pv_power_pipeline()` (DC+AC), `pv_power_ensemble()` (all combinations)

### Clear-Sky Pathway (Estimated Irradiance)
1. **Clear-sky models** - `ineichen_clearsky()` (solar geometry + turbidity) or `haurwitz_clearsky()` (solar geometry only)
2. **Clear-sky pipelines** - `pv_clearsky_dc_pipeline()` and `pv_clearsky_power_pipeline()` with model selection
3. **Performance metrics** - `clearsky_index()` and `clearsky_performance_ratio()` for monitoring

### Data Flow Pipeline (Main Pathway)

```
Input (time, lat, lon, GHI, T_air, wind, tilt, azimuth)
    ↓
Step 1: Decomposition (GHI → DNI/DHI) [if needed]
    ├─ erbs_decomposition()
    └─ boland_decomposition()
    ↓
Step 2: Transposition (GHI/DNI/DHI → G_poa)
    ├─ olmo_transposition() [direct GHI → POA]
    ├─ haydavies_transposition() [uses DNI/DHI]
    ├─ reindl_transposition() [uses DNI/DHI]
    └─ perez_transposition() [uses DNI/DHI]
    ↓
Step 3: Cell Temperature (G_poa, T_air, wind → T_cell)
    ├─ skoplaki_cell_temperature()
    └─ faiman_cell_temperature()
    ↓
Step 4: DC Power (G_poa, T_cell → P_dc)
    └─ pvwatts_dc()
    ↓
Step 5: AC Conversion (P_dc → P_ac)
    └─ pv_ac_simple_clipping()
```

### Data Flow Pipeline (Clear-Sky Pathway)

```
Input (time, lat, lon, T_air, wind, tilt, azimuth, linke_turbidity, altitude)
    ↓
Step 1: Clear-Sky Irradiance
    ├─ ineichen_clearsky() → (GHI, DNI, DHI, airmass) [requires turbidity]
    └─ haurwitz_clearsky() → (GHI, DNI, DHI, airmass) [simpler, no turbidity]
    ↓
Step 2-5: Same as Main Pathway
    └─ Uses GHI from clear-sky model instead of measurements
```

### Key Model Combinations for Ensemble Analysis

`pv_power_ensemble()` and `pv_dc_ensemble()` run all 8 combinations:

| Identifier | Transposition | Cell Temperature |
|------------|---------------|------------------|
| haydavies_skoplaki | Hay-Davies | Skoplaki (default) |
| haydavies_faiman | Hay-Davies | Faiman |
| reindl_skoplaki | Reindl | Skoplaki |
| reindl_faiman | Reindl | Faiman |
| perez_skoplaki | Perez | Skoplaki |
| perez_faiman | Perez | Faiman |
| olmo_skoplaki | Olmo | Skoplaki |
| olmo_faiman | Olmo | Faiman |

**Ensemble utility functions** in `R/ensemble_utils.R`:
- `ensemble_summary()` - Calculate mean, SD, min, max across models
- `ensemble_wide()` - Convert long format to wide (one column per model)
- `ensemble_rank()` - Rank models by performance metric
- `pv_spread()` - Calculate ensemble spread (max - min)
- `ensemble_plot_data()` - Prepare data for plotting with uncertainty bands

## Common Development Commands

```bash
# Generate/refresh roxygen2 documentation (ALWAYS run after editing function docs)
R -e "roxygen2::roxygenise()"

# Install package locally (without vignettes, faster)
R CMD INSTALL .

# Install with vignettes (slower, for testing documentation)
R -e "devtools::install(build_vignettes = TRUE)"

# Build package tarball
R CMD build .

# Check package (CRAN-style checks, even though not submitting to CRAN)
R CMD check .

# Run test scripts manually
Rscript test_clearsky.R
```

**Note:** There is no automated test suite (`tests/` directory). Use manual test scripts like `test_clearsky.R` for validation.

## Key Conventions

### Function Naming
- Pipeline functions: `pv_[level]_[type]` (e.g., `pv_power_pipeline`, `pv_dc_pipeline`, `pv_clearsky_dc_pipeline`)
- Model functions: `[model]_[phenomenon]` (e.g., `olmo_transposition`, `skoplaki_cell_temperature`, `ineichen_clearsky`)
- Legacy functions: `pv_dc_[models]_pvwatts` (maintained for backward compatibility)
- Utility functions: `[purpose]` or `ensemble_[purpose]` (e.g., `clearsky_index`, `ensemble_summary`)

### Parameter Patterns
- `transposition_model`: `c("haydavies", "reindl", "perez", "olmo")` - haydavies is first/default
- `decomposition_model`: `c("erbs", "boland")` - erbs is first/default
- `cell_temp_model`: `c("skoplaki", "faiman")` - skoplaki is first/default
- `match.arg()` is used internally for argument validation
- First element in `c()` vectors is always the default

### Default Site Parameters (De Aar, South Africa)
- Latitude: -30.6279°
- Longitude: 24.0054°
- Altitude: 1233 m (Mulilo De Aar PV plant)
- Linke turbidity: 3.0 (clean rural conditions)

### Return Value Structure
All pipeline functions return a data frame with:
- **Core columns**: `time`, `GHI`, `G_poa`, `T_air`, `wind`, `T_cell`, `P_dc`, `zenith`, `incidence`
- **Metadata**: `transposition`, `cell_temp`
- **Model-specific**:
  - Olmo: `sun_azimuth`, `k_t`, `I_0`
  - Hay-Davies/Reindl/Perez: `azimuth`, `DNI`, `DHI`, `ai`, `rb`
  - Clear-sky: `ghi_clearsky`, `dni_clearsky`, `dhi_clearsky`, `airmass`
- **IAM**: `iam` column added if `iam_exp` is enabled
- **AC conversion**: `P_ac`, `clipped`, `P_ac_rated` (if using `pv_power_pipeline()` or `pv_clearsky_power_pipeline()`)

## Timezone Handling

The package uses `lubridate` for consistent timezone handling:

1. **Input flexibility**: `time` parameter accepts POSIXct, POSIXlt, character, or numeric timestamps
2. **Automatic UTC conversion**: All solar position calculations use UTC internally
3. **Timezone preservation**: Output timestamps are returned in the original input timezone
4. **Default behavior**: If no timezone is specified, UTC is assumed

**Key functions**: `prepare_time_utc()` and `restore_time_tz()` in `R/time_utils.R`

**Example**:
```r
# Local time input (SAST = UTC+2)
time_local <- as.POSIXct("2026-01-15 12:00", tz = "Africa/Johannesburg")

# Package converts to UTC (10:00 UTC), computes solar position, returns SAST
result <- erbs_decomposition(time = time_local, lat = -30.6, lon = 24.0, GHI = 1000)
attr(result$time, "tzone")  # "Africa/Johannesburg"
```

## Important Model Details

### Transposition Models
- **Hay-Davies** (recommended, default): Requires decomposition first (GHI → DNI/DHI)
- **Reindl**: Similar to Hay-Davies with additional horizon brightening term
- **Perez**: Most sophisticated anisotropic model, uses 8-bin sky classification
- **Olmo**: Direct GHI → POA without decomposition, but only validated for Granada, Spain

**Note:** Hay-Davies, Reindl, and Perez all require decomposition (Erbs or Boland-Ridley) to split GHI into DNI/DHI first. Olmo bypasses this step.

### Decomposition Models
- **Erbs** (default): Piecewise polynomial diffuse fraction model
- **Boland-Ridley**: Logistic function model with configurable coefficients for different time resolutions

### Cell Temperature Models
- **Skoplaki** (default): NOCT-based, two wind model variants
  - Model 1 (default): `h_w = 8.91 + 2.00*v_f`
  - Model 2: `h_w = 5.7 + 3.8*v_w` where `v_w = 0.68*v_f - 0.5`
- **Faiman**: Empirical model from IEC 61853 standard

### Clear-Sky Models
- **Ineichen-Perez**: Estimates clear-sky GHI, DNI, DHI based on:
  - Solar zenith angle (computed from time, lat, lon)
  - Linke turbidity coefficient (atmospheric clarity)
  - Altitude-based atmospheric pressure correction
  - Optional Perez enhancement factor for very clear conditions
- **Haurwitz**: Simple GHI model based only on solar zenith angle
  - Uses formula: `I = 1098 * cos(z) * exp(-0.059 / cos(z))`
  - Estimates DNI and DHI using atmospheric transmission approximations
  - No turbidity or altitude parameters required
  - Suitable for quick estimates when atmospheric data is unavailable
- **Helper functions**: `kasten_young_airmass()`, `atm_pressure_altitude_correction()`, `simple_linke_turbidity()`

### Default Module Parameters (Trina TSM-PC05)
- `P_dc0 = 230` W (nameplate DC power)
- `gamma = -0.0043` /K (temperature coefficient)
- `T_NOCT = 45` °C
- `eta_STC = 0.141` (efficiency)
- `tau_alpha = 0.9` (optical product)

### IAM (Incidence Angle Modifier)
- Power-law model: `cos(θ)^b` where `b = 0.05` (default)
- Set `iam_exp = NA` or `FALSE` to disable
- Corrects for optical losses at high incidence angles

## File Structure Notes

### Core Files
- `R/insol_functions.R` - Solar position calculations (ported from archived `insol` package, GPL-2)
- `R/time_utils.R` - Timezone handling with lubridate (`prepare_time_utc()`, `restore_time_tz()`)

### Model Implementation Files (Ported from pvlib-python)
- `R/erbs_decomposition.R` - Erbs diffuse fraction model
- `R/boland_decomposition.R` - Boland-Ridley logistic model
- `R/haydavies_transposition.R` - Hay-Davies anisotropic sky model
- `R/reindl_transposition.R` - Reindl transposition model
- `R/perez_transposition.R` - Perez 8-bin transposition model
- `R/olmo_transposition.R` - Olmo clearness index method (not from pvlib)
- `R/skoplaki_cell_temperature.R` - Skoplaki NOCT-based model
- `R/faiman_cell_temperature.R` - Faiman IEC 61853 model
- `R/pvwatts_dc.R` - PVWatts DC power model
- `R/ineichen_clearsky.R` - Ineichen-Perez clear-sky model
- `R/haurwitz_clearsky.R` - Haurwitz clear-sky model

### Pipeline Files
- `R/pv_dc_pipeline.R` - Main modular DC pipeline with model selection
- `R/pv_power_pipeline.R` - DC + AC pipeline
- `R/pv_clearsky_pipeline.R` - Clear-sky pipelines and performance metrics
- `R/pv_ac_simple_clipping.R` - Simple inverter model

### Ensemble Files
- `R/pv_dc_ensemble.R` - Run all 8 DC model combinations
- `R/pv_power_ensemble.R` - Run all 8 AC model combinations
- `R/ensemble_utils.R` - Utility functions for ensemble analysis

### Legacy Files (Backward Compatibility)
- `R/pv_dc_olmo_skoplaki_pvwatts.R` - Legacy convenience function
- `R/pv_dc_haydavies_faiman_pvwatts.R` - Legacy convenience function

### Documentation
- `vignettes/de_aar.Rmd` - Main vignette with De Aar solar plant example
- `test_clearsky.R` - Manual test script for clear-sky implementation

## Code Attribution and Licensing

**Critical:** Much of this package is ported from [pvlib-python](https://github.com/pvlib/pvlib-python) (BSD 3-Clause License). Key components include:
- Ineichen-Perez clear-sky model
- Perez, Reindl, Hay-Davies transposition models
- Erbs and Boland-Ridley decomposition models
- PVWatts DC power model
- Kasten-Young airmass formula

Solar position calculations are from the archived `insol` package by Javier G. Corripio (GPL-2).

When implementing new models, check pvlib-python first for validated implementations.

## Known Limitations

1. **Olmo model**: Significant errors (RMSE 21-52%) outside Granada, Spain - recommend Hay-Davies for production
2. **No automated tests**: Package has no `tests/` directory. Use manual test scripts for validation
3. **GitHub-only**: Not intended for CRAN submission
4. **No Linke turbidity database**: Use `simple_linke_turbidity()` for rough estimates or provide site-specific values
