# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

pvflux is an R package for solar PV power forecasting with a modular pipeline architecture. Users can independently select transposition models (Hay-Davies, Olmo) and cell temperature models (Skoplaki, Faiman) and run ensemble analysis across all combinations.

**Important:** The Olmo transposition model has known validation issues (RMSE of 21-52%) outside its calibration region (Granada, Spain). Default is now "haydavies" - see `?olmo_transposition` for details.

## Architecture

The package uses a **three-level abstraction**:

1. **Individual model functions** - Atomic components for each physical model
2. **Modular pipeline** - `pv_dc_pipeline()` allows independent model selection
3. **Convenience functions** - `pv_power_pipeline()` (DC+AC), `pv_power_ensemble()` (all combinations)

### Data Flow Pipeline

```
Input (time, lat, lon, GHI, T_air, wind, tilt, azimuth)
    ↓
Step 1: Transposition (GHI → G_poa)
    ├─ olmo_transposition()
    └─ haydavies_transposition() (requires erbs_decomposition first)
    ↓
Step 2: Cell Temperature (G_poa, T_air, wind → T_cell)
    ├─ skoplaki_cell_temperature()
    └─ faiman_cell_temperature()
    ↓
Step 3: DC Power (G_poa, T_cell → P_dc)
    └─ pvwatts_dc()
    ↓
Step 4: AC Conversion (P_dc → P_ac)
    └─ pv_ac_simple_clipping()
```

### Key Model Combinations

| Identifier | Transposition | Cell Temperature |
|------------|---------------|------------------|
| haydavies_skoplaki | Hay-Davies | Skoplaki (default) |
| haydavies_faiman | Hay-Davies | Faiman |
| olmo_skoplaki | Olmo | Skoplaki |
| olmo_faiman | Olmo | Faiman |

## Common Development Commands

```bash
# Generate/refresh roxygen2 documentation
R -e "roxygen2::roxygenise()"

# Build package
R CMD build .

# Check package (runs tests)
R CMD check .

# Install locally with vignettes
R -e "devtools::install(build_vignettes = TRUE)"

# Install from current directory
R CMD INSTALL .
```

## Key Conventions

### Function Naming
- Pipeline functions: `pv_[level]_[type]` (e.g., `pv_power_pipeline`, `pv_dc_pipeline`)
- Model functions: `[model]_[phenomenon]` (e.g., `olmo_transposition`, `skoplaki_cell_temperature`)
- Legacy functions: `pv_dc_[models]_pvwatts` (maintained for backward compatibility)

### Parameter Patterns
- `transposition_model`: `c("haydavies", "olmo")` - haydavies is first/default
- `cell_temp_model`: `c("skoplaki", "faiman")` - skoplaki is first/default
- `match.arg()` is used internally for argument validation

### Return Value Structure
All pipeline functions return a data frame with:
- **Core columns**: `time`, `GHI`, `G_poa`, `T_air`, `wind`, `T_cell`, `P_dc`, `zenith`, `incidence`
- **Metadata**: `transposition`, `cell_temp`
- **Model-specific**: `sun_azimuth` (olmo), `azimuth/DNI/DHI/ai/rb` (haydavies)
- **IAM**: `iam` column added if `iam_exp` is enabled

## Important Model Details

### Transposition Models
- **Hay-Davies** (recommended): Requires `erbs_decomposition()` first (GHI → DNI/DHI)
- **Olmo**: Direct GHI → POA, but only validated for Granada, Spain

### Cell Temperature Models
- **Skoplaki**: NOCT-based, two wind model variants
  - Model 1: `h_w = 8.91 + 2.00*v_f`
  - Model 2: `h_w = 5.7 + 3.8*v_w` where `v_w = 0.68*v_f - 0.5`
- **Faiman**: Empirical model from IEC 61853 standard

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

- `R/insol_functions.R` - Solar position calculations (ported from the archived `insol` package, used under GPL-2)
- `vignettes/de_aar.Rmd` - Main vignette with De Aar solar plant example
- Legacy functions (`pv_dc_olmo_skoplaki_pvwatts`, `pv_dc_haydavies_faiman_pvwatts`) maintained for backward compatibility

## Known Limitations

1. **Olmo model**: Significant errors outside southern Spain - recommend Hay-Davies for production use
2. **No tests**: Package currently has no automated test suite
3. **GitHub-only**: Not intended for CRAN submission
