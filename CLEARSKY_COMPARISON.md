# Clear-Sky Model Comparison: Ineichen-Perez vs Haurwitz

## Overview

This document compares the two clear-sky models implemented in pvflux:

1. **Ineichen-Perez Model**: Complex atmospheric model requiring turbidity data
2. **Haurwitz Model**: Simple geometric model requiring only solar position

## Key Findings

### Accuracy and Sensitivity

| Aspect | Ineichen-Perez | Haurwitz |
|--------|----------------|----------|
| **GHI (Peak)** | 1043-1170 W/m² | 1011 W/m² |
| **Daily Energy** | 26.0-29.4 kWh | 26.0 kWh |
| **Turbidity Sensitivity** | High (126 W/m² variation) | None |
| **RMSE vs Ineichen** | N/A | 62.8 W/m² |
| **Energy Difference** | N/A | -0.1 to -3.4 kWh (-0.1% to -11.6%) |

### Model Characteristics

#### Ineichen-Perez Model
- **Pros**:
  - Most accurate for specific locations with known atmospheric conditions
  - Accounts for altitude effects
  - Optional Perez enhancement for very clear conditions
  - Well-validated against measurements

- **Cons**:
  - Requires Linke turbidity input
  - More complex (requires turbidity database for auto-lookup)
  - Computationally more intensive
  - Significant sensitivity to atmospheric conditions

#### Haurwitz Model
- **Pros**:
  - Simple implementation (only needs solar position)
  - No atmospheric parameters required
  - Faster computation
  - Good for initial estimates and quick calculations
  - Consistent results regardless of location/time

- **Cons**:
  - Less accurate (does not account for atmospheric conditions)
  - Fixed formula regardless of location/time of year
  - May overestimate/underestimate in specific conditions
  - No altitude correction

## Performance Comparison

### Under Different Atmospheric Conditions

| Linke Turbidity | Ineichen GHI (W/m²) | Haurwitz GHI (W/m²) | Difference (W/m²) |
|-----------------|-------------------|-------------------|------------------|
| Very Clean (2.0) | 1170 | 1011 | +159 |
| Clean (3.0) | 1137 | 1011 | +126 |
| Urban (4.0) | 1105 | 1011 | +94 |
| Polluted (5.0) | 1073 | 1011 | +63 |
| Very Polluted (6.0) | 1043 | 1011 | +32 |

### Daily Energy Production

| Atmospheric Condition | Ineichen Energy (kWh) | Haurwitz Energy (kWh) | Difference |
|----------------------|----------------------|---------------------|-----------|
| Very Clean | 29.4 | 26.0 | -3.4 (-11.6%) |
| Clean | 28.5 | 26.0 | -2.5 (-8.9%) |
| Urban | 27.7 | 26.0 | -1.7 (-6.0%) |
| Polluted | 26.9 | 26.0 | -0.8 (-3.1%) |
| Very Polluted | 26.0 | 26.0 | -0.0 (-0.1%) |

## When to Use Each Model

### Use Ineichen-Perez When:
- You have accurate Linke turbidity data for the location
- Performing detailed resource assessment or feasibility studies
- Historical performance analysis with known atmospheric conditions
- High accuracy is critical for financial modeling
- Working at high altitudes (>1000m)

### Use Haurwitz When:
- Quick estimates are needed without atmospheric data
- Preliminary site screening or assessment
- Educational purposes or demonstrations
- When turbidity data is unavailable
- For regions with consistently clear atmospheric conditions

## Recommendations

1. **For Professional Studies**: Use Ineichen-Perez with appropriately determined turbidity values
2. **For Preliminary Analysis**: Use Haurwitz for quick estimates
3. **For Performance Monitoring**: Use Ineichen-Perez with measured turbidity or database lookup
4. **For Educational Purposes**: Both models provide good insights, Haurwitz is simpler to explain

## Implementation Notes

Both models return identical data structures, making it easy to switch between them:

```r
# Ineichen-Perez (default, more accurate)
result <- pv_clearsky_power_pipeline(
  time, lat, lon, T_air, wind, tilt, azimuth,
  clearsky_model = "ineichen",
  linke_turbidity = 3.0
)

# Haurwitz (simpler, no turbidity needed)
result <- pv_clearsky_power_pipeline(
  time, lat, lon, T_air, wind, tilt, azimuth,
  clearsky_model = "haurwitz"
)
```

## Future Improvements

The Haurwitz model could be enhanced by:
1. Adding altitude corrections based on atmospheric pressure
2. Incorporating simple atmospheric clarity estimates
3. Providing seasonal adjustments based on typical atmospheric conditions

However, the simplicity of the current implementation is its strength for use cases where detailed atmospheric data is unavailable.