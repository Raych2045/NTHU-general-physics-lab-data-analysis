# Physics Lab Analysis

Python scripts for data analysis in undergraduate physics experiments.

## Requirements
- Python 3.12
- NumPy
- Pandas
- Matplotlib
- SciPy
- openpyxl (needed by Pandas to read `.xlsx` files)

## Workflow

1. Fill in each lab's `data.xlsx` (and `data2.xlsx`, etc., where a lab has more than one) with measurements. Column meanings are documented per lab below.
2. Run that lab's `lab.py`. It can be run from anywhere — `BASE_DIR` is derived from the script's own location (`Path(__file__).resolve().parent.parent`), so the current working directory doesn't matter.
3. `lab.py` reads the spreadsheet(s), recomputes every derived quantity, then:
   - clears and rewrites a single `report.tex` per lab with all LaTeX table fragments (so reruns are idempotent, not appended);
   - regenerates every figure into `figures/`, overwriting the old ones.
4. Copy the relevant tables from `report.tex` and figures from `figures/` into the matching `Overleaf Projects/post-lab <N>/main.tex`.
5. Shared helpers (`pm`, `percent`, `df_to_latex`, `typeA_uncertainty`, `linear_model`, `compute_r_squared`) live in the root-level `package.py` and are imported by every `lab.py` via `sys.path.insert(0, str(Path(__file__).resolve().parent.parent))`.

## data.xlsx reference

### Lab 13 — Galvanometer and self-made meters

**Columns**

| Column | Unit | Used in | Meaning |
|---|---|---|---|
| `subexp_id` | – | all | Sub-experiment symbol, A–E |
| `V_multimeter` | V | A, B, D, E | Multimeter voltage reading (ground truth) |
| `I_Galvanometer` | µA | all | Raw galvanometer deflection current |
| `I_multimeter` | mA | C | Multimeter current reading (ground truth) |
| `R_load` | Ω | all | External/load resistor value |
| `uncertainty_R_load(%)` | % | A | Tolerance of `R_load` |
| `R_used` | Ω | C, D, E | Resistor used to build the self-made meter |
| `uncertainty_R_used(%)` | % | C, D, E | Tolerance of `R_used` |
| `R_load_multimeter` | Ω | E | Multimeter-measured load resistance (ground truth) |
| `V_power_supply` | V | C, D | Power-supply voltage setting |

**Sub-experiments**

| Symbol | Meaning |
|---|---|
| A | Calibrate the galvanometer: relate raw deflection current to true current via a known `R_load` |
| B | Measure the galvanometer's internal resistance $R_G$ |
| C | Self-made ammeter — compared against a multimeter |
| D | Self-made voltmeter — compared against a multimeter |
| E | Self-made ohmmeter — compared against a multimeter |

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 22 | `u_IG` | $2/\sqrt{12}$ | µA | all | Type-B uncertainty of the galvanometer current reading (2 µA reading resolution) |
| 146 | literal `[1.65, √227, √226]` (assigned to `u_R_used`) | $1.65, \sqrt{15^2+1+1}, \sqrt{15^2+1}$ | Ω | C | Per-row Type-B uncertainty of `R_used`, one value per row of sub-exp C. Order-dependent on the row order in `data.xlsx` — check it after modifying/adding/reordering rows |
| 205, 209 | literal `0.00005` | $5\times10^{-5}$ | A | E | Galvanometer full-scale deflection current, used in the self-made-ohmmeter formula |

### Lab 14 — Lorentz force and current balance

**Columns**

| Column | Unit | Used in | Meaning |
|---|---|---|---|
| `subexp` | – | all | Sub-experiment symbol, A–D |
| `mass` | g | all | Mass on the balance pan (force $F = mg$) |
| `board_number` | – | A | Circuit-board identifier |
| `outer_length` | mm | A | Outer-edge length of the current loop |
| `inner_length` | mm | A | Inner-edge length of the current loop |
| `magnet_number` | – | B | Number of magnets used to vary $B$ |
| `B_values` | T | B | Magnetic field strength |
| `I_power_supply` | A | C, D | Power-supply current |

**Sub-experiments**

| Symbol | Meaning |
|---|---|
| A | Lorentz force vs. wire length — find $B$ from the slope |
| B | Lorentz force vs. number of magnets — find supply current from the slope |
| C | Current balance: force vs. $I^2$ — find the permeability of free space $\mu_0$ |
| D | Predict the mass of an unknown object from the $F$–$I^2$ fit in C |

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 19 | `g` | 9.80665 | m/s² | all | Standard gravity, converts mass to force ($F=mg$) |
| 20 | `I_const_A` | 3.04 | A | A | Power-supply current reading used to compute $B_{\text{obs}}=s/I$ |
| 21 | `B_ref` | $0.09677-0.00064=0.09613$ | T | A | Reference magnetic field value (probe reading corrected for offset, if you use geomagnetism as the zero point and didn't calibrate gaussmeter by inserting the probe into the zero magnetic field tube, you can ignore the offset term.) |
| 23 | `D` | 0.0177 | m | C | Light spot displacement on scale |
| 24 | `a_param` | 0.2115 | m | C | Mirror-to-moving-wire distance |
| 25 | `b_param` | 1.0065 | m | C | Mirror-to-scale distance |
| 26 | `R` | 0.0031 | m | C | Moving-wire diameter $r_0$ |
| 27 | `L` | 0.279 | m | C | Moving-wire length |
| 28 | `mu0` | $4\pi\times10^{-7}$ | N/A² | C | Permeability of free space — theoretical reference for $\mu_{\text{est}}$ |
| 30 | `I_ref_B` | 3.04 | A | A | Power-supply current reading as reference value |
| 128, 129 | literal `0.038` | 0.038 | m | B | Inner length of circuit board #38 — a magic number; update it if you switch the board or you have your own measured value. Remember to check line 233 and 254 also. |

### Lab 17 — I-V curves and Planck's constant

**`data.xlsx` columns**

| Column | Unit | Used in | Meaning |
|---|---|---|---|
| `subexp_id` | – | all | Sub-experiment symbol, A–C |
| `V` | V (A, C); mV (B) | all | Applied voltage |
| `I` | mA | all | Measured current |
| `color` | text | C | LED color — must be `"red"`, `"yellow"`, `"green"`, or `"blue"` (lowercase only) |

**`turnonvoltage.xlsx` columns** (cross-class turn-on voltages, all LED colors)

| Column | Unit | Meaning |
|---|---|---|
| `lambda` | nm | LED wavelength |
| `VF` | V | Turn-on (forward) voltage |

**Sub-experiments**

| Symbol | Meaning |
|---|---|
| A | I–V curve of a resistor |
| B | I–V curve of a silicon diode (reverse-bias included) |
| C | I–V curve of an LED — the turn-on voltage is used to estimate Planck's constant |

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 20 | `R_used` | 1974 | Ω | A | Reference resistor value, compared against $R_{\text{obs}}=1/s$ |
| 22 | `lambda_red` | 635 | nm | C | Reference wavelength of the red LED |
| 23 | `lambda_yel` | 590 | nm | C | Reference wavelength of the yellow LED |
| 24 | `lambda_gre` | 525 | nm | C | Reference wavelength of the green LED |
| 25 | `lambda_blu` | 460 | nm | C | Reference wavelength of the blue LED |
| 27 | `e_charge` | $1.602176634\times10^{-19}$ | C (coulomb) | C | Elementary charge |
| 28 | `c_light` | $2.99792458\times10^{8}$ | m/s | C | Speed of light |
| 29 | `h_true` | $6.62607015\times10^{-34}$ | J·s | C | Planck's constant (reference value for comparison) |

### Lab 18 — Polarization and diffraction

`data.xlsx` has one sheet per sub-experiment.

**Sheet `A`** — acrylic plate: refractive index from circular area diameter

| Column | Unit | Meaning |
|---|---|---|
| `description` | text | Label shown as the table's column header |
| `d_mm` | mm, `;`-separated | Measured plate thickness |
| `2R_mm` | mm | Circular area diameter |
| `R2_note` | text | `"暗圈"` or `"亮圈"` |

**Sheet `B`** — acrylic rod: refractive index from rangefinder distances

| Column | Unit | Meaning |
|---|---|---|
| `D_acry_m` | m | Rangefinder distance through the acrylic rod |
| `D_air_m` | m | Rangefinder distance through air |
| `l_mm` | mm | Rangefinder length |

**Sheet `C`** — Brewster's angle: refractive index from the polarizing angle

| Column | Unit | Meaning |
|---|---|---|
| `theta_B_deg` | degree | Observed Brewster's angle |

**Sheets `malus_two` / `malus_three`** — Malus' law with two/three polarizers

| Column | Unit | Meaning |
|---|---|---|
| `theta_deg` | degree | Angle between polarizer axes |
| `Iobs` | lux, `;`-separated | Repeated raw intensity readings for that angle, e.g. `"12.3;12.5;12.1"` |

**Sheet `single_slit`** — single-slit diffraction: wavelength from fringe width

| Column | Unit | Meaning |
|---|---|---|
| `a_mm` | mm | Slit width |
| `W_mm` | mm | Central fringe width |

**Sheet `double_slit`** — double-slit diffraction: slit width/spacing

| Column | Unit | Meaning |
|---|---|---|
| `label` | text | Slit identifier |
| `W_mm` | mm | Central fringe width |
| `span_mm` | mm | Total span measured across `N_fringes` fringes |
| `N_fringes` | count | Number of fringes spanned by `span_mm` |
| `a_ref_mm` | mm | Manufacturer-stated slit width (reference value) |
| `d_ref_mm` | mm | Manufacturer-stated slit spacing (reference value) |

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 19 | `n_ref` | 1.49 | – | A, B, C | Reference refractive index of acrylic |
| 20 | `lambda_ref_nm` | 532 | nm | single_slit, double_slit | Reference laser wavelength |
| 61 | literal `0.02/√12` | 0.02 mm resolution | mm | A | Type-B uncertainty of the thickness $d$ measurement (vernier caliper resolution) |
| 109 | literal `0.001/√12` | 1 mm resolution | m | B | Type-B uncertainty of $D_{\text{acry}}$ (rangefinder resolution) |
| 110 | literal `0.001/√12` | 1 mm resolution | m | B | Type-B uncertainty of $D_{\text{air}}$ (rangefinder resolution) |
| 111 | literal `0.02/√12` | 0.02 mm resolution | mm | B | Type-B uncertainty of $l$ (vernier caliper resolution) |
| 268 | literal `1.0/√12` | 1° resolution | degree | C | Type-B uncertainty of $\theta_B$ (rotation-stage resolution) |
| 302 | `L_ss` | 77.85 | cm | single_slit | Slit-to-screen distance |
| 336 | `L_ds` | 86.5 | cm | double_slit | Slit-to-screen distance |

### Lab 19 — RC and RLC circuits

**`data1.xlsx` columns** (read by the Gain/frequency-response part)

| Column | Unit | Meaning |
|---|---|---|
| `subexp_id` | – | Sub-experiment symbol, A–C |
| `T` | × $\tau$ | Square-wave period, as a multiple of the RC time constant $\tau$ |
| `Vin` | V (Vpp) | Input peak-to-peak voltage |
| `Vout` | V (Vpp) | Output peak-to-peak voltage |

**Sub-experiments**

| Symbol | Meaning |
|---|---|
| A | RC low-pass filter — sinusoidal input, measure capacitor voltage |
| B | RC high-pass filter — sinusoidal input, measure resistor voltage |
| C | Square-wave input, measure capacitor voltage |

**`data2.xlsx` columns** (RLC resonance — single experiment, no sub-experiment split)

| Column | Unit | Meaning |
|---|---|---|
| `f` | Hz | Driving frequency |
| `Vout` | V (Vpp) | Output peak-to-peak amplitude |
| `phase` | degree | Phase of the output relative to the input |

Note: the damped-oscillation fit (peak time vs. voltage to find $\beta$ and $\omega_1$) uses hardcoded arrays in line 57 and 58 in `lab.py` rather than a spreadsheet, since it's a single short dataset.

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 21 | `tau` | $1.005\times10^{-3}$ | s | A, B, C | RC time constant $\tau=R_{\text{sum}}C$ |
| 22 | `C_cap` | $1\times10^{-9}$ | F | RLC damped and RLC resonance | RLC Circuit capacitance $C$ |
| 25 | `R_damp` | $100+50+48.1=198.1$ | Ω | damped | Total series resistance for the damped-oscillation (RLC) circuit (including signal generator's and inductor's own resistances) |
| 26 | `L_damp` | $22\times10^{-3}$ | H | damped | Inductance for the damped-oscillation circuit |
| 31 | `R_reson` | $100+50+25.8=175.8$ | Ω | resonance | Total series resistance for the resonance (RLC) circuit |
| 32 | `L_reson` | $10\times10^{-3}$ | H | resonance | Inductance for the resonance circuit |
| 33 | `Vin_reson` | 1.02 | V (Vpp) | resonance | Input amplitude for the resonance experiment |

### Lab 21 — Helmholtz coil

**Columns**

| Column | Unit | Used in | Meaning |
|---|---|---|---|
| `subexp_id` | – | all | Sub-experiment symbol, A–E |
| `d` | cm | all | Pick-up coil position along the axis |
| `Vpp` | mV | all | Peak-to-peak voltage measured on CH1 |

**Sub-experiments**

| Symbol | Meaning |
|---|---|
| A | Single coil — magnetic field distribution along the axis |
| B | Coil pair, currents flow in same direction, separation = 5 cm |
| C | Helmholtz coil pair, separation = 10.5 cm (= coil radius) |
| D | Coil pair, currents flow in same direction, separation = 15 cm |
| E | Anti-Helmholtz coil pair, separation = 10.5 cm |

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 17 | `R` | 0.105 | m | all | Field-coil radius, used in the on-axis single-coil magnetic field formula |
| 28, 31, 34, 37 | literal `d` (0.05, 0.105, 0.15, 0.105) | m | B, C, D, E | Coil-pair separation for each sub-experiment's simulated field |

### Lab 22 — Michelson interferometer

**Columns**

| Column | Unit | Used in | Meaning |
|---|---|---|---|
| `subexpid` | – | all | Sub-experiment symbol, A–C |
| `initial` | mm (A); degree (B); kPa (C) | all | Micrometer/angle/pressure reading before the count |
| `final` | mm (A); degree (B); kPa (C) | all | Micrometer/angle/pressure reading after the count |
| `deltaN` | count | all | Number of fringes counted between `initial` and `final` |

**Sub-experiments**

| Symbol | Meaning |
|---|---|
| A | Measure laser wavelength from micrometer travel vs. fringe count |
| B | Measure glass refractive index from rotation angle vs. fringe count |
| C | Measure air refractive index from pressure change vs. fringe count, extrapolated to 1 atm |

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 14 | `lambda_ref` | $657\times10^{-9}$ | m | all | Reference laser wavelength |
| 15 | `t` | $6\times10^{-3}$ | m | B | Glass plate thickness |
| 16 | `n_glass_ref` | 1.46 | – | B | Reference refractive index of glass |
| 17 | `h` | $3\times10^{-2}$ | m | C | Air-chamber thickness |
| 18 | `p_atm` | 101.325 | kPa | C | Atmospheric pressure, used to extrapolate $n_{\text{air}}$ from the fit |
| 93 | `theta_resolution_deg` | 1 | degree | B | Resolution of the rotation-stage angle reading, used to calculate Type-B uncertainty |

### Lab 23 — Spectrum of atomic hydrogen and Planck constant

**Columns**

| Column | Unit | Used in | Meaning |
|---|---|---|---|
| `subexp_id` | – | all | Sub-experiment symbol, A–B |
| `m` | – (order) | A | Diffraction order |
| `ni` | – | all | Initial-state principal quantum number of the Balmer line |
| `theta_main` | degree | all | Main-scale reading at the spectral line |
| `theta_aux` | arcminute | all | Vernier reading at the spectral line |
| `theta_main_ref` | degree | all | Main-scale reading at the zero-order (reference) position |
| `theta_aux_ref` | arcminute | all | Vernier reading at the zero-order (reference) position |

**Sub-experiments**

| Symbol | Meaning |
|---|---|
| A | Diffraction grating — measure Balmer-series wavelengths, fit Rydberg constant and Planck's constant |
| B | Prism — measure minimum-deviation angle per line to find the prism's refractive indices |

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 14 | `d` | $10000/6$ | nm | A | Grating spacing (600 lines/mm) |
| 15–18 | `lambda_ref_3`…`lambda_ref_6` | 656.3, 486.1, 434, 410.2 | nm | all | Reference Balmer-series wavelengths for $n_i=3,4,5,6$ |
| 19 | `alpha` | $\pi/3$ | rad | B | Prism apex angle |
| 21 | `h_ref` | $6.62607015\times10^{-34}$ | J·s | A | Reference Planck's constant |
| 22 | `c` | $2.99792458\times10^{8}$ | m/s | A | Speed of light |
| 23 | `e` | $1.602176634\times10^{-19}$ | C | A | Elementary charge |
| 24 | `m_e` | $9.1093837\times10^{-31}$ | kg | A | Electron mass |
| 25 | `eps0` | $8.854187817\times10^{-12}$ | F/m | A | Vacuum permittivity |
| 26 | `m_p` | $1.67262192\times10^{-27}$ | kg | A | Proton mass, used to compute the reduced mass $\mu$ |
