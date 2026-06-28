# Physics Lab Analysis

Python scripts for data analysis in undergraduate physics experiments.

## Requirements for running `lab.py`
- Python 3.10 or later versions
- NumPy
- Pandas
- Matplotlib
- SciPy
- openpyxl (needed by Pandas to read `.xlsx` files)

## Workflow

1. Fill in each lab's `data.xlsx` (and `data2.xlsx`, etc., if a lab has more than one) with measurements. Column meanings are documented per lab below.
2. Run that lab's `lab.py`. It can be run from anywhere. The current working directory doesn't matter.
3. `lab.py` reads the spreadsheet(s), recomputes every derived quantity, then:
   - clears and rewrites a single `report.tex` per lab with all LaTeX table fragments (so reruns are idempotent, not appended);
   - regenerates every figure into `figures/`, overwriting the old ones.

### Details

1. `BASE_DIR` is derived from the script's own location (`Path(__file__).resolve().parent.parent`)
2. Shared helpers (`pm`, `percent`, `df_to_latex`, `typeA_uncertainty`, `linear_model`, `compute_r_squared`) live in the root-level `package.py` and are imported by every `lab.py` via `sys.path.insert(0, str(Path(__file__).resolve().parent.parent))`.

## How to use `report.tex`?

It only contains LaTeX tables, and direct compilation will not produce any results.

The recommended approach is to import the images from `figures/` and `report.tex` together into `latexplate.tex`, and then compile it locally or on Overleaf. You should then get a PDF file.

Of course, you can also use your preferred preamble to format the document, but you must include the following three LaTeX packages: `booktabs`, `siunitx`, and `xeCJK`.

The compiled results of all data table can be found in `example/`, with one PDF file per lab.

## data.xlsx reference

### Lab 1 — Measurement of density

`data.xlsx` has three sheets, one per test object, each storing 20 raw trials. You may add new rows if needed.

**Sheet `pvc`** — PVC pipe

| Column | Unit | Meaning |
|---|---|---|
| `inner_diameter_cm` | cm | Inner diameter (vernier caliper) |
| `outer_diameter_cm` | cm | Outer diameter (vernier caliper) |
| `length_cm` | cm | Tube length (vernier caliper) |
| `mass_g` | g | Mass (electronic balance) |
| `buoyancy_gw` | gw | Balance reading while fully submerged — by Archimedes' principle (water density is assumed to be 1 g/cm³) this numerically equals the tube's volume in cm³ |

**Sheet `cylinder`** — solid cylinder

| Column | Unit | Meaning |
|---|---|---|
| `diameter_cm` | cm | Base diameter (vernier caliper) |
| `height_cm` | cm | Height (vernier caliper) |
| `mass_g` | g | Mass (electronic balance) |

**Sheet `steel_ball`** — steel ball

| Column | Unit | Meaning |
|---|---|---|
| `zero_reading_mm` | mm | Micrometer zero-point offset |
| `diameter_raw_mm` | mm | Raw micrometer diameter reading, before zero-point correction |
| `mass_g` | g | Mass (electronic balance) |

`lab.py` computes the corrected diameter (`diameter_raw_mm - zero_reading_mm`) dynamically rather than using stored values.

**Parts**

| Part | Meaning |
|---|---|
| PVC pipe | Inner/outer radius → circle areas → base area (outer − inner) → volume → density = mass/volume; cross-checked against a second density from Archimedes' principle (mass / buoyancy reading) |
| Cylinder | Radius → base area → volume → density = mass/volume |
| Steel ball | Radius (from the corrected diameter) → volume ($4\pi r^3/3$) → density = mass/volume |

Tables 1–3 (per-object raw-measurement statistics) report the mean to 5 decimal places and all other columns (std. dev., Type A/B/combined uncertainty) to 3 significant figures; final rounding is left to users.

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 9 | `VERNIER_RES_CM` | 0.005 | cm | all (lengths) | Vernier caliper least count (0.05 mm), Type-B uncertainty source |
| 10 | `MICROMETER_RES_MM` | 0.01 | mm | Steel ball | Micrometer caliper least count, Type-B uncertainty source |
| 11 | `BALANCE_RES_G` | 0.01 | g | all (mass/buoyancy) | Electronic balance least count, Type-B uncertainty source |

### Lab 2 — Newton's second law and forces

`data.xlsx` has two sheets.

**Sheet `incline`** — raw acceleration of the cart on the incline, at each hanging mass

| Column | Unit | Meaning |
|---|---|---|
| `mass_g` | g | Hanging mass |
| `acceleration_cm_s2` | cm/s², `,`-separated | The 3 raw trial readings at that mass, e.g. `"35.60,35.46,35.24"` |
| `group` | – | `a_minus` (mass range where the cart's true signed acceleration is negative, before the zero crossing) or `a_plus` (after) — explicit column rather than inferred from row position, so reordering/editing rows won't silently break the split |

**Sheet `horizontal`** — raw acceleration of the cart on a leveled (horizontal) track, no hanging mass, used as an independent friction measurement

| Column | Unit | Meaning |
|---|---|---|
| `left_cm_s2` | cm/s² | Acceleration, cart sliding left |
| `right_cm_s2` | cm/s² | Acceleration, cart sliding right |

**Parts**

| Part | Meaning |
|---|---|
| $a_{-}(m)$ / $a_{+}(m)$ | At each hanging mass $m$, the cart's incline acceleration is recorded as an unsigned magnitude; on the unified $a$-$m$ graph the $a_{-}$ group's values are negated (sign-flipped) to recover the true signed acceleration. Both groups are fit linearly ($a_{+}$ restricted to its `N_FIT_A_PLUS_LINEAR` smaller-$m$ points — the larger masses break the linear approximation) to get each line's slope and x-intercept ($m_{-}$ and $m_{+}$) |
| Slope-method results | $M=g/s_{-}$, $\|\mathbf{F_f}\|=(m_{+}-m_{-})g/2$, $\theta=\arcsin\dfrac{(m_{+}+m_{-})s_{-}}{2g}$ ($s_-$ = the $a_-(m)$ slope), compared against an independent method: $M$ from a direct balance reading, $\|\mathbf{F_f}\|$ from $M_{\text{ref}}\times(\text{mean horizontal-track deceleration})$, $\theta$ from $\arcsin(\text{height difference}/\text{incline length})$ |
| Nonlinear curve fit | Both groups (now using **all** of $a_+$'s points, since the nonlinear model — unlike the linear fit — already captures the $m+M$ denominator's curvature) are fit to $\dfrac{gm-g\sin\theta\, M\pm F}{m+M}$ ($\theta$ fixed at the independent trig-method value; $+F$ for $a_-$, $-F$ for $a_+$), then $M$ and $F$ from the two fits are averaged and compared in the same way |

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 11 | `G_CGS` | 980.665 | cm/s² | all | Standard gravity — this lab works throughout in cm-g-s (CGS) units, not SI |
| 12 | `M_CART_REF_G` | 382.70 | g | comparison tables | Cart mass from a direct balance reading (the "other method" reference) |
| 13 | `SLOPE_RULER_LENGTH_CM` | 113.92 | cm | comparison tables, curve fit | Ruler-measured incline length |
| 14 | `SLOPE_HEIGHT_DIFF_CM` | 4.55 | cm | comparison tables, curve fit | Ruler-measured height difference between the incline's two ends — together with the line above, gives the independent trig-method $\theta$ used both as a comparison target and as the fixed $\theta$ inside the curve-fit model |
| 16 | `N_FIT_A_PLUS_LINEAR` | 15 | count | $a_+(m)$ linear fit only | Number of points corresponding to a relatively small $m$ kept for the **linear** fit (the curve fit uses all) |
| 17 | `ACCEL_RES_CM_S2` | 0.01 | cm/s² | all (acceleration) | Photogates' acceleration-mode reading resolution, Type-B uncertainty source (combined with Type A for every acceleration mean) |

Notes:
- The x-intercept uncertainty in the linear-fit tables includes the slope/intercept covariance term (`pcov[0,1]` from `curve_fit`).
- All fits read directly from the raw data in `data.xlsx`, not from already-rounded table values.

### Lab 3 — Centripetal force

`data.xlsx` has one sheet per experiment, each storing raw (ungrouped) period readings.

**Sheet `Exp1`**

| Column | Unit | Meaning |
|---|---|---|
| `m_g` | g | Hanging mass as a equivalent of centripetal force, varied across rows; multiple raw rows per mass |
| `T_s` | s | Raw measured period of rotation |

**Sheet `Exp2`**

| Column | Unit | Meaning |
|---|---|---|
| `r_cm` | cm | Rotation radius, varied across rows; multiple raw rows per radius |
| `T_s` | s | Raw measured period of rotation |

**Experiments**

| Symbol | Meaning |
|---|---|
| Exp1 | $M$ and $r$ fixed, hanging mass $m$ varied — infer $m$ from each measured period ($mg=4\pi^2Mr/T^2$), then fit $T^{-2}$ vs. $m$ to estimate $M$ from the slope |
| Exp2 | $M$ and $m$ fixed, rotation radius $r$ varied — fit $T^2$ vs. $r$ to estimate $m$ from the slope; cross-checked by substituting Exp2's period at the radius matching Exp1's fixed $r$ into Exp1's $T^{-2}$–$m$ fit line |

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 14 | `g` | 9.8 | m/s² | both | Standard gravity |
| 15 | `M_ref` | 209.2e-3 | kg | both | Rotating-body mass (balance reading), fixed in both experiments |
| 16 | `r1` | 0.1200 | m | both | Rotation radius fixed in Experiment 1; also the radius in Experiment 2 used for the cross-check |
| 17 | `m_ref` | 55.34e-3 | kg | Exp2 | Hanging-mass balance reading, fixed in Experiment 2 |
| 23 | `T_TYPEB_S` | $0.010/\sqrt{12}$ | s | both | Type-B uncertainty of $T$ from the 10 ms stopwatch resolution, combined with the Type-A (std/√N) uncertainty per mass/radius group |.

### Lab 5 — Moments of inertia

`data.xlsx` has three sheets.

**Sheet `part1_3_timestamps`**

Each column is one f-t regression's raw stopwatch-reading chain (ascending, ms). Consecutive entries form a $(T_1,T_2)$ pair (`T1=chain[:-1]`, `T2=chain[1:]`).

| Column | Meaning |
|---|---|
| `hori_loaded`, `hori_friction` | Part 1, horizontal disk — loaded and friction (unloaded) runs |
| `verti_loaded`, `verti_friction` | Part 2, vertical disk |
| `d10_loaded`/`_friction` … `d14_loaded`/`_friction` | Part 3, off-axis disk at each `d` condition |
| `platform_loaded`, `platform_friction` | Part 3, rotating platform + counterweight only (no disk) |

**Sheet `part1_3_conditions`** (indexed by `group`, matching the column-name prefixes above)

| Column | Unit | Meaning |
|---|---|---|
| `group` | – | Condition identifier — matches the timestamp-column prefix; only used as labels/lookup keys, **the values 10,11,...,14 itself carry no meaning**. You shouldn't change anything in this column. |
| `mass_g` | g | Hanging mass for this condition's f-t fit. |
| `d_cm` | cm | Disk offset distance (Part 3 only; NaN for `hori`/`verti`/`platform`) — drives both the parallel-axis-theorem calculation and the table/figure caption text, so **changing it here is sufficient**, no code edit needed |

**Sheet `part4_raw`**

| Column | Unit | Meaning |
|---|---|---|
| `config` | – | `disk` (calibration disk with known $I$) / `hole1` / `hole2` / `hole3` (black-box suspension holes) |
| `T1_s`, `T2_s` | s | Start/end stopwatch reading; $30T=T_2-T_1$ |

**Parts**

| Part | Meaning |
|---|---|
| 1 | Horizontal disk about its central axis — $I_{\text{obs}}$ vs. $I_{\text{ref}}=MR^2/2$ |
| 2 | Vertical disk about its diameter axis — $I_{\text{obs}}$ vs. $I_{\text{ref}}=MR^2/4+ML^2/12$ |
| 3 | Parallel axis theorem — off-axis disk at five $d$ values plus a no-disk platform baseline; $I_{\text{off}}$ vs. $d^2$ regression compared against the disk+adapter mass (slope) and $I_{\text{hori,ref}}$ (intercept $-\,I_{\text{plat}}$) |
| 4 | Torsion pendulum ($T=2\pi\sqrt{I/\kappa}$) — calibrate $\kappa$ from a disk of known $I$, then solve for the unknown black box's $I_{xx}$, $I_{yy}$, $I_{zz}$ |

First theree parts use $\alpha,\alpha_f=2\pi s$, the slopes of weighted f-t linear fits ($f=0.1/(T_2-T_1)$, $t=(T_1+T_2)/2$, loaded run gives $\alpha$, friction/unloaded run gives $\alpha_f$), fed into $I=\dfrac{mr(g-r\alpha)}{\alpha-\alpha_f}$.

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 14 | `g` | 9.8 | m/s² | Parts 1–3 | Standard gravity |
| 17 | `r_pulley` | 0.0125 | m | Parts 1–3 | step-pulley radius (2r = 25.00 mm) |
| 18 | `M_disk` | 1.4453 | kg | Parts 1–3 | Disk mass |
| 19 | `R_disk` | 0.1142 | m | Parts 1–3 | Disk radius (2R = 22.84 cm) |
| 20 | `L_disk` | 0.0254 | m | Part 2 | Disk thickness |
| 28 | `mass_ref_kg` | 1.4757 | kg | Part 3 | Disk + adapter mass (weighed together, not separately) — compared against the $I_{\text{off}}$–$d^2$ fit's slope |
| 30 | `T_TYPEB_MS` | $1/\sqrt{12}$ | ms | Parts 1–3 | Type-B uncertainty of each raw stopwatch reading $T_1$/$T_2$  (1 ms resolution), propagated into each f-t data point's uncertainty before fitting |
| 33 | `M_ref4` | 0.6606 | kg | Part 4 | Calibration disk mass |
| 34 | `R_ref4` | 0.06015 | m | Part 4 | Calibration disk radius (2R = 120.30 mm) |

Note: remember to include the package `makecell` when compiling raw data latex tables.

### Lab 7 — One-dimensional standing wave

`data.xlsx` has five sheets.

**Sheet `density`** — raw mass/length used to compute the string's linear density

| Column | Unit | Meaning |
|---|---|---|
| `mass_g` | g | Mass of the string segment |
| `length_cm` | cm | Length of the string segment |

**Sheet `resonance`** — raw data of standing wave on the string at different hanging masses

| Column | Unit | Meaning |
|---|---|---|
| `mass_g` | g | Hanging mass, varied across groups; 3 raw rows per mass |
| `n` | – | Antinode count |
| `f_Hz` | Hz | Resonance frequency at that antinode count |
| `halfwidth_cm` | cm | Measured antinode width |

**Sheet `ring`** — raw resonance frequencies of the metal ring

| Column | Unit | Meaning |
|---|---|---|
| `n` | – | Antinode count (3, 4, 5, 7, 9, 11) — only $n=3,5,7,9$ feed the $\log f$-$\log n$/$f$-$n^2$ regressions; $n=4,11$ are held out to test the $f$-$n^2$ fit's prediction |
| `f_Hz` | Hz | Resonance frequency |

**Sheet `spring_K`** — raw scale readings used to measure the spring constant

| Column | Unit | Meaning |
|---|---|---|
| `mass_g` | g | Hanging mass |
| `L_high_cm` | cm | Scale reading at the upper support loop's upper edge |
| `L_low_cm` | cm | Scale reading at the lower support loop's lower edge |

**Sheet `spring_n`** — raw resonance frequencies of the spring's longitudinal standing wave

| Column | Unit | Meaning |
|---|---|---|
| `n` | – | Antinode count (7–11) |
| `f_Hz` | Hz | Resonance frequency |

**Parts**

| Part | Meaning |
|---|---|
| 1 | Transverse wave on a string, formula: $v=\sqrt{T/\mu}$ — wave speed (from resonance $f$ and antinode width $\lambda/2$) vs. tension $T=mg$; the regression slope's reciprocal $\mu_{\text{obs}}$ is compared with $\mu_{\text{ref}}$ measured directly from the string's mass and length |
| 2 | Ring standing wave — $\log f$-$\log n$ slope compared with theoretical exponent 2; $f$-$n^2$ fit used to predict $f$ at $n=4,11$, compared with the actual measurements |
| 3 | Spring longitudinal wave — $F$-$x$ fit gives the spring constant $K$ (Hooke's-law); $f$-$n$ slope compared with $\sqrt{K/M}/2$; $\log f$-$\log n$ slope compared with theoretical exponent 1 |

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 11 | `g` | 9.8 | m/s² | Part 1 and 3 | Standard gravity |
| 12 | `L0_cm` | 7.50 | cm | Part 3 | Spring's natural (unloaded, vertically hung) length |
| 13 | `M_g` | 11.46 | g | Part 3 | Spring mass |

Note: `mu_ref` (Part 1) is computed dynamically from the `density` sheet, not hardcoded.

### Lab 9 — Simple harmonic oscillation

`data.xlsx` has eleven sheets. The four short tabular sheets hold the raw spring-constant data; the seven long sheets hold the recorded motion traces (Tracker exports). Just overwrite the sheets with your own data.

**Sheets `static_small` / `static_medium`** — static (Hooke's-law) method, small-k / medium-k spring

| Column | Unit | Meaning |
|---|---|---|
| `m_g` | g | Hanging mass |
| `x_cm` | cm | Height of the spring's lower support coil above the floor |

**Sheets `dynamic_small` / `dynamic_medium`** — dynamic method, small-k / medium-k spring

| Column | Unit | Meaning |
|---|---|---|
| `m_g` | g | Hanging mass |
| `T1_s` | s | Stopwatch reading at the start of the 1st period |
| `T2_s` | s | Stopwatch reading at the end of the 30th period ($30T=T_2-T_1$) |

**Sheet `energy_yt`** — vertical SHM trace

| Column | Unit | Meaning |
|---|---|---|
| `t_s` | s | Time |
| `y_m` | m | Displacement **about the equilibrium point** |

**Sheet `accel`** — acceleration of the unloaded cart

| Column | Unit | Meaning |
|---|---|---|
| `right_cm_s2` | cm/s² | Acceleration, cart sliding rightward |
| `left_cm_s2` | cm/s² | Acceleration, cart sliding leftward |

**Sheets `friction_level` / `friction_tilt` / `friction_blue`** — cart $x$-$t$ traces (Correspond to Part 3 rail-level, Part 3 rail-tilted, Part 4 equal-mass blue cart, resp.)

| Column | Unit | Meaning |
|---|---|---|
| `t_s` | s | Time |
| `x_m` | m | Displacement |

**Sheets `damped_high` / `damped_low`** — underdamped $x$-$t$ traces, magnet high / low (Part 4)

| Column | Unit | Meaning |
|---|---|---|
| `t_s` | s | Time |
| `x_m` | m | Displacement |

**Parts**

| Part | Meaning |
|---|---|
| 1 | Spring constants — **static** (Hooke's law, slope of $\Delta F$ vs. $\Delta x$ gives $k_s$) and **dynamic** ($T^2$ vs. $m$, slope gives $k_d=4\pi^2/s$, intercept gives the observed effective spring mass $m_{s,\text{obs}}=3b/s$ compared with the weighed $m_{s,\text{ref}}$); done for the low-k and medium-k springs |
| 2 | Energy conservation of the vertical SHM — kinetic/potential/total-energy plots (figures only; constants hardcoded, no error propagation) |
| 3 | SHM under constant friction — detect cart turning points, get period $T$ and the per-cycle amplitude drop $D$, then $f=k_{\text{eff}}D/4$ compared with $f'=M\|a\|$ or $M\|a\|\cos\theta$ (rail level and rail tilted) |
| 4 | Underdamped oscillation — fit $x=Ae^{Et}\cos(Bt+C)+D$ and derive $\omega_0$, damping ratio $\zeta$, mass $M$, damping coefficient $c$, and quasi-period $T$; magnet high, magnet low, and the equal-mass blue cart |

**Choosing which spring constant feeds each calculation**

Part 1 fits four spring constants and leaves them as global variables — `k_s1`, `k_d1` (small-k static/dynamic) and `k_s2`, `k_d2` (medium-k static/dynamic), each with an uncertainty `u_ks1`, `u_kd1`, `u_ks2`, `u_kd2`. Which estimate flows into Parts 2–4 is decided in a few lines of the main flow, all easy to swap:

| Line | Setting | What it controls | How to change |
|---|---|---|---|
| 405 | `energy_part(k_energy=k_s2)` | Part 2 potential energy $PE=\tfrac12 k_{\text{energy}}y^2$ | Pass a different fitted constant (`k_d2`, `k_s1`, …) or a literal number |
| 397, 398 | `K_SMALL, U_SMALL = k_d1, u_kd1` and `K_MEDIUM, U_MEDIUM = k_s2, u_ks2` | the small-k and medium-k estimates that build `k_eff` for Part 3 | Swap each to its static/dynamic counterpart, e.g. `k_s1, u_ks1` or `k_d2, u_kd2` (keep value and uncertainty paired) |
| 399 | `k_eff_two = K_SMALL + K_MEDIUM` | Part 3 effective constant (one small-k + one medium-k spring) | Built from the two lines above; edit the formula only if the spring arrangement differs |
| 401 | `k_eff_dbl = 2 * k_s2` | Part 4 effective constant (two medium springs) | Replace `k_s2` with `k_d2` (and `u_ks2`→`u_kd2` on the next line) to use the dynamic estimate |

The defaults are `k_energy = k_s2`, `k_eff_two = k_d1 + k_s2`,and `k_eff_dbl = 2·k_s2`.

Note: If you change the way `k_eff_two`/`k_eff_dbl` is calculated, please remember to update manually the text on line 445/491.

**Constants** (hardcoded in `lab.py`)

| Line | Variable | Value | Unit | Used in | Meaning |
|---|---|---|---|---|---|
| 43 | `g` | 9.8 | m/s² | Part 1 | Standard gravity, converts hanging mass to restoring force $\Delta F=(m-m_{\text{ref}})g$ |
| 195 | `m` (inside `energy_part`) | 0.08036 | kg | Part 2 | Hanging mass of the vertical oscillator |
| 196 | `m_s` (inside `energy_part`) | 0.01316 | kg | Part 2 | Medium spring mass; its one-third $m_s/3$ is the effective spring mass added in the kinetic energy |
| 381, 384 | `ms_ref_g` | 10.40 / 13.16 | g | Part 1 | Weighed spring masses (small-k / medium-k), compared with $m_{s,\text{obs}}=3b/s$ |
| 408 | `M_BLUE_LEVEL` | 0.58718 | kg | Part 3 (level) | Blue cart + 200 g weights, total oscillating mass $M$ |
| 409 | `M_BLUE_TILT` | 0.59986 | kg | Part 3 (tilt) | Blue cart + 200 g weights + light-blocker, total mass $M'$ |
| 410 | `M_EFF_TWO` | 0.02338 | kg | **Part 3** | Weighed mass of the two springs (one small and one medium), $m_{\text{eff}}$ |
| 411 | `THETA_DEG` | 43.7631 | degree | Part 3 (tilt) | Rail tilt angle $\theta$ ($\delta\theta$ neglected in the propagation) |
| 480 | `M_RED` | 0.44664 | kg | Part 4 (high/low) | Red cart mass |
| 481 | `M_EFF_D` | 0.02632 | kg | **Part 4** | Weighed mass of the two medium springs |
| 530 | `M_BLUE_EQ` | 0.44776 | kg | Part 4 (blue) | Equal-mass blue cart total mass $M''$ |

Notes:
- Turning points (Part 3/4) are found automatically by `turning_points()` via `scipy.signal.find_peaks` with a prominence/spacing heuristic; if your trace is noisier or has a turning point right at the recording boundary, tweak the `prominence`/`distance` settings there.
   - [What is prominence?](https://www.mathworks.com/help/signal/ug/prominence.html)
- The underdamped fits (lines ~487 and ~512) pass an explicit `p0` initial guess tuned to the sample traces. If `curve_fit` fails to converge on your data, update those guesses (especially `B`≈angular frequency and `E`≈decay rate).

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

# About me

## Grades

### 11410 General Physics Laboratory (I)

1st Midterm Exam: 73/100

2nd Midterm Exam: 85/100

Final Exam: 考卷我丟了

Grade: A+

Ranking: 1/32

T Score: 65.22

### 11420 General Physics Laboratory (II)

Midterm Exam: 102/122

Final Exam: 111/120

Grade: A+

Ranking: 1/24

T Score: 72.36

**Details: Report Scores**

| lab | pre | post |
|---|---|---|
| 1 | 93 | 90 |
| 2 | 97 | 90 |
| 3 | 100 | 91 |
| 5 | 100 | 97.14 |
| 7 | 94 | 95 |
| 9 | 97 | 98 |
| 10 | 92 | 95 |
| 13 | 90 | 100 |
| 14 | 94 | 99 |
| 15 | 98 | 93 |
| 17 | 95 | 95 |
| 18 | 93 | 93 |
| 19 | 95 | 100 |
| 21 | 100 | 95 |
| 22 | 90 | 93 |
| 23 | 83 | 86 |

## Project Background and Precautions

The most frustrating thing about *General Physics Laboratory* courses is the disproportionate time invested compared to the final report scores.

Since it's only one credit, we should focus on what's truly important.

After understanding the procedures and formulas used, data processing becomes nothing but repetitive work, so I had ChatGPT write a bunch of reusable Python scripts.

During the summer break, I rewrote these scripts using Claude Code, automating data reading, image storage, and LaTeX table output.

Codex should be able to assist in modifying the code should the lab content change in the future.

I asked Claude Code to package their project processing experience into a skill, which includes my formatting preferences. I may release `SKILL.md` later, which you can then modify to your liking. My goal is to help others replicate this workflow in their labs.

Please note: **DO NOT promote this tool**. If the NTHU General Physics Lab were aware that the existence of such a tool would deprive students of opportunities to train their data analysis skills, they might immediately revise their lab content.