# COBOL Orbital Mechanics Computation System

A comprehensive COBOL program that performs complex orbital mechanics calculations, demonstrating both classical physics principles and mainframe-era programming techniques.

All forumalas are taken from NASA's webiste with thanks. 

## Overview

This program simulates and analyzes spacecraft orbital motion around Earth, including orbit transfers and perturbation effects. It implements fundamental celestial mechanics equations using mainframe-compatible COBOL syntax.

**Computational Intensity:** The program performs extensive iterative calculations including:
- Up to 100,000 iterations for Kepler's equation convergence
- 10,000 orbital propagation time steps simulating ~7 days
- Millions of floating-point operations for trajectory analysis
- Demonstrates COBOL's capability for long-running scientific computations

---

## Physics Concepts and Code Sections

### 1. **Keplerian Orbital Elements**
**Physics:** The classical six orbital elements that uniquely define an orbit in space:
- **Semi-major axis (a)**: Half the longest diameter of the orbital ellipse; determines orbit size and energy
- **Eccentricity (e)**: Shape of the orbit (0 = circular, 0 < e < 1 = elliptical)
- **Inclination (i)**: Tilt of the orbital plane relative to the equator
- **Right Ascension of Ascending Node (RAAN/Ω)**: Where the orbit crosses the equator going northward
- **Argument of Periapsis (ω)**: Orientation of the ellipse within the orbital plane
- **True Anomaly (ν)**: Current position of the satellite along its orbit

**Code Section:** `PRIMARY-ORBITAL-ELEMENTS` (lines 66-73)
```cobol
01  PRIMARY-ORBITAL-ELEMENTS.
    05  WS-SEMI-MAJOR-AXIS      PIC 9(8)V9(6).
    05  WS-ECCENTRICITY         PIC 9V9(15).
    05  WS-INCLINATION          PIC 9(3)V9(12).
    05  WS-RAAN                 PIC 9(3)V9(12).
    05  WS-ARG-PERIAPSIS        PIC 9(3)V9(12).
    05  WS-TRUE-ANOMALY         PIC 9(3)V9(12).
```

**Initialization:** `1000-INITIALIZE-SYSTEM` (lines 194-220)

---

### 2. **Gravitational Parameters**
**Physics:** The standard gravitational parameter (μ = GM) combines the gravitational constant and mass of a celestial body. Used in all orbital calculations:
- **Earth μ**: 398,600.4418 km³/s²
- **Moon μ**: 4,902.8001 km³/s²
- **Sun μ**: 132,712,440,018 km³/s²

The vis-viva equation uses μ to relate orbital velocity to position:
```
v² = μ(2/r - 1/a)
```

**Code Section:** `GRAVITATIONAL-PARAMETERS` (lines 49-58)
```cobol
01  GRAVITATIONAL-PARAMETERS.
    05  WS-MU-EARTH             PIC 9(6)V9(9) VALUE 398600.441800000.
    05  WS-MU-MOON              PIC 9(4)V9(9) VALUE 4902.800066000.
    05  WS-MU-SUN               PIC 9(11)V9(9) VALUE 132712440018.000000000.
```

---

### 3. **Orbital Period (Kepler's Third Law)**
**Physics:** The time it takes for one complete orbit, derived from Kepler's Third Law:
```
T = 2π√(a³/μ)
```
This relates orbital period to semi-major axis. Larger orbits take longer to complete.

**Code Section:** `2000-COMPUTE-PRIMARY-ORBIT` (lines 222-279)
```cobol
*    COMPUTE ORBITAL PERIOD
     COMPUTE WS-TEMP1 = (WS-SEMI-MAJOR-AXIS ** 3) / WS-MU-EARTH.
     PERFORM 2100-COMPUTE-SQUARE-ROOT.
     COMPUTE WS-ORBITAL-PERIOD = WS-TWO-PI * WS-SQRT-VALUE.
```

---

### 4. **Periapsis and Apoapsis**
**Physics:** 
- **Periapsis**: Closest point to Earth in the orbit
- **Apoapsis**: Farthest point from Earth in the orbit

Formulas:
```
r_periapsis = a(1 - e)
r_apoapsis = a(1 + e)
```

At periapsis, the satellite moves fastest; at apoapsis, it moves slowest (conservation of angular momentum).

**Code Section:** `2000-COMPUTE-PRIMARY-ORBIT` (lines 234-237)
```cobol
*    COMPUTE PERIAPSIS AND APOAPSIS RADII
     COMPUTE WS-PERIAPSIS-RADIUS = WS-SEMI-MAJOR-AXIS * (1.0 - WS-ECCENTRICITY).
     COMPUTE WS-APOAPSIS-RADIUS = WS-SEMI-MAJOR-AXIS * (1.0 + WS-ECCENTRICITY).
```

**Velocities Section:** `3000-COMPUTE-DERIVED-PARAMETERS` (lines 306-318)
```cobol
*    COMPUTE VELOCITIES AT PERIAPSIS AND APOAPSIS
     COMPUTE WS-TEMP1 = WS-MU-EARTH * 
         ((2.0 / WS-PERIAPSIS-RADIUS) - (1.0 / WS-SEMI-MAJOR-AXIS)).
```

---

### 5. **Specific Orbital Energy**
**Physics:** The total mechanical energy per unit mass of the orbiting object:
```
ε = -μ/(2a)
```
Negative values indicate bound (elliptical) orbits. Zero energy means parabolic escape trajectory. This is a conserved quantity in two-body problems.

**Code Section:** `2000-COMPUTE-PRIMARY-ORBIT` (lines 246-248)
```cobol
*    COMPUTE SPECIFIC ORBITAL ENERGY
     COMPUTE WS-SPECIFIC-ENERGY = 
         (0.0 - WS-MU-EARTH) / (2.0 * WS-SEMI-MAJOR-AXIS).
```

---

### 6. **Angular Momentum**
**Physics:** The rotational momentum of the orbit, conserved in Keplerian motion:
```
h = √(μp)
where p = a(1 - e²) is the semi-latus rectum
```
Angular momentum determines the "flatness" of the orbit and is perpendicular to the orbital plane.

**Code Section:** `2000-COMPUTE-PRIMARY-ORBIT` (lines 252-255)
```cobol
*    COMPUTE ANGULAR MOMENTUM
     COMPUTE WS-TEMP1 = WS-MU-EARTH * WS-SEMI-LATUS-RECTUM.
     PERFORM 2100-COMPUTE-SQUARE-ROOT.
     MOVE WS-SQRT-VALUE TO WS-ANGULAR-MOMENTUM.
```

---

### 7. **Kepler's Equation**
**Physics:** Relates position in orbit (true anomaly) to time through the mean and eccentric anomalies:
```
M = E - e sin(E)
```
Where:
- **M** = Mean anomaly (proportional to time)
- **E** = Eccentric anomaly (geometric parameter)
- **e** = Eccentricity

This transcendental equation has no closed-form solution and must be solved iteratively using Newton-Raphson method:
```
E(n+1) = E(n) - (E - e sin(E) - M)/(1 - e cos(E))
```

The program allows up to 100,000 iterations for high-precision convergence, demonstrating the computational intensity of orbital mechanics calculations.

**Code Section:** `4000-SOLVE-KEPLER-EQUATION` (lines 363-408)
```cobol
*    SOLVE KEPLER'S EQUATION ITERATIVELY
     PERFORM UNTIL WS-CONVERGENCE-FLAG = 1 OR
                   WS-KEPLER-ITERATIONS >= WS-KEPLER-MAX-ITER
         COMPUTE WS-TEMP1 = WS-ECCENTRIC-ANOMALY -
             (WS-ECCENTRICITY * FUNCTION SIN(WS-ECCENTRIC-ANOMALY)) -
             WS-MEAN-ANOMALY
         
         COMPUTE WS-TEMP2 = 1.0 - 
             (WS-ECCENTRICITY * FUNCTION COS(WS-ECCENTRIC-ANOMALY))
         
         COMPUTE WS-TEMP3 = WS-ECCENTRIC-ANOMALY - (WS-TEMP1 / WS-TEMP2)
```

---

### 8. **State Vectors (Position & Velocity)**
**Physics:** Complete description of orbital state at any moment:
- **Position vector (r)**: Location in 3D space [x, y, z]
- **Velocity vector (v)**: Speed and direction [vx, vy, vz]

Computed from eccentric anomaly in the orbital plane:
```
x = a(cos(E) - e)
y = a√(1-e²) sin(E)
vx = -√(μa)/r sin(E)
vy = √(μa)/r √(1-e²) cos(E)
```

These vectors can be rotated to the inertial frame using the orientation angles (i, Ω, ω).

**Code Section:** `5000-COMPUTE-STATE-VECTORS` (lines 410-455)
```cobol
*    COMPUTE POSITION IN ORBITAL PLANE
     COMPUTE WS-POSITION-X = WS-SEMI-MAJOR-AXIS * 
         (FUNCTION COS(WS-TEMP1) - WS-ECCENTRICITY).
     COMPUTE WS-POSITION-Y = WS-SEMI-MAJOR-AXIS * 
         FUNCTION SQRT(1.0 - (WS-ECCENTRICITY ** 2)) *
         FUNCTION SIN(WS-TEMP1).
     
*    COMPUTE VELOCITY IN ORBITAL PLANE
     COMPUTE WS-VELOCITY-X = 0.0 - (WS-TEMP2 * FUNCTION SIN(WS-TEMP1)).
     COMPUTE WS-VELOCITY-Y = WS-TEMP2 * 
         FUNCTION SQRT(1.0 - (WS-ECCENTRICITY ** 2)) *
         FUNCTION COS(WS-TEMP1).
```

---

### 9. **Hohmann Transfer Orbit**
**Physics:** The most fuel-efficient two-impulse maneuver to transfer between two circular coplanar orbits.

Process:
1. Apply **Δv₁** at initial orbit to enter elliptical transfer orbit
2. Coast along transfer orbit (half period)
3. Apply **Δv₂** at target orbit to circularize

Formulas:
```
a_transfer = (r₁ + r₂)/2
Δv₁ = √(μ/r₁) * (√(2r₂/(r₁+r₂)) - 1)
Δv₂ = √(μ/r₂) * (1 - √(2r₁/(r₁+r₂)))
t_transfer = π√(a³_transfer/μ)
```

The phase angle is calculated to ensure the target is at the correct position when the spacecraft arrives.

**Code Section:** `6000-COMPUTE-TRANSFER-ORBIT` (lines 457-529)
```cobol
*    COMPUTE TRANSFER ORBIT SEMI-MAJOR AXIS
     COMPUTE WS-TRANSFER-SMA = (WS-INITIAL-ORBIT-R + WS-TARGET-ORBIT-R) / 2.0.
     
*    COMPUTE FIRST DELTA-V
     COMPUTE WS-DELTA-V1 = WS-TEMP3 - WS-TEMP2.
     
*    COMPUTE SECOND DELTA-V
     COMPUTE WS-DELTA-V2 = WS-TEMP2 - WS-TEMP3.
     
*    COMPUTE TOTAL DELTA-V
     COMPUTE WS-TOTAL-DELTA-V = WS-DELTA-V1 + WS-DELTA-V2.
     
*    COMPUTE TRANSFER TIME
     COMPUTE WS-TRANSFER-TIME = WS-PI * WS-SQRT-VALUE.
```

Example in program: LEO (400 km altitude) to GEO (35,786 km altitude).

---

### 10. **J2 Perturbations (Oblateness Effect)**
**Physics:** Earth is not a perfect sphere—it bulges at the equator. The J2 coefficient (≈0.00108263) quantifies this oblateness and causes two main orbital perturbations:

**Nodal Precession (Regression of Nodes):**
```
dΩ/dt = -1.5 n J₂ (Rₑ/a)² cos(i)/(1-e²)²
```
The orbital plane rotates around Earth's axis. Retrograde for prograde orbits (i < 90°).

**Apsidal Precession (Rotation of Line of Apsides):**
```
dω/dt = 0.75 n J₂ (Rₑ/a)² (5cos²(i)-1)/(1-e²)²
```
The orientation of the ellipse within its plane rotates. Zero for critical inclination (≈63.4°).

**Code Section:** `7000-ANALYZE-PERTURBATIONS` (lines 531-586)
```cobol
*    COMPUTE J2 PERTURBATION EFFECTS
*    NODAL PRECESSION RATE
     COMPUTE WS-NODAL-PRECESSION = 
         -1.5 * WS-J2-COEFFICIENT * WS-MEAN-MOTION *
         ((WS-EARTH-RADIUS / WS-SEMI-MAJOR-AXIS) ** 2) *
         WS-TEMP3 / (WS-TEMP1 ** 2).
     
*    APSIDAL PRECESSION RATE
     COMPUTE WS-APSIDAL-PRECESSION = 
         0.75 * WS-J2-COEFFICIENT * WS-MEAN-MOTION *
         ((WS-EARTH-RADIUS / WS-SEMI-MAJOR-AXIS) ** 2) *
         WS-TEMP4 / (WS-TEMP1 ** 2).
```

---

### 11. **Third-Body Perturbations**
**Physics:** The Sun and Moon exert gravitational forces that perturb Earth orbits. These are especially important for high-altitude orbits (GEO, lunar missions).

The perturbing acceleration is proportional to:
```
a_perturb ≈ (μ_body/d³) * r_satellite
```
Where d is the distance to the perturbing body.

**Effects:**
- **Solar perturbations**: Cause long-term drift in semi-major axis and eccentricity
- **Lunar perturbations**: Similar effects but with monthly period instead of yearly

**Code Section:** `7000-ANALYZE-PERTURBATIONS` (lines 568-580)
```cobol
*    COMPUTE THIRD-BODY PERTURBATIONS
*    SOLAR PERTURBATION (SIMPLIFIED)
     COMPUTE WS-TEMP1 = WS-MU-SUN / (WS-MOON-DISTANCE ** 3).
     COMPUTE WS-SOLAR-PERTURBATION = WS-TEMP1 * WS-RADIUS-MAGNITUDE.
     
*    LUNAR PERTURBATION (SIMPLIFIED)
     COMPUTE WS-TEMP1 = WS-MU-MOON / (WS-MOON-DISTANCE ** 3).
     COMPUTE WS-LUNAR-PERTURBATION = WS-TEMP1 * WS-RADIUS-MAGNITUDE.
```

---

### 12. **Orbit Propagation**
**Physics:** Numerical integration of the equations of motion to predict future orbital positions. The program uses a simplified propagation that:

1. Updates mean anomaly based on mean motion
2. Applies perturbation effects to RAAN and argument of periapsis
3. Solves Kepler's equation for new position
4. Updates state vectors

This is a simplified numerical propagator. Real systems use more sophisticated integrators (Runge-Kutta, Cowell's method) that account for all perturbations simultaneously.

**Time steps:** The program uses 60-second steps to propagate the orbit over approximately 166.67 hours (10,000 iterations), demonstrating long-running computational capability.

**Code Section:** `8000-PROPAGATE-ORBIT` (lines 588-664)
```cobol
PERFORM UNTIL WS-CURRENT-TIME >= WS-SIMULATION-DURATION
    *    UPDATE MEAN ANOMALY
         COMPUTE WS-MEAN-ANOMALY = WS-MEAN-ANOMALY + 
             (WS-MEAN-MOTION * WS-TIME-STEP * WS-RAD-TO-DEG).
         
    *    APPLY PERTURBATIONS
         COMPUTE WS-RAAN = WS-RAAN + (WS-NODAL-PRECESSION * WS-TIME-STEP).
         COMPUTE WS-ARG-PERIAPSIS = WS-ARG-PERIAPSIS + 
             (WS-APSIDAL-PRECESSION * WS-TIME-STEP).
         
    *    SOLVE KEPLER'S EQUATION FOR NEW TIME STEP
         [... iterative solution ...]
         
    *    UPDATE STATE VECTORS
         [... position and velocity recalculation ...]
END-PERFORM.
```

---

### 13. **Mathematical Utilities**

#### Newton-Raphson Square Root
**Physics/Math:** Iterative method to compute √x by repeatedly improving a guess:
```
x(n+1) = (x(n) + a/x(n))/2
```

**Code Section:** `2100-COMPUTE-SQUARE-ROOT` (lines 281-304)
```cobol
PERFORM UNTIL WS-SQRT-ITERATIONS > 20
    COMPUTE WS-TEMP2 = (WS-SQRT-GUESS + (WS-TEMP1 / WS-SQRT-GUESS)) / 2.0
    IF FUNCTION ABS(WS-TEMP2 - WS-SQRT-GUESS) < WS-TOLERANCE
        MOVE WS-TEMP2 TO WS-SQRT-VALUE
        GO TO 2100-EXIT
    END-IF
```

#### Sine and Cosine (Taylor Series)
**Math:** Approximate trigonometric functions using power series:
```
sin(x) = x - x³/3! + x⁵/5! - x⁷/7! + ...
cos(x) = 1 - x²/2! + x⁴/4! - x⁶/6! + ...
```

**Code Section:** `3100-COMPUTE-SINE-COSINE` (lines 337-361)
```cobol
*    COMPUTE SINE USING TAYLOR SERIES (FIRST 5 TERMS)
     COMPUTE WS-SINE-VALUE = WS-ANGLE-WORK
         - ((WS-ANGLE-WORK ** 3) / 6.0)
         + ((WS-ANGLE-WORK ** 5) / 120.0)
         - ((WS-ANGLE-WORK ** 7) / 5040.0)
         + ((WS-ANGLE-WORK ** 9) / 362880.0).
```

---

## Program Structure

### Data Division Sections

1. **Mathematical Constants** (lines 41-47): π, 2π, conversion factors
2. **Gravitational Parameters** (lines 49-58): μ values for celestial bodies
3. **Orbital Elements** (lines 66-73): Six Keplerian elements
4. **Derived Parameters** (lines 78-87): Period, energy, velocities, etc.
5. **State Vectors** (lines 92-99): Position and velocity components
6. **Transfer Orbit Parameters** (lines 104-111): Hohmann transfer data
7. **Perturbation Variables** (lines 116-125): J2, drag, third-body effects
8. **Control Variables** (lines 130-139): Iteration counters, time steps
9. **Temporary Variables** (lines 144-154): Calculation workspace
10. **Report Variables** (lines 159-189): Output formatting

### Procedure Division Flow

```
MAIN-CONTROL (0000)
├── INITIALIZE-SYSTEM (1000)
├── COMPUTE-PRIMARY-ORBIT (2000)
│   └── COMPUTE-SQUARE-ROOT (2100)
├── COMPUTE-DERIVED-PARAMETERS (3000)
│   └── COMPUTE-SINE-COSINE (3100)
├── SOLVE-KEPLER-EQUATION (4000)
├── COMPUTE-STATE-VECTORS (5000)
├── COMPUTE-TRANSFER-ORBIT (6000)
├── ANALYZE-PERTURBATIONS (7000)
├── PROPAGATE-ORBIT (8000)
├── GENERATE-REPORT (9000)
└── TERMINATE-PROGRAM (9999)
```

---

## Compilation and Execution

### Prerequisites
- GnuCOBOL (cobc) compiler installed
- POSIX-compatible system (Linux, macOS, Unix)

### Compile
```bash
cobc -x orbital-mechanics.cob -o orbital-mechanics
```

Flags:
- `-x`: Create executable

### Run
```bash
./orbital-mechanics
```

### Output
1. **Console output**: Real-time progress and key calculations
2. **Report file**: `ORBITRPT.TXT` - Comprehensive formatted report

---

## Example Scenario

The program simulates a **LEO to GEO transfer mission with extended orbit propagation**:

- **Initial Orbit**: Low Earth Orbit at 400 km altitude
  - Semi-major axis: 6,778.137 km
  - Period: ~92.6 minutes
  - Velocity: ~7.67 km/s

- **Target Orbit**: Geostationary Earth Orbit at 35,786 km altitude
  - Semi-major axis: 42,164 km
  - Period: ~24 hours
  - Velocity: ~3.07 km/s

- **Transfer Maneuver**:
  - Δv₁ ≈ 2.44 km/s (LEO to transfer orbit)
  - Δv₂ ≈ 1.47 km/s (transfer to GEO)
  - Total Δv ≈ 3.91 km/s
  - Transfer time ≈ 5.25 hours

- **Orbit Propagation**:
  - Simulation duration: 600,000 seconds (~166.67 hours / ~6.94 days)
  - Time steps: 60 seconds
  - Total iterations: 10,000
  - Demonstrates long-running computational capability of COBOL for scientific applications

---

## Physics References

1. **Vallado, D.A.** (2013). *Fundamentals of Astrodynamics and Applications*. Microcosm Press.
2. **Bate, R.R., Mueller, D.D., White, J.E.** (1971). *Fundamentals of Astrodynamics*. Dover.
3. **Curtis, H.D.** (2013). *Orbital Mechanics for Engineering Students*. Butterworth-Heinemann.
4. **Prussing, J.E., Conway, B.A.** (2012). *Orbital Mechanics*. Oxford University Press.

---

## Technical Notes

### Numerical Precision
- Uses COBOL `PIC 9V9(15)` for high-precision floating-point (15 decimal places)
- Convergence tolerance: 10⁻¹⁵
- Suitable for orbital mechanics requiring km-level accuracy
- Maximum iterations: 100,000 for Kepler's equation solver
- Orbit propagation: 10,000 time steps over ~7 days of simulated time
- Total computational operations: millions of floating-point calculations

### Coordinate Systems
- Orbital plane coordinates (perifocal frame)
- Could be extended to inertial frame (ECI) with rotation matrices

### Limitations
- Two-body problem assumption (primary perturbations included but simplified)
- Constant perturbation rates (actual rates vary with position)
- No atmospheric drag model integration
- Simplified third-body perturbations (actual geometry not modeled)

### Potential Extensions
- Full ECI coordinate transformations
- Runge-Kutta 4th-order integrator
- Ground track visualization
- Eclipse prediction
- Station keeping analysis
- Multi-revolution Lambert solver

---

## License

This code is provided for educational purposes to demonstrate both orbital mechanics principles and COBOL programming capabilities. All rights reserved by moshix


  
Moshix   
November 2025

