       IDENTIFICATION DIVISION.
       PROGRAM-ID. ORBITAL-MECHANICS.
       AUTHOR. ORBITAL DYNAMICS LABORATORY.
       DATE-WRITTEN. 2025-11-14.
      *****************************************************************
      * COPYRIGHT 2025 BY MOSHIX. ALL RIGHTS RESERVED.                *
      *****************************************************************
      * ORBITAL MECHANICS COMPUTATION SYSTEM                          *
      * PERFORMS COMPLEX ORBITAL CALCULATIONS INCLUDING:              *
      * - KEPLERIAN ORBITAL ELEMENTS                                  *
      * - ORBITAL VELOCITY AND ENERGY CALCULATIONS                    *
      * - KEPLER'S EQUATION SOLVER (ITERATIVE)                        *
      * - POSITION AND VELOCITY VECTORS IN ORBITAL PLANE              *
      * - HOHMANN TRANSFER ORBIT CALCULATIONS                         *
      * - MULTI-BODY PERTURBATION ANALYSIS                            *
      * - LONG-DURATION ORBITAL PROPAGATION                           *
      *****************************************************************
       
       ENVIRONMENT DIVISION.
       CONFIGURATION SECTION.
       SOURCE-COMPUTER. IBM-370.
       OBJECT-COMPUTER. IBM-370.
       
       INPUT-OUTPUT SECTION.
       FILE-CONTROL.
           SELECT ORBIT-REPORT-FILE ASSIGN TO "ORBITRPT.TXT"
               ORGANIZATION IS LINE SEQUENTIAL.
       
       DATA DIVISION.
       FILE SECTION.
       FD  ORBIT-REPORT-FILE
           RECORDING MODE IS F
           BLOCK CONTAINS 0 RECORDS.
       01  ORBIT-REPORT-RECORD         PIC X(132).
       
       WORKING-STORAGE SECTION.
      *****************************************************************
      * MATHEMATICAL CONSTANTS                                        *
      *****************************************************************
       01  MATHEMATICAL-CONSTANTS.
           05  WS-PI                   PIC 9(1)V9(15) VALUE 
               3.141592653589793.
           05  WS-TWO-PI               PIC 9(2)V9(15) VALUE 
               6.283185307179586.
           05  WS-HALF-PI              PIC 9(1)V9(15) VALUE 
               1.570796326794897.
           05  WS-DEG-TO-RAD           PIC 9(1)V9(15) VALUE 
               0.017453292519943.
           05  WS-RAD-TO-DEG           PIC 9(2)V9(15) VALUE 
               57.295779513082320.
           05  WS-TOLERANCE            PIC 9V9(15) VALUE 
               0.000000000000001.
       
      *****************************************************************
      * GRAVITATIONAL PARAMETERS (KM^3/S^2)                           *
      *****************************************************************
       01  GRAVITATIONAL-PARAMETERS.
           05  WS-MU-EARTH             PIC 9(6)V9(9) VALUE 
               398600.441800000.
           05  WS-MU-MOON              PIC 9(4)V9(9) VALUE 
               4902.800066000.
           05  WS-MU-SUN               PIC 9(12)V9(9) VALUE 
               132712440018.000000000.
           05  WS-EARTH-RADIUS         PIC 9(4)V9(6) VALUE 
               6378.137000.
           05  WS-MOON-DISTANCE        PIC 9(6)V9(6) VALUE 
               384400.000000.
       
      *****************************************************************
      * ORBITAL ELEMENTS - PRIMARY ORBIT                              *
      *****************************************************************
       01  PRIMARY-ORBITAL-ELEMENTS.
           05  WS-SEMI-MAJOR-AXIS      PIC 9(8)V9(6).
           05  WS-ECCENTRICITY         PIC 9V9(15).
           05  WS-INCLINATION          PIC 9(3)V9(12).
           05  WS-RAAN                 PIC 9(3)V9(12).
           05  WS-ARG-PERIAPSIS        PIC 9(3)V9(12).
           05  WS-TRUE-ANOMALY         PIC 9(3)V9(12).
           05  WS-MEAN-ANOMALY         PIC 9(3)V9(12).
           05  WS-ECCENTRIC-ANOMALY    PIC 9(3)V9(12).
       
      *****************************************************************
      * ORBITAL PARAMETERS - DERIVED                                  *
      *****************************************************************
       01  DERIVED-ORBITAL-PARAMETERS.
           05  WS-ORBITAL-PERIOD       PIC 9(8)V9(6).
           05  WS-MEAN-MOTION          PIC 9(3)V9(12).
           05  WS-PERIAPSIS-RADIUS     PIC 9(8)V9(6).
           05  WS-APOAPSIS-RADIUS      PIC 9(8)V9(6).
           05  WS-PERIAPSIS-VELOCITY   PIC 9(5)V9(10).
           05  WS-APOAPSIS-VELOCITY    PIC 9(5)V9(10).
           05  WS-SPECIFIC-ENERGY      PIC S9(5)V9(10).
           05  WS-ANGULAR-MOMENTUM     PIC 9(10)V9(6).
           05  WS-SEMI-LATUS-RECTUM    PIC 9(8)V9(6).
       
      *****************************************************************
      * POSITION AND VELOCITY VECTORS                                 *
      *****************************************************************
       01  STATE-VECTORS.
           05  WS-POSITION-X           PIC S9(8)V9(6).
           05  WS-POSITION-Y           PIC S9(8)V9(6).
           05  WS-POSITION-Z           PIC S9(8)V9(6).
           05  WS-VELOCITY-X           PIC S9(5)V9(10).
           05  WS-VELOCITY-Y           PIC S9(5)V9(10).
           05  WS-VELOCITY-Z           PIC S9(5)V9(10).
           05  WS-RADIUS-MAGNITUDE     PIC 9(8)V9(6).
           05  WS-VELOCITY-MAGNITUDE   PIC 9(5)V9(10).
       
      *****************************************************************
      * HOHMANN TRANSFER ORBIT PARAMETERS                             *
      *****************************************************************
       01  TRANSFER-ORBIT-PARAMETERS.
           05  WS-INITIAL-ORBIT-R      PIC 9(8)V9(6).
           05  WS-TARGET-ORBIT-R       PIC 9(8)V9(6).
           05  WS-TRANSFER-SMA         PIC 9(8)V9(6).
           05  WS-DELTA-V1             PIC 9(5)V9(10).
           05  WS-DELTA-V2             PIC 9(5)V9(10).
           05  WS-TOTAL-DELTA-V        PIC 9(5)V9(10).
           05  WS-TRANSFER-TIME        PIC 9(8)V9(6).
           05  WS-PHASE-ANGLE          PIC 9(3)V9(12).
       
      *****************************************************************
      * PERTURBATION ANALYSIS VARIABLES                               *
      *****************************************************************
       01  PERTURBATION-VARIABLES.
           05  WS-J2-COEFFICIENT       PIC 9V9(15) VALUE 
               0.001082630000000.
           05  WS-NODAL-PRECESSION     PIC S9(3)V9(12).
           05  WS-APSIDAL-PRECESSION   PIC S9(3)V9(12).
           05  WS-SOLAR-PERTURBATION   PIC S9(5)V9(10).
           05  WS-LUNAR-PERTURBATION   PIC S9(5)V9(10).
           05  WS-DRAG-COEFFICIENT     PIC 9V9(15) VALUE 
               0.000000100000000.
           05  WS-ATMOSPHERIC-DENSITY  PIC 9V9(15) VALUE 
               0.000000000001000.
       
      *****************************************************************
      * ITERATION AND COMPUTATION CONTROL                             *
      *****************************************************************
       01  COMPUTATION-CONTROL.
           05  WS-ITERATION-COUNTER    PIC 9(10) VALUE ZERO.
           05  WS-MAX-ITERATIONS       PIC 9(10) VALUE 100000000.
           05  WS-TIME-STEP            PIC 9(5)V9(6) VALUE 60.000000.
           05  WS-CURRENT-TIME         PIC 9(12)V9(6) VALUE ZERO.
           05  WS-SIMULATION-DURATION  PIC 9(12)V9(6) VALUE 
               600000.000000.
           05  WS-KEPLER-ITERATIONS    PIC 9(7) VALUE ZERO.
           05  WS-KEPLER-MAX-ITER      PIC 9(7) VALUE 100000.
           05  WS-CONVERGENCE-FLAG     PIC 9 VALUE ZERO.
       
      *****************************************************************
      * TEMPORARY CALCULATION VARIABLES                               *
      *****************************************************************
       01  TEMP-CALCULATION-VARS.
           05  WS-TEMP1                PIC S9(15)V9(15).
           05  WS-TEMP2                PIC S9(15)V9(15).
           05  WS-TEMP3                PIC S9(15)V9(15).
           05  WS-TEMP4                PIC S9(15)V9(15).
           05  WS-TEMP5                PIC S9(15)V9(15).
           05  WS-SINE-VALUE           PIC S9V9(15).
           05  WS-COSINE-VALUE         PIC S9V9(15).
           05  WS-SQRT-VALUE           PIC 9(15)V9(15).
           05  WS-SQRT-GUESS           PIC 9(15)V9(15).
           05  WS-SQRT-ITERATIONS      PIC 9(3).
           05  WS-ANGLE-WORK           PIC S9(5)V9(15).
       
      *****************************************************************
      * REPORT FORMATTING VARIABLES                                   *
      *****************************************************************
       01  REPORT-VARIABLES.
           05  WS-LINE-COUNTER         PIC 9(5) VALUE ZERO.
           05  WS-PAGE-COUNTER         PIC 9(5) VALUE 1.
           05  WS-LINES-PER-PAGE       PIC 9(3) VALUE 60.
       
       01  REPORT-HEADER-1.
           05  FILLER                  PIC X(45) VALUE SPACES.
           05  FILLER                  PIC X(42) VALUE
               "ORBITAL MECHANICS COMPUTATION REPORT".
           05  FILLER                  PIC X(45) VALUE SPACES.
       
       01  REPORT-HEADER-2.
           05  FILLER                  PIC X(10) VALUE "PAGE: ".
           05  RPT-PAGE-NUM            PIC Z,ZZ9.
           05  FILLER                  PIC X(114) VALUE SPACES.
       
       01  REPORT-SECTION-HEADER.
           05  FILLER                  PIC X(5) VALUE SPACES.
           05  RPT-SECTION-NAME        PIC X(60).
           05  FILLER                  PIC X(67) VALUE SPACES.
       
       01  REPORT-DETAIL-LINE.
           05  FILLER                  PIC X(10) VALUE SPACES.
           05  RPT-PARAMETER-NAME      PIC X(40).
           05  FILLER                  PIC X(5) VALUE " = ".
           05  RPT-PARAMETER-VALUE     PIC X(25).
           05  FILLER                  PIC X(5) VALUE SPACES.
           05  RPT-PARAMETER-UNIT      PIC X(20).
           05  FILLER                  PIC X(27) VALUE SPACES.
       
       01  REPORT-BLANK-LINE           PIC X(132) VALUE SPACES.
       
       01  REPORT-SEPARATOR-LINE.
           05  FILLER                  PIC X(5) VALUE SPACES.
           05  FILLER                  PIC X(80) VALUE ALL "=".
           05  FILLER                  PIC X(47) VALUE SPACES.
       
      *****************************************************************
      * DISPLAY FORMATTING VARIABLES                                  *
      *****************************************************************
       01  DISPLAY-FORMATS.
           05  DISP-NUMERIC-1          PIC Z,ZZZ,ZZZ,ZZZ,ZZ9.999999.
           05  DISP-NUMERIC-2          PIC Z,ZZZ,ZZZ,ZZ9.999999999999.
           05  DISP-NUMERIC-3          PIC Z,ZZZ,ZZZ,ZZ9.999.
       
       PROCEDURE DIVISION.
       
      *****************************************************************
      * MAIN PROGRAM CONTROL                                          *
      *****************************************************************
       0000-MAIN-CONTROL.
           DISPLAY "========================================".
           DISPLAY "ORBITAL MECHANICS COMPUTATION SYSTEM".
           DISPLAY "INITIALIZING COMPLEX ORBITAL ANALYSIS".
           DISPLAY "(c) 2025 by moshix. All rights reserved".
           DISPLAY "========================================".
           DISPLAY " ".
           
           PERFORM 1000-INITIALIZE-SYSTEM THRU 1000-EXIT.
           PERFORM 2000-COMPUTE-PRIMARY-ORBIT THRU 2000-EXIT.
           PERFORM 3000-COMPUTE-DERIVED-PARAMETERS THRU 3000-EXIT.
           PERFORM 4000-SOLVE-KEPLER-EQUATION THRU 4000-EXIT.
           PERFORM 5000-COMPUTE-STATE-VECTORS THRU 5000-EXIT.
           PERFORM 6000-COMPUTE-TRANSFER-ORBIT THRU 6000-EXIT.
           PERFORM 7000-ANALYZE-PERTURBATIONS THRU 7000-EXIT.
           PERFORM 8000-PROPAGATE-ORBIT THRU 8000-EXIT.
           PERFORM 9000-GENERATE-REPORT THRU 9000-EXIT.
           PERFORM 9999-TERMINATE-PROGRAM.
           
           STOP RUN.
       
      *****************************************************************
      * INITIALIZE ORBITAL PARAMETERS                                 *
      *****************************************************************
       1000-INITIALIZE-SYSTEM.
           DISPLAY "INITIALIZING ORBITAL PARAMETERS...".
           
      *    INITIALIZE PRIMARY ORBITAL ELEMENTS
      *    LOW EARTH ORBIT (LEO) TO GEOSTATIONARY ORBIT (GEO) SCENARIO
           COMPUTE WS-SEMI-MAJOR-AXIS = 
               WS-EARTH-RADIUS + 400.000000.
           MOVE 0.001500000000000 TO WS-ECCENTRICITY.
           MOVE 28.500000000000 TO WS-INCLINATION.
           MOVE 45.000000000000 TO WS-RAAN.
           MOVE 30.000000000000 TO WS-ARG-PERIAPSIS.
           MOVE 0.000000000000 TO WS-TRUE-ANOMALY.
           
      *    INITIALIZE TRANSFER ORBIT PARAMETERS
           COMPUTE WS-INITIAL-ORBIT-R = WS-SEMI-MAJOR-AXIS.
           COMPUTE WS-TARGET-ORBIT-R = 
               WS-EARTH-RADIUS + 35786.000000.
           
           DISPLAY "  SEMI-MAJOR AXIS: " WS-SEMI-MAJOR-AXIS " KM".
           DISPLAY "  ECCENTRICITY: " WS-ECCENTRICITY.
           DISPLAY "  INCLINATION: " WS-INCLINATION " DEG".
           DISPLAY "INITIALIZATION COMPLETE.".
           DISPLAY " ".
       1000-EXIT.
           EXIT.
       
      *****************************************************************
      * COMPUTE PRIMARY ORBITAL CHARACTERISTICS                       *
      *****************************************************************
       2000-COMPUTE-PRIMARY-ORBIT.
           DISPLAY "COMPUTING PRIMARY ORBITAL ELEMENTS...".
           
      *    COMPUTE PERIAPSIS AND APOAPSIS RADII
           COMPUTE WS-PERIAPSIS-RADIUS = 
               WS-SEMI-MAJOR-AXIS * (1.0 - WS-ECCENTRICITY).
           COMPUTE WS-APOAPSIS-RADIUS = 
               WS-SEMI-MAJOR-AXIS * (1.0 + WS-ECCENTRICITY).
           
           DISPLAY "  PERIAPSIS RADIUS: " WS-PERIAPSIS-RADIUS " KM".
           DISPLAY "  APOAPSIS RADIUS: " WS-APOAPSIS-RADIUS " KM".
           
      *    COMPUTE SEMI-LATUS RECTUM
           COMPUTE WS-TEMP1 = 1.0 - (WS-ECCENTRICITY ** 2).
           COMPUTE WS-SEMI-LATUS-RECTUM = 
               WS-SEMI-MAJOR-AXIS * WS-TEMP1.
           
      *    COMPUTE SPECIFIC ORBITAL ENERGY
           COMPUTE WS-SPECIFIC-ENERGY = 
               (0.0 - WS-MU-EARTH) / (2.0 * WS-SEMI-MAJOR-AXIS).
           
           DISPLAY "  SPECIFIC ENERGY: " WS-SPECIFIC-ENERGY 
               " KM^2/S^2".
           
      *    COMPUTE ANGULAR MOMENTUM
           COMPUTE WS-TEMP1 = WS-MU-EARTH * WS-SEMI-LATUS-RECTUM.
           PERFORM 2100-COMPUTE-SQUARE-ROOT THRU 2100-EXIT.
           MOVE WS-SQRT-VALUE TO WS-ANGULAR-MOMENTUM.
           
           DISPLAY "  ANGULAR MOMENTUM: " WS-ANGULAR-MOMENTUM 
               " KM^2/S".
           
      *    COMPUTE ORBITAL PERIOD
           COMPUTE WS-TEMP1 = (WS-SEMI-MAJOR-AXIS ** 3) / WS-MU-EARTH.
           PERFORM 2100-COMPUTE-SQUARE-ROOT THRU 2100-EXIT.
           COMPUTE WS-ORBITAL-PERIOD = WS-TWO-PI * WS-SQRT-VALUE.
           
           DISPLAY "  ORBITAL PERIOD: " WS-ORBITAL-PERIOD " SECONDS".
           MOVE WS-ORBITAL-PERIOD TO DISP-NUMERIC-3.
           DISPLAY "                 (" DISP-NUMERIC-3 " SEC)".
           
      *    COMPUTE MEAN MOTION
           COMPUTE WS-MEAN-MOTION = WS-TWO-PI / WS-ORBITAL-PERIOD.
           
           DISPLAY "PRIMARY ORBIT COMPUTATION COMPLETE.".
           DISPLAY " ".
       2000-EXIT.
           EXIT.
       
      *****************************************************************
      * COMPUTE SQUARE ROOT USING NEWTON-RAPHSON METHOD               *
      *****************************************************************
       2100-COMPUTE-SQUARE-ROOT.
           IF WS-TEMP1 <= 0.0
               MOVE 0.0 TO WS-SQRT-VALUE
               GO TO 2100-EXIT
           END-IF.
           
           MOVE WS-TEMP1 TO WS-SQRT-GUESS.
           MOVE 0 TO WS-SQRT-ITERATIONS.
           
           PERFORM 2110-SQRT-ITERATION THRU 2110-DONE
               UNTIL WS-SQRT-ITERATIONS > 20.
           
           MOVE WS-SQRT-GUESS TO WS-SQRT-VALUE.
       2100-EXIT.
           EXIT.
       
      *****************************************************************
      * SQUARE ROOT ITERATION HELPER                                  *
      *****************************************************************
       2110-SQRT-ITERATION.
           COMPUTE WS-TEMP2 = 
               (WS-SQRT-GUESS + (WS-TEMP1 / WS-SQRT-GUESS)) / 2.0.
           IF FUNCTION ABS(WS-TEMP2 - WS-SQRT-GUESS) < WS-TOLERANCE
               MOVE WS-TEMP2 TO WS-SQRT-VALUE
               MOVE 21 TO WS-SQRT-ITERATIONS
               GO TO 2110-DONE
           END-IF.
           MOVE WS-TEMP2 TO WS-SQRT-GUESS.
           ADD 1 TO WS-SQRT-ITERATIONS.
       2110-DONE.
           EXIT.
       
      *****************************************************************
      * COMPUTE DERIVED ORBITAL PARAMETERS                            *
      *****************************************************************
       3000-COMPUTE-DERIVED-PARAMETERS.
           DISPLAY "COMPUTING DERIVED PARAMETERS...".
           
      *    COMPUTE VELOCITIES AT PERIAPSIS AND APOAPSIS
           COMPUTE WS-PERIAPSIS-VELOCITY = 
               FUNCTION SQRT(WS-MU-EARTH * 
               ((2.0 / WS-PERIAPSIS-RADIUS) - 
               (1.0 / WS-SEMI-MAJOR-AXIS))).
           
           COMPUTE WS-APOAPSIS-VELOCITY = 
               FUNCTION SQRT(WS-MU-EARTH * 
               ((2.0 / WS-APOAPSIS-RADIUS) - 
               (1.0 / WS-SEMI-MAJOR-AXIS))).
           
           DISPLAY "  PERIAPSIS VELOCITY: " WS-PERIAPSIS-VELOCITY 
               " KM/S".
           DISPLAY "  APOAPSIS VELOCITY: " WS-APOAPSIS-VELOCITY 
               " KM/S".
           DISPLAY "DERIVED PARAMETERS COMPLETE.".
           DISPLAY " ".
       3000-EXIT.
           EXIT.
       
      *****************************************************************
      * SOLVE KEPLER'S EQUATION USING NEWTON-RAPHSON                  *
      *****************************************************************
       4000-SOLVE-KEPLER-EQUATION.
           DISPLAY "SOLVING KEPLER'S EQUATION...".
           
      *    COMPUTE ECCENTRIC ANOMALY FROM TRUE ANOMALY
           COMPUTE WS-TEMP1 = 
               FUNCTION ATAN(
                   FUNCTION SQRT(1.0 - (WS-ECCENTRICITY ** 2)) *
                   FUNCTION SIN(WS-TEMP2) /
                   (WS-ECCENTRICITY + FUNCTION COS(WS-TEMP2))).
           COMPUTE WS-TEMP1 = WS-TEMP1 * 2.0.
           
      *    COMPUTE MEAN ANOMALY
           COMPUTE WS-MEAN-ANOMALY = WS-TEMP1 - 
               (WS-ECCENTRICITY * FUNCTION SIN(WS-TEMP1)).
           
      *    SOLVE KEPLER'S EQUATION ITERATIVELY
           MOVE WS-MEAN-ANOMALY TO WS-ECCENTRIC-ANOMALY.
           MOVE 0 TO WS-KEPLER-ITERATIONS.
           MOVE 0 TO WS-CONVERGENCE-FLAG.
           
           PERFORM UNTIL WS-CONVERGENCE-FLAG = 1 OR
                         WS-KEPLER-ITERATIONS >= WS-KEPLER-MAX-ITER
               ADD 1 TO WS-KEPLER-ITERATIONS
               
               COMPUTE WS-TEMP1 = WS-ECCENTRIC-ANOMALY -
                   (WS-ECCENTRICITY * 
                   FUNCTION SIN(WS-ECCENTRIC-ANOMALY)) -
                   WS-MEAN-ANOMALY
               
               COMPUTE WS-TEMP2 = 1.0 - 
                   (WS-ECCENTRICITY * 
                   FUNCTION COS(WS-ECCENTRIC-ANOMALY))
               
               COMPUTE WS-TEMP3 = WS-ECCENTRIC-ANOMALY - 
                   (WS-TEMP1 / WS-TEMP2)
               
               IF FUNCTION ABS(WS-TEMP3 - WS-ECCENTRIC-ANOMALY) < 
                   WS-TOLERANCE
                   MOVE 1 TO WS-CONVERGENCE-FLAG
               END-IF
               
               MOVE WS-TEMP3 TO WS-ECCENTRIC-ANOMALY
           END-PERFORM.
           
           COMPUTE WS-ECCENTRIC-ANOMALY = 
               WS-ECCENTRIC-ANOMALY * WS-RAD-TO-DEG.
           COMPUTE WS-MEAN-ANOMALY = WS-MEAN-ANOMALY * WS-RAD-TO-DEG.
           
           DISPLAY "  MEAN ANOMALY: " WS-MEAN-ANOMALY " DEG".
           DISPLAY "  ECCENTRIC ANOMALY: " WS-ECCENTRIC-ANOMALY " DEG".
           DISPLAY "  ITERATIONS: " WS-KEPLER-ITERATIONS.
           DISPLAY "KEPLER'S EQUATION SOLVED.".
           DISPLAY " ".
       4000-EXIT.
           EXIT.
       
      *****************************************************************
      * COMPUTE POSITION AND VELOCITY STATE VECTORS                   *
      *****************************************************************
       5000-COMPUTE-STATE-VECTORS.
           DISPLAY "COMPUTING STATE VECTORS...".
           
      *    CONVERT ECCENTRIC ANOMALY TO RADIANS
           COMPUTE WS-TEMP1 = WS-ECCENTRIC-ANOMALY * WS-DEG-TO-RAD.
           
      *    COMPUTE POSITION IN ORBITAL PLANE
           COMPUTE WS-POSITION-X = WS-SEMI-MAJOR-AXIS * 
               (FUNCTION COS(WS-TEMP1) - WS-ECCENTRICITY).
           COMPUTE WS-POSITION-Y = WS-SEMI-MAJOR-AXIS * 
               FUNCTION SQRT(1.0 - (WS-ECCENTRICITY ** 2)) *
               FUNCTION SIN(WS-TEMP1).
           MOVE 0.0 TO WS-POSITION-Z.
           
      *    COMPUTE RADIUS MAGNITUDE
           COMPUTE WS-RADIUS-MAGNITUDE = 
               FUNCTION SQRT((WS-POSITION-X ** 2) + 
               (WS-POSITION-Y ** 2) + (WS-POSITION-Z ** 2)).
           
      *    COMPUTE VELOCITY IN ORBITAL PLANE
           COMPUTE WS-TEMP2 = FUNCTION SQRT(WS-MU-EARTH * 
               WS-SEMI-MAJOR-AXIS) / WS-RADIUS-MAGNITUDE.
           COMPUTE WS-VELOCITY-X = 
               0.0 - (WS-TEMP2 * FUNCTION SIN(WS-TEMP1)).
           COMPUTE WS-VELOCITY-Y = WS-TEMP2 * 
               FUNCTION SQRT(1.0 - (WS-ECCENTRICITY ** 2)) *
               FUNCTION COS(WS-TEMP1).
           MOVE 0.0 TO WS-VELOCITY-Z.
           
      *    COMPUTE VELOCITY MAGNITUDE
           COMPUTE WS-VELOCITY-MAGNITUDE = 
               FUNCTION SQRT((WS-VELOCITY-X ** 2) + 
               (WS-VELOCITY-Y ** 2) + (WS-VELOCITY-Z ** 2)).
           
           DISPLAY "  POSITION: (" WS-POSITION-X ", " 
               WS-POSITION-Y ", " WS-POSITION-Z ") KM".
           DISPLAY "  VELOCITY: (" WS-VELOCITY-X ", " 
               WS-VELOCITY-Y ", " WS-VELOCITY-Z ") KM/S".
           DISPLAY "  RADIUS: " WS-RADIUS-MAGNITUDE " KM".
           DISPLAY "  VELOCITY MAG: " WS-VELOCITY-MAGNITUDE " KM/S".
           DISPLAY "STATE VECTORS COMPUTED.".
           DISPLAY " ".
       5000-EXIT.
           EXIT.
       
      *****************************************************************
      * COMPUTE HOHMANN TRANSFER ORBIT                                *
      *****************************************************************
       6000-COMPUTE-TRANSFER-ORBIT.
           DISPLAY "COMPUTING HOHMANN TRANSFER ORBIT...".
           
      *    COMPUTE TRANSFER ORBIT SEMI-MAJOR AXIS
           COMPUTE WS-TRANSFER-SMA = 
               (WS-INITIAL-ORBIT-R + WS-TARGET-ORBIT-R) / 2.0.
           
           DISPLAY "  INITIAL ORBIT RADIUS: " WS-INITIAL-ORBIT-R " KM".
           DISPLAY "  TARGET ORBIT RADIUS: " WS-TARGET-ORBIT-R " KM".
           DISPLAY "  TRANSFER ORBIT SMA: " WS-TRANSFER-SMA " KM".
           
      *    COMPUTE INITIAL CIRCULAR VELOCITY
           COMPUTE WS-TEMP1 = WS-MU-EARTH / WS-INITIAL-ORBIT-R.
           PERFORM 2100-COMPUTE-SQUARE-ROOT THRU 2100-EXIT.
           MOVE WS-SQRT-VALUE TO WS-TEMP2.
           
      *    COMPUTE VELOCITY AT PERIAPSIS OF TRANSFER ORBIT
           COMPUTE WS-TEMP1 = WS-MU-EARTH * 
               ((2.0 / WS-INITIAL-ORBIT-R) - (1.0 / WS-TRANSFER-SMA)).
           PERFORM 2100-COMPUTE-SQUARE-ROOT THRU 2100-EXIT.
           MOVE WS-SQRT-VALUE TO WS-TEMP3.
           
      *    COMPUTE FIRST DELTA-V
           COMPUTE WS-DELTA-V1 = WS-TEMP3 - WS-TEMP2.
           
           DISPLAY "  DELTA-V1: " WS-DELTA-V1 " KM/S".
           
      *    COMPUTE TARGET CIRCULAR VELOCITY
           COMPUTE WS-TEMP1 = WS-MU-EARTH / WS-TARGET-ORBIT-R.
           PERFORM 2100-COMPUTE-SQUARE-ROOT THRU 2100-EXIT.
           MOVE WS-SQRT-VALUE TO WS-TEMP2.
           
      *    COMPUTE VELOCITY AT APOAPSIS OF TRANSFER ORBIT
           COMPUTE WS-TEMP1 = WS-MU-EARTH * 
               ((2.0 / WS-TARGET-ORBIT-R) - (1.0 / WS-TRANSFER-SMA)).
           PERFORM 2100-COMPUTE-SQUARE-ROOT THRU 2100-EXIT.
           MOVE WS-SQRT-VALUE TO WS-TEMP3.
           
      *    COMPUTE SECOND DELTA-V
           COMPUTE WS-DELTA-V2 = WS-TEMP2 - WS-TEMP3.
           
           DISPLAY "  DELTA-V2: " WS-DELTA-V2 " KM/S".
           
      *    COMPUTE TOTAL DELTA-V
           COMPUTE WS-TOTAL-DELTA-V = WS-DELTA-V1 + WS-DELTA-V2.
           
           DISPLAY "  TOTAL DELTA-V: " WS-TOTAL-DELTA-V " KM/S".
           
      *    COMPUTE TRANSFER TIME
           COMPUTE WS-TEMP1 = (WS-TRANSFER-SMA ** 3) / WS-MU-EARTH.
           PERFORM 2100-COMPUTE-SQUARE-ROOT THRU 2100-EXIT.
           COMPUTE WS-TRANSFER-TIME = WS-PI * WS-SQRT-VALUE.
           
           MOVE WS-TRANSFER-TIME TO DISP-NUMERIC-3.
           DISPLAY "  TRANSFER TIME: " DISP-NUMERIC-3 " SECONDS".
           COMPUTE WS-TEMP1 = WS-TRANSFER-TIME / 3600.0.
           MOVE WS-TEMP1 TO DISP-NUMERIC-2.
           DISPLAY "                (" DISP-NUMERIC-2 " HOURS)".
           
      *    COMPUTE PHASE ANGLE
           COMPUTE WS-TEMP1 = WS-MU-EARTH / (WS-TARGET-ORBIT-R ** 3).
           PERFORM 2100-COMPUTE-SQUARE-ROOT THRU 2100-EXIT.
           COMPUTE WS-PHASE-ANGLE = WS-PI - 
               (WS-SQRT-VALUE * WS-TRANSFER-TIME).
           COMPUTE WS-PHASE-ANGLE = WS-PHASE-ANGLE * WS-RAD-TO-DEG.
           
           DISPLAY "  REQUIRED PHASE ANGLE: " WS-PHASE-ANGLE " DEG".
           DISPLAY "TRANSFER ORBIT COMPUTED.".
           DISPLAY " ".
       6000-EXIT.
           EXIT.
       
      *****************************************************************
      * ANALYZE ORBITAL PERTURBATIONS                                 *
      *****************************************************************
       7000-ANALYZE-PERTURBATIONS.
           DISPLAY "ANALYZING ORBITAL PERTURBATIONS...".
           
      *    COMPUTE J2 PERTURBATION EFFECTS
      *    NODAL PRECESSION RATE
           COMPUTE WS-TEMP1 = 1.0 - (WS-ECCENTRICITY ** 2).
           COMPUTE WS-TEMP1 = FUNCTION SQRT(WS-TEMP1).
           COMPUTE WS-TEMP2 = WS-INCLINATION * WS-DEG-TO-RAD.
           COMPUTE WS-TEMP3 = FUNCTION COS(WS-TEMP2).
           
           COMPUTE WS-NODAL-PRECESSION = 
               -1.5 * WS-J2-COEFFICIENT * WS-MEAN-MOTION *
               ((WS-EARTH-RADIUS / WS-SEMI-MAJOR-AXIS) ** 2) *
               WS-TEMP3 / (WS-TEMP1 ** 2).
           
           COMPUTE WS-NODAL-PRECESSION = 
               WS-NODAL-PRECESSION * WS-RAD-TO-DEG.
           
           DISPLAY "  J2 NODAL PRECESSION: " WS-NODAL-PRECESSION 
               " DEG/S".
           
      *    APSIDAL PRECESSION RATE
           COMPUTE WS-TEMP4 = 5.0 * (WS-TEMP3 ** 2) - 1.0.
           
           COMPUTE WS-APSIDAL-PRECESSION = 
               0.75 * WS-J2-COEFFICIENT * WS-MEAN-MOTION *
               ((WS-EARTH-RADIUS / WS-SEMI-MAJOR-AXIS) ** 2) *
               WS-TEMP4 / (WS-TEMP1 ** 2).
           
           COMPUTE WS-APSIDAL-PRECESSION = 
               WS-APSIDAL-PRECESSION * WS-RAD-TO-DEG.
           
           DISPLAY "  J2 APSIDAL PRECESSION: " WS-APSIDAL-PRECESSION 
               " DEG/S".
           
      *    COMPUTE THIRD-BODY PERTURBATIONS
      *    SOLAR PERTURBATION (SIMPLIFIED)
           COMPUTE WS-TEMP1 = WS-MU-SUN / (WS-MOON-DISTANCE ** 3).
           COMPUTE WS-SOLAR-PERTURBATION = 
               WS-TEMP1 * WS-RADIUS-MAGNITUDE.
           
           DISPLAY "  SOLAR PERTURBATION: " WS-SOLAR-PERTURBATION 
               " KM/S^2".
           
      *    LUNAR PERTURBATION (SIMPLIFIED)
           COMPUTE WS-TEMP1 = WS-MU-MOON / (WS-MOON-DISTANCE ** 3).
           COMPUTE WS-LUNAR-PERTURBATION = 
               WS-TEMP1 * WS-RADIUS-MAGNITUDE.
           
           DISPLAY "  LUNAR PERTURBATION: " WS-LUNAR-PERTURBATION 
               " KM/S^2".
           
           DISPLAY "PERTURBATION ANALYSIS COMPLETE.".
           DISPLAY " ".
       7000-EXIT.
           EXIT.
       
      *****************************************************************
      * PROPAGATE ORBIT OVER TIME                                     *
      *****************************************************************
       8000-PROPAGATE-ORBIT.
           DISPLAY "PROPAGATING ORBIT OVER TIME...".
           DISPLAY "  SIMULATION DURATION: " WS-SIMULATION-DURATION 
               " SECONDS".
           COMPUTE WS-TEMP1 = WS-SIMULATION-DURATION / 3600.0.
           MOVE WS-TEMP1 TO DISP-NUMERIC-2.
           DISPLAY "                      (" DISP-NUMERIC-2 
               " HOURS)".
           DISPLAY "  TIME STEP: " WS-TIME-STEP " SECONDS".
           
           MOVE 0.0 TO WS-CURRENT-TIME.
           MOVE 0 TO WS-ITERATION-COUNTER.
           
           PERFORM 8100-PROPAGATE-TIME-STEP
               UNTIL WS-CURRENT-TIME >= WS-SIMULATION-DURATION.
           
           DISPLAY "  TOTAL ITERATIONS: " WS-ITERATION-COUNTER.
           DISPLAY "  FINAL TIME: " WS-CURRENT-TIME " SECONDS".
           COMPUTE WS-TEMP1 = WS-CURRENT-TIME / 3600.0.
           MOVE WS-TEMP1 TO DISP-NUMERIC-2.
           DISPLAY "             (" DISP-NUMERIC-2 " HOURS)".
           DISPLAY "  FINAL MEAN ANOMALY: " WS-MEAN-ANOMALY " DEG".
           DISPLAY "  FINAL RAAN: " WS-RAAN " DEG".
           DISPLAY "  FINAL ARG PERIAPSIS: " WS-ARG-PERIAPSIS " DEG".
           DISPLAY "ORBIT PROPAGATION COMPLETE.".
           DISPLAY " ".
       8000-EXIT.
           EXIT.
       
      *****************************************************************
      * PROPAGATE SINGLE TIME STEP                                    *
      *****************************************************************
       8100-PROPAGATE-TIME-STEP.
           ADD 1 TO WS-ITERATION-COUNTER.
           
      *    UPDATE MEAN ANOMALY
           COMPUTE WS-MEAN-ANOMALY = WS-MEAN-ANOMALY + 
               (WS-MEAN-MOTION * WS-TIME-STEP * WS-RAD-TO-DEG).
           
      *    NORMALIZE MEAN ANOMALY
           PERFORM 8110-NORMALIZE-MEAN-ANOMALY
               UNTIL WS-MEAN-ANOMALY < 360.0.
           
      *    APPLY PERTURBATIONS
           COMPUTE WS-RAAN = WS-RAAN + 
               (WS-NODAL-PRECESSION * WS-TIME-STEP).
           COMPUTE WS-ARG-PERIAPSIS = WS-ARG-PERIAPSIS + 
               (WS-APSIDAL-PRECESSION * WS-TIME-STEP).
           
      *    NORMALIZE ANGLES
           PERFORM 8120-NORMALIZE-RAAN
               UNTIL WS-RAAN < 360.0.
           PERFORM 8130-NORMALIZE-ARG-PERIAPSIS
               UNTIL WS-ARG-PERIAPSIS < 360.0.
           
      *    SOLVE KEPLER'S EQUATION FOR NEW TIME STEP
           COMPUTE WS-TEMP1 = WS-MEAN-ANOMALY * WS-DEG-TO-RAD.
           MOVE WS-TEMP1 TO WS-ECCENTRIC-ANOMALY.
           MOVE 0 TO WS-KEPLER-ITERATIONS.
           
           PERFORM 8140-SOLVE-KEPLER-STEP
               UNTIL WS-KEPLER-ITERATIONS >= 10.
           
      *    UPDATE STATE VECTORS
           COMPUTE WS-POSITION-X = WS-SEMI-MAJOR-AXIS * 
               (FUNCTION COS(WS-ECCENTRIC-ANOMALY) - 
               WS-ECCENTRICITY).
           COMPUTE WS-POSITION-Y = WS-SEMI-MAJOR-AXIS * 
               FUNCTION SQRT(1.0 - (WS-ECCENTRICITY ** 2)) *
               FUNCTION SIN(WS-ECCENTRIC-ANOMALY).
           
           COMPUTE WS-RADIUS-MAGNITUDE = 
               FUNCTION SQRT((WS-POSITION-X ** 2) + 
               (WS-POSITION-Y ** 2)).
           
           COMPUTE WS-CURRENT-TIME = 
               WS-CURRENT-TIME + WS-TIME-STEP.
       
      *****************************************************************
      * NORMALIZE MEAN ANOMALY HELPER                                 *
      *****************************************************************
       8110-NORMALIZE-MEAN-ANOMALY.
           COMPUTE WS-MEAN-ANOMALY = WS-MEAN-ANOMALY - 360.0.
       
      *****************************************************************
      * NORMALIZE RAAN HELPER                                         *
      *****************************************************************
       8120-NORMALIZE-RAAN.
           COMPUTE WS-RAAN = WS-RAAN - 360.0.
       
      *****************************************************************
      * NORMALIZE ARGUMENT OF PERIAPSIS HELPER                        *
      *****************************************************************
       8130-NORMALIZE-ARG-PERIAPSIS.
           COMPUTE WS-ARG-PERIAPSIS = WS-ARG-PERIAPSIS - 360.0.
       
      *****************************************************************
      * KEPLER EQUATION SOLVER ITERATION STEP                         *
      *****************************************************************
       8140-SOLVE-KEPLER-STEP.
           ADD 1 TO WS-KEPLER-ITERATIONS.
           COMPUTE WS-TEMP2 = WS-ECCENTRIC-ANOMALY -
               (WS-ECCENTRICITY * 
               FUNCTION SIN(WS-ECCENTRIC-ANOMALY)) - WS-TEMP1.
           COMPUTE WS-TEMP3 = 1.0 - (WS-ECCENTRICITY * 
               FUNCTION COS(WS-ECCENTRIC-ANOMALY)).
           COMPUTE WS-ECCENTRIC-ANOMALY = 
               WS-ECCENTRIC-ANOMALY - (WS-TEMP2 / WS-TEMP3).
       
      *****************************************************************
      * GENERATE COMPREHENSIVE REPORT                                 *
      *****************************************************************
       9000-GENERATE-REPORT.
           DISPLAY "GENERATING COMPREHENSIVE REPORT...".
           
           OPEN OUTPUT ORBIT-REPORT-FILE.
           
      *    WRITE REPORT HEADER
           MOVE 1 TO WS-PAGE-COUNTER.
           MOVE WS-PAGE-COUNTER TO RPT-PAGE-NUM.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-HEADER-1.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-BLANK-LINE.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-HEADER-2.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-SEPARATOR-LINE.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-BLANK-LINE.
           
      *    ORBITAL ELEMENTS SECTION
           MOVE "KEPLERIAN ORBITAL ELEMENTS" TO RPT-SECTION-NAME.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-SECTION-HEADER.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-BLANK-LINE.
           
           MOVE "SEMI-MAJOR AXIS" TO RPT-PARAMETER-NAME.
           MOVE WS-SEMI-MAJOR-AXIS TO DISP-NUMERIC-1.
           MOVE DISP-NUMERIC-1 TO RPT-PARAMETER-VALUE.
           MOVE "KM" TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           MOVE "ECCENTRICITY" TO RPT-PARAMETER-NAME.
           MOVE WS-ECCENTRICITY TO DISP-NUMERIC-2.
           MOVE DISP-NUMERIC-2 TO RPT-PARAMETER-VALUE.
           MOVE " " TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           MOVE "INCLINATION" TO RPT-PARAMETER-NAME.
           MOVE WS-INCLINATION TO DISP-NUMERIC-2.
           MOVE DISP-NUMERIC-2 TO RPT-PARAMETER-VALUE.
           MOVE "DEGREES" TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           MOVE "RAAN" TO RPT-PARAMETER-NAME.
           MOVE WS-RAAN TO DISP-NUMERIC-2.
           MOVE DISP-NUMERIC-2 TO RPT-PARAMETER-VALUE.
           MOVE "DEGREES" TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           MOVE "ARGUMENT OF PERIAPSIS" TO RPT-PARAMETER-NAME.
           MOVE WS-ARG-PERIAPSIS TO DISP-NUMERIC-2.
           MOVE DISP-NUMERIC-2 TO RPT-PARAMETER-VALUE.
           MOVE "DEGREES" TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           WRITE ORBIT-REPORT-RECORD FROM REPORT-BLANK-LINE.
           
      *    DERIVED PARAMETERS SECTION
           MOVE "DERIVED ORBITAL PARAMETERS" TO RPT-SECTION-NAME.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-SECTION-HEADER.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-BLANK-LINE.
           
           MOVE "ORBITAL PERIOD" TO RPT-PARAMETER-NAME.
           MOVE WS-ORBITAL-PERIOD TO DISP-NUMERIC-1.
           MOVE DISP-NUMERIC-1 TO RPT-PARAMETER-VALUE.
           MOVE "SECONDS" TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           MOVE "PERIAPSIS RADIUS" TO RPT-PARAMETER-NAME.
           MOVE WS-PERIAPSIS-RADIUS TO DISP-NUMERIC-1.
           MOVE DISP-NUMERIC-1 TO RPT-PARAMETER-VALUE.
           MOVE "KM" TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           MOVE "APOAPSIS RADIUS" TO RPT-PARAMETER-NAME.
           MOVE WS-APOAPSIS-RADIUS TO DISP-NUMERIC-1.
           MOVE DISP-NUMERIC-1 TO RPT-PARAMETER-VALUE.
           MOVE "KM" TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           MOVE "PERIAPSIS VELOCITY" TO RPT-PARAMETER-NAME.
           MOVE WS-PERIAPSIS-VELOCITY TO DISP-NUMERIC-2.
           MOVE DISP-NUMERIC-2 TO RPT-PARAMETER-VALUE.
           MOVE "KM/S" TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           MOVE "APOAPSIS VELOCITY" TO RPT-PARAMETER-NAME.
           MOVE WS-APOAPSIS-VELOCITY TO DISP-NUMERIC-2.
           MOVE DISP-NUMERIC-2 TO RPT-PARAMETER-VALUE.
           MOVE "KM/S" TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           WRITE ORBIT-REPORT-RECORD FROM REPORT-BLANK-LINE.
           
      *    TRANSFER ORBIT SECTION
           MOVE "HOHMANN TRANSFER ORBIT ANALYSIS" TO RPT-SECTION-NAME.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-SECTION-HEADER.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-BLANK-LINE.
           
           MOVE "DELTA-V BURN 1" TO RPT-PARAMETER-NAME.
           MOVE WS-DELTA-V1 TO DISP-NUMERIC-2.
           MOVE DISP-NUMERIC-2 TO RPT-PARAMETER-VALUE.
           MOVE "KM/S" TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           MOVE "DELTA-V BURN 2" TO RPT-PARAMETER-NAME.
           MOVE WS-DELTA-V2 TO DISP-NUMERIC-2.
           MOVE DISP-NUMERIC-2 TO RPT-PARAMETER-VALUE.
           MOVE "KM/S" TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           MOVE "TOTAL DELTA-V" TO RPT-PARAMETER-NAME.
           MOVE WS-TOTAL-DELTA-V TO DISP-NUMERIC-2.
           MOVE DISP-NUMERIC-2 TO RPT-PARAMETER-VALUE.
           MOVE "KM/S" TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           MOVE "TRANSFER TIME" TO RPT-PARAMETER-NAME.
           MOVE WS-TRANSFER-TIME TO DISP-NUMERIC-1.
           MOVE DISP-NUMERIC-1 TO RPT-PARAMETER-VALUE.
           MOVE "SECONDS" TO RPT-PARAMETER-UNIT.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-DETAIL-LINE.
           
           WRITE ORBIT-REPORT-RECORD FROM REPORT-BLANK-LINE.
           WRITE ORBIT-REPORT-RECORD FROM REPORT-SEPARATOR-LINE.
           
           CLOSE ORBIT-REPORT-FILE.
           
           DISPLAY "REPORT GENERATED: ORBITRPT.TXT".
           DISPLAY " ".
       9000-EXIT.
           EXIT.
       
      *****************************************************************
      * TERMINATE PROGRAM                                             *
      *****************************************************************
       9999-TERMINATE-PROGRAM.
           DISPLAY "========================================".
           DISPLAY "ORBITAL MECHANICS COMPUTATION COMPLETE".
           DISPLAY "(c) 2025 by moshix. All right reserved".
           DISPLAY "========================================".
           DISPLAY " ".
           DISPLAY "SUMMARY OF COMPUTATIONS:".
           DISPLAY "  - ORBITAL ELEMENTS CALCULATED".
           DISPLAY "  - KEPLER'S EQUATION SOLVED".
           DISPLAY "  - STATE VECTORS COMPUTED".
           DISPLAY "  - TRANSFER ORBIT ANALYZED".
           DISPLAY "  - PERTURBATIONS EVALUATED".
           COMPUTE WS-TEMP1 = WS-ITERATION-COUNTER.
           MOVE WS-TEMP1 TO DISP-NUMERIC-3.
           DISPLAY "  - " DISP-NUMERIC-3 " ORBIT PROPAGATIONS".
           DISPLAY " ".
           DISPLAY "ALL CALCULATIONS COMPLETE.".

