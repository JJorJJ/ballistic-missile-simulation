# Ballistic Missile Interception Simulation

**Course Project: DRSP** **Language:** MATLAB  

## Overview
This project simulates the trajectory of a ballistic missile and evaluates different methods for interceptor guidance. It models the missile's flight path under noisy measurement conditions and compares two prediction algorithms—**Polynomial Fitting** and the **Wiener Filter**—to determine the optimal launch time for an Arrow interceptor.

## Key Features
* **Trajectory Simulation:** Generates a 2D/3D ballistic arc with simulated AR(1) measurement noise.
* **Spectral Analysis:** Analyzes noise characteristics using Periodogram, Welch, and Yule-Walker PSD estimation.
* **Descent Detection:** Uses a Matched Filter to automatically detect the start of the missile's descent phase.
* **Prediction Algorithms:**
    * **Polynomial Fit:** Fits a 2nd-degree polynomial to the descent trajectory.
    * **Wiener Filter:** Uses statistical signal processing to predict future trajectory points.
* **Monte Carlo Validation:** Runs 100 trials to statistically compare the accuracy and variance of both prediction methods.
* **Visualization:** Includes 3D plotting of intercept geometry and a time-stepped flight animation.

## How to Run
1.  Open MATLAB.
2.  Navigate to the project folder.
3.  Run the script:
    ```matlab
    missile_interception
    ```

## Results Summary
The simulation compares the "Clean" trajectory against the "Noisy" sensor data.
* **Polynomial Fit:** Generally provides a more stable intercept point but is sensitive to the window size chosen for fitting.
* **Wiener Prediction:** Theoretically optimal for stationary noise but showed higher sensitivity to the non-stationary nature of the ballistic descent in this specific configuration.

## Author
Jordan Jacob