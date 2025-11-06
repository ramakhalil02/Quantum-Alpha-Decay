# Numerical Modeling of Alpha Decay: Quantum Tunneling Through the Coulomb Barrier

This project contains the computational implementation used to model and investigate the fundamental phenomenon of **quantum tunneling** that governs alpha decay in heavy nuclei.

The primary objective was to numerically verify the **exponential nature** of the half-life's dependence on the alpha particle's emitted energy, which can vary by **more than 20 orders of magnitude** despite minimal changes in kinetic energy.

---

## Table of Contents

* [About The Project](#about-the-project)
* [Key Technical Achievements & Results](#key-technical-achievements--results)
* [Technology Stack](#technology-stack)
* [Getting Started](#getting-started)
* [Full Project Report & Context](#full-project-report--context)
* [Contact](#contact)

---

## About The Project

This simulation successfully solves the Time-Independent Schrödinger Equation (TISE) across a segmented Coulomb potential. The core challenge addressed is demonstrating the strong, non-linear relationship between a minor variance in alpha particle energy and the massive resulting difference in decay half-life, which spans from nanoseconds to billions of years.

The numerical approach models the continuous Coulomb barrier using the **piecewise constant potential approximation** across three discrete regions: inside the nucleus, under the barrier, and in free space.

---

## Key Technical Achievements & Results

The solver calculated the tunneling probability for two contrasting isotopes:

* **Isotopes Modeled:** Polonium-212 ($\text{Po}^{212}$) and Uranium-238 ($\text{U}^{238}$).
* **Result Validation:** Calculated the **Transmission Coefficient ($T$)** and the corresponding half-lives ($t_{1/2}$) for both isotopes.
* **Computational Proof:** The simulation demonstrated that the relatively small difference in alpha energy between $\text{Po}^{212}$ and $\text{U}^{238}$ resulted in the wavefunction amplitude decaying by **over 35 orders of magnitude** across the barrier. This quantitatively confirms the vast difference in observed half-lives.

---

## Technology Stack

| Category | Tools & Libraries | Competency Demonstrated |
| :--- | :--- | :--- |
| **Language** | Python | Efficient development and handling of complex number mathematics. |
| **Numerical** | NumPy, SciPy | Advanced array manipulation and solving systems of linear equations for the TISE. |
| **Visualization** | Matplotlib | Generating high-quality, physics-focused data plots. |

---

## Getting Started

### Execution

To run the simulation and generate the results and visualizations, execute the core solver script:

```bash
python Alpha_Decay.py 
```
---

## Full Project Report

For a complete breakdown of the theoretical derivations and full results, please see the final project report:

[**Full Alpha Decay Project Report (PDF)**](AlphaDecay_report.pdf)

---

## Contact

Rama Khalil

