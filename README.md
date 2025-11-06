Quantum Tunneling Simulation

⚛️ Project Overview: Computational Quantum Dynamics

This project presents a robust numerical simulation demonstrating the phenomenon of quantum tunneling through the Coulomb barrier, a core mechanism governing alpha decay in heavy atomic nuclei like Polonium-212 and Uranium-238.

This work serves as a practical application of numerical methods to solve the time-independent Schrödinger equation (TISE) across piecewise constant potentials, directly calculating the transmission coefficient and half-life variations for isotopes with slight differences in decay energy.

Key Achievements

Core Physics: Successfully modeled the strong exponential dependence of half-life on emitted alpha particle energy, spanning over 20 orders of magnitude between isotopes.

Numerical Methodology: Implemented a piecewise constant potential approximation to model the complex Coulomb barrier. The TISE was numerically solved across three critical regions (inside the nucleus, under the barrier, and in the free region).

Output Metrics: Calculated the Transmission Coefficient ($T$) and the resulting half-life ($t_{1/2}$) for both ${}^{212}\text{Po}$ and ${}^{238}\text{U}$, demonstrating mastery of quantum mechanics and numerical integration.

Code Verification: Utilized established boundary conditions to ensure continuity of the wave function ($\psi$) and its derivative ($\psi'$) across all segments.

Technology Stack

Language: Python

Libraries: NumPy, SciPy

Advanced matrix manipulation, linear algebra, and complex number handling.

Visualization: Matplotlib

Generation of log-probability density plots to visualize quantum behavior.

Repository Structure

Code: The primary solver script is located in the root directory.

Report: The full theoretical and methodological details are available in the dedicated documentation folder.

Full Project Report

The complete analysis, including the theoretical framework, derivation of boundary conditions, and detailed results are presented in the report. 
