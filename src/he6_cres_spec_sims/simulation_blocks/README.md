# simulation blocks

This directory contains the core of the Monte Carlo simulation, starting with generating fake betas (according to user-set probability distributions), computing the resultant trajectories, and optionally writing fake .spec(k) files for further analysis.

## Blocks

1. **Physics**:
   - Generates distributions of beta kinematic parameters including position and velocity.

2. **Event Builder**:
   - Selects particles that remain trapped, avoiding collisions with waveguide walls.
   - Applies a minimum angle threshold (θ_0 ≥ θ_min).

3. **Segment Builder**:
   - Computes main track (segment) parameters including:
     - **Durations** of track segments.
     - **Start Frequencies** (average magnetic field)
     - **Powers** and **slopes**.
     - **Axial frequencies**.
     - Assigns start and end times to all tracks, including sidebands.

4. **Band Builder**:
   - Adds sidebands to main tracks and calculates their power contributions.

6. **DMTrack Builder**:
   - Downmixes track frequencies to DAQ bandwidth

7. **DAQ**:
   - Passes tracks through a simulated DAQ system for validation.
   - Generates a time-series chirp with noise and applies FFT.
   - Outputs the result to a `.spec(k)` file for further offline analysis with Katydid.
