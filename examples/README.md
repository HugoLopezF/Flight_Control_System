# Examples

These scripts are the canonical runnable entry points for the project. Each one is small and single-purpose so you can focus on a specific workflow (aircraft initialization, frequency analysis, SAS design, autopilot design, etc.).

**Quick start**
1. Run from the repo root so imports resolve correctly.
2. Prefer module execution:

```bash
python -m examples.01_initialize_aircraft
```

Direct execution also works if you prefer:

```bash
python examples/01_initialize_aircraft.py
```

**Configuration**
- Shared defaults live in `examples/common.py`.
- Update `MODEL`, sweep ranges, components, and gains there to affect all examples consistently.

**Scripts**
1. `01_initialize_aircraft.py` - load aircraft + linearized system.
2. `02_frequency_analysis.py` - Bode/Nichols + channel metrics.
3. `03_time_response_aircraft.py` - aircraft time response.
4. `04_compare_components.py` - actuator/sensor comparison.
5. `05_sensitivity_aircraft.py` - stability coefficient sweep.
6. `06_sas_sensitivity.py` - SAS gain sweep + pzmap (+ optional Nichols).
7. `07_sas_build_response.py` - build SAS and plot response.
8. `08_autopilot_sweep_pid.py` - PID sweep for AP + pzmap (+ optional Nichols).
9. `09_autopilot_build_response.py` - build AP and plot response.

**Notes**
- The default aircraft model is set by `MODEL` in `examples/common.py`.
- Most scripts display plots; they do not save figures unless you add that explicitly.
