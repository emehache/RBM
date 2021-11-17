# RBM

## Simulated data

#### 1 Simulate data.  `Rscript simulations.R n_rep=50 N=1e5`
* (default values) generates trajectories of length `N` and step `h = {0.001, 0.002, 0.003}`. Simulated data is stored in `objetos_r/trayectorias/trayectoria_%N_%h.rds`. In every `.rds` file there are stored `n_rep` trajectories.

#### 2 Compute Hausdorff distances

* `compute_Hausdorff_dist_entire.R` generates `r_objects/hausdorff_distances/entire_model-%s-%s.rds'`
* `compute_Hausdorff_dist_onoff.R` generates `r_objects/hausdorff_distances/onoff_model-%s-%s.rds`

### 4 Compute measure distances

* `compute_measure_dist_entire.R` generates `r_objects/measure_distances/entire_model-%s.rds'`
* `compute_measure_dist_onoff.R` generates `r_objects/measure_distances/onoff_model-%s.rds'`


## Auxiliar functions

* `rbmd.cpp` Simulates a Reflected Brownian Motion with Drift in the holled ellipse.
* `dist_medida.R` Computes distance measure between the r-convex hull of the trajectory and the set.
* `distPaQ.cpp` Computes an aproximation of the Hausdorff distance between the trajectory and the set.
* `on_off.R` Applies on/off model to a trajectory.
