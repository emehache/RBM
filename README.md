# RBM

## Simulated data

#### 1 Simulate data.  `Rscript simulations.R n_rep=50 N=1e5`
* (default values) generates trajectories of length `N` and step `h = {0.001, 0.002, 0.003}`. Simulated data is stored in `objetos_r/trayectorias/trayectoria_%N_%h.rds`. In every `.rds` file there are stored `n_rep` trajectories.

#### 2 Cálculo de distancias de Hausdorff

* `calculo_distancias_H_enteras_comparable.R` genera `objetos_r/distH/dist_entera_%s-%s.rds'`
* `calculo_distancias_H_on_off_replicas2.R` genera `objetos_r/distH/cadena_%s-%s.rds`

* Esto no lo usamos en este proyecto `calculo_distancias_H_replicas2.rds` genera `objetos_r/distH/cadena_entera_%s-%s.rds`

### 3 Cierre r-convexo de las trayectorias

### 4 Cálculo de distancia en medida

* `calculo_distancias_medida_enteras_comparable.R` genera `objetos_r/dist_medida/dist_entera-%s.rds'`
* `calculo_distancias_medida_on_off.R` genera `objetos_r/dist_medida/cadena_on_off_%s.rds`



## Scripts con funciones

* `rbmd.cpp` Simula el Reflected Brownian Motion with Drift en la elipse agujereada.
* `dist_medida.R` Estima la distancia en medida entre el cierre r-convexo de una trayectoria del browniano y el conjunto.
* `distPaQ.cpp` Aproxima la distancia de Hausdorff entre una trayectoria del browniano y el conjunto.
* `on_off.R` Aplica el modelo de observación on/off. Recibe como input una trayectoria, y devuelve la trayectoria indicando los tiempos en que se observa.
* `graficar.R` Grafica una trayectoria.



## Objetos

* En `objetos_r/trayectorias/trayectoria_%s, %s.rds` están las trayectorias con `N = %s` y `h = %s`. En cada objeto hay `n_rep` trayectorias guardadas.

[//]: # * `objetos_r/sim_replicas.rds` tiene lo anterior, todo junto en un solo objeto.
[//]: # * `objetos_r/sim.rds` tiene lo anterior sin repeticiones.

* `objetos_r/cr_tray_sim.rds` tiene el cierre $r$-convexo generado por el script `cierre_r_convexo.R`

