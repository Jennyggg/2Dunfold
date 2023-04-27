#!/bin/bash
for obs in spherocity thrust transverse_spherocity transverse_thrust broadening isotropy; do
    python evaluate_iterations.py --config config/dataset_v7_MCEPOS_genweightCP1_optimize_1p3M.json    --method multifold --migiter 0 1 --obs nparticle,${obs} --step1
done
    python plot_iter.py --obs nparticle,spherocity nparticle,thrust nparticle,transverse_spherocity nparticle,transverse_thrust nparticle,broadening nparticle,isotropy --config config/plot_1d_v7_MCEPOS_genweightCP1.json

