#!/bin/bash
for obs in spherocity thrust transverse_spherocity transverse_thrust broadening isotropy; do
    python evaluate_iterations.py --config config/dataset_v7_MCEPOS_unfolddata_optimize_1p3M.json    --method multifold --migiter 0 1 --obs nparticle,${obs} --eff-acc
done
python merge_bs.py --input config/merge_unfolddata.txt --output results_finebin_v7_MCEPOS_unfolddata_1d_1p3M_eff_acc_ensemble4/merge.csv 
python plot_iter.py --obs nparticle,spherocity nparticle,thrust nparticle,transverse_spherocity nparticle,transverse_thrust nparticle,broadening nparticle,isotropy --config config/plot_1d_v7_MCEPOS_unfolddata.json 

