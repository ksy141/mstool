#!/bin/bash

rm -rf output.dtr

garden with -m desmond-configure-for/1.36c7/bin desmond_configure_for.py \
	input.dms \
    --ensemble npt    \
    --dt 0.01         \
    --last-time 0.0   \
    --temp      310   \
    --pressure  1.0   \
    --barostat  MTK   \
    --output-dms output.dms \
    --cfg force.nonbonded.type=cutoff \
    --cfg force.nonbonded.far=none \
    --cfg force.near.average_dispersion=0 \
    --cfg force.nonbonded.r_cut=11 \
    --cfg force.nonbonded.near.average_dispersion=0 \
    --output-ark sim.cfg \
    --add-plugin energy_groups \
    --add-plugin eneseq \
    --add-plugin trajectory \
    --add-plugin randomize_velocities \
    --cfg mdsim.plugin.trajectory.interval=10  \
    --cfg mdsim.plugin.trajectory.name=output.dtr \

garden with -m gdesmond/1.9.6.1-01c7/bin gdesmond --include sim.cfg
cat egrp

