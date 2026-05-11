#!/bin/bash
python analyze_traj.py \
  --traj ../o.dump_all.lammpstrj \
  --dt 0.002 \
  --nbeads 8 \
  --temperature 1.0 \
  --gamma 10.0 \
  --b 1.0 \
  --msd-fit-tmin 0 \
  --msd-fit-tmax 200 \

