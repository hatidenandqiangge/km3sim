#!/usr/bin/env zsh

../build/km3sim \
  -p INPUTParametersRun_3500_scatter_WPD3.3_p0.0075_1.0 \
  -d geo_orca_120strings_av23min20mhorizontal_18OMs_alt9mvertical_v1.gdml \
  -i e_E_10_3000.evt \
  -o e_E_10_3000.out.evt
