#!/bin/sh

# Apply MJO-filter to daily NCEP dataset
ncl 'model="NCEP"' 'expname="historical"' ymdStrt=19790101 ymdLast=20081231 get_anom.ncl

# Get MJO EOF 1 and 2 in NCEP
ncl 'model="NCEP"' 'expname="historical"' ymdStrt=19790101 ymdLast=20081231 get_EOF.ncl

# Apply MJO-filter to daily outputs in CNRM-CM5
ncl 'model="CNRM-CM5"' 'expname="historical"' ymdStrt=19790101 ymdLast=20041231 get_anom.ncl
ncl 'model="CNRM-CM5"' 'expname="historical"' ymdStrt=19790101 ymdLast=20041231 get_anom_3d.ncl

# Get MJO Phase in CNRM-CM5 bases on EOF of NCEP
ncl 'model="CNRM-CM5"' 'expname="historical"' ymdStrt=19790101 ymdLast=20041231 get_RMM.ncl

# Get MJO-Phase-Composited Fields
ncl 'model="CNRM-CM5"' 'expname="historical"' ymdStrt=19790101 ymdLast=20041231 get_phase_composite.ncl
