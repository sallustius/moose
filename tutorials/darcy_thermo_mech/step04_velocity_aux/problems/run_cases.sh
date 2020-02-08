#!/bin/bash

../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=true -ksp_gmres_restart 1000 -snes_rtol 1e-12 -snes_mf_operator 0 2>&1 | tee no-scaling-yes-preset-newton-amg.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=true GlobalParams/preset=true -ksp_gmres_restart 1000 -snes_rtol 1e-12 -snes_mf_operator 0 2>&1 | tee auto-scaling-yes-preset-newton-amg.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=true Variables/pressure/scaling=1e12 -ksp_gmres_restart 1000 -snes_rtol 1e-12 -snes_mf_operator 0 2>&1 | tee over-scaling-yes-preset-newton-amg.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=false -ksp_gmres_restart 1000 -snes_rtol 1e-12 -snes_mf_operator 0 2>&1 | tee no-scaling-no-preset-newton-amg.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=true GlobalParams/preset=false -ksp_gmres_restart 1000 -snes_rtol 1e-12 -snes_mf_operator 0 2>&1 | tee auto-scaling-no-preset-newton-amg.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=false Variables/pressure/scaling=1e12 -ksp_gmres_restart 1000 -snes_rtol 1e-12 -snes_mf_operator 0 2>&1 | tee over-scaling-no-preset-newton-amg.log

../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=true -ksp_gmres_restart 1000 -snes_rtol 1e-12 -snes_mf_operator 1 2>&1 | tee no-scaling-yes-preset-pjfnk-amg.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=true GlobalParams/preset=true -ksp_gmres_restart 1000 -snes_rtol 1e-12 -snes_mf_operator 1 2>&1 | tee auto-scaling-yes-preset-pjfnk-amg.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=true Variables/pressure/scaling=1e12 -ksp_gmres_restart 1000 -snes_rtol 1e-12 -snes_mf_operator 1 2>&1 | tee over-scaling-yes-preset-pjfnk-amg.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=false -ksp_gmres_restart 1000 -snes_rtol 1e-12 -snes_mf_operator 1 2>&1 | tee no-scaling-no-preset-pjfnk-amg.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=true GlobalParams/preset=false -ksp_gmres_restart 1000 -snes_rtol 1e-12 -snes_mf_operator 1 2>&1 | tee auto-scaling-no-preset-pjfnk-amg.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=false Variables/pressure/scaling=1e12 -ksp_gmres_restart 1000 -snes_rtol 1e-12 -snes_mf_operator 1 2>&1 | tee over-scaling-no-preset-pjfnk-amg.log

../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=true -ksp_gmres_restart 1000 -snes_rtol 1e-12 -pc_type ilu -snes_mf_operator 0 2>&1 | tee no-scaling-yes-preset-newton-ilu.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=true GlobalParams/preset=true -ksp_gmres_restart 1000 -snes_rtol 1e-12 -pc_type ilu -snes_mf_operator 0 2>&1 | tee auto-scaling-yes-preset-newton-ilu.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=true Variables/pressure/scaling=1e12 -ksp_gmres_restart 1000 -snes_rtol 1e-12 -pc_type ilu -snes_mf_operator 0 2>&1 | tee over-scaling-yes-preset-newton-ilu.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=false -ksp_gmres_restart 1000 -snes_rtol 1e-12 -pc_type ilu -snes_mf_operator 0 2>&1 | tee no-scaling-no-preset-newton-ilu.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=true GlobalParams/preset=false -ksp_gmres_restart 1000 -snes_rtol 1e-12 -pc_type ilu -snes_mf_operator 0 2>&1 | tee auto-scaling-no-preset-newton-ilu.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=false Variables/pressure/scaling=1e12 -ksp_gmres_restart 1000 -snes_rtol 1e-12 -pc_type ilu -snes_mf_operator 0 2>&1 | tee over-scaling-no-preset-newton-ilu.log

../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=true -ksp_gmres_restart 1000 -snes_rtol 1e-12 -pc_type ilu -snes_mf_operator 1 2>&1 | tee no-scaling-yes-preset-pjfnk-ilu.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=true GlobalParams/preset=true -ksp_gmres_restart 1000 -snes_rtol 1e-12 -pc_type ilu -snes_mf_operator 1 2>&1 | tee auto-scaling-yes-preset-pjfnk-ilu.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=true Variables/pressure/scaling=1e12 -ksp_gmres_restart 1000 -snes_rtol 1e-12 -pc_type ilu -snes_mf_operator 1 2>&1 | tee over-scaling-yes-preset-pjfnk-ilu.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=false -ksp_gmres_restart 1000 -snes_rtol 1e-12 -pc_type ilu -snes_mf_operator 1 2>&1 | tee no-scaling-no-preset-pjfnk-ilu.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=true GlobalParams/preset=false -ksp_gmres_restart 1000 -snes_rtol 1e-12 -pc_type ilu -snes_mf_operator 1 2>&1 | tee auto-scaling-no-preset-pjfnk-ilu.log
../darcy_thermo_mech-oprof -i step4.i Executioner/automatic_scaling=false GlobalParams/preset=false Variables/pressure/scaling=1e12 -ksp_gmres_restart 1000 -snes_rtol 1e-12 -pc_type ilu -snes_mf_operator 1 2>&1 | tee over-scaling-no-preset-pjfnk-ilu.log
