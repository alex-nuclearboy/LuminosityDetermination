#!/bin/bash
mkdir ./input
scp wasa:/data7/users/khreptak/OUTPUT/DATA/LUMIN/DATA_ppn_qf_offset.root ./input
scp wasa:/data7/users/khreptak/OUTPUT/MC/LUMIN/MC-ppn_qf-PARIS-x6.root ./input
scp wasa:/data7/users/khreptak/OUTPUT/MC/LUMIN/MC-ppn_qf-CDBONN-x6.root ./input
scp wasa:/data7/users/khreptak/OUTPUT/MC/LUMIN/MC-ppn_qf-CDBONN-x1.root ./input
scp wasa:/data7/users/khreptak/OUTPUT/MC/LUMIN/MC-pd-momcut.root ./input
scp wasa:/data7/users/khreptak/OUTPUT/MC/LUMIN/MC-dnpi_plus_thetacut-x6.root ./input
mkdir ./output
mkdir ./output/plots
mkdir ./output/files
mkdir ./ThetaCD_fit/plots
mkdir ./ThetaCD_fit/files
scp wasa:/data7/users/khreptak/SIMULATION_OUTPUT/PLUTO_OUTPUT/ppn_qf/ProtonVariables/ProtonVariables-PARIS-1.dat ./input
scp wasa:/data7/users/khreptak/SIMULATION_OUTPUT/PLUTO_OUTPUT/ppn_qf/ProtonVariables/ProtonVariables-CDBONN-1.dat ./input
