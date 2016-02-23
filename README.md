# EXO DiPhoton Analysis Code

## Instructions

export SCRAM_ARCH=slc6_amd64_gcc493 (bash)

cmsrel CMSSW_7_6_3_patch2 
cd CMSSW_7_6_3_patch2/src 
cmsenv 
git clone git@github.com:abuccilli/diphoton-analysis 
cd diphoton-analysis 
scram b -j 8
