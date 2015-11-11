# ZFinder

## Setting up a working area

To set up a working area, run the following commands:

```bash
scramv1 project -n CMSSW_5_3_28_ZPhysics CMSSW CMSSW_5_3_28
cd CMSSW_5_3_28_ZPhysics/

git clone -b JPsi --single-branch https://github.com/jturkewitz/ZFinder.git src/
cd src
cmsenv
scram build -j 4
```
Then to run the analysis on for example a test file
cmsRun ZFinder/Event/zfinder_cfg.py
Note that you'll have to change the filename to the file you would like to run over.

The heart of the analysis is in ZFinder/Event/
Specifically ZFinder/Event/src/ZFinderEvent.cc applies the selection criteria to search for Z and JPsi particles decaying leptonically in the event. Also important are ZFinder/Event/src/ZFinder.cc and  ZFinder/Event/src/ZFinderPlotter.cc which together create histograms at various cut levels (e.g. a Z in the event, a Z and a JPsi in the event). 

These histograms are then used by the scripts in RooScripts/scripts/.
```bash
cd Rooscripts/scripts
make
./RooFitCombined.exe dimuon8_2012b_crab.root zjpsi_DoubleMuonFile.root zjpsi_DoubleElectronFile.root output_directory_of_images
```
where the files were produced by running zfinder over the corresponding datasets. For the dimuon8 file which fits the parameters of the fit by using inclusive JPsi dataset:
condor_filelist.perl ZFinder/Event/zfinder_dimuon8_cfg.py Metadata/filelists/MuOniaParked2012B.txt --batch=1 --jobname=jpsiTest401_MuOniaPartial2012B_dimuon8_eta24_ptsub25 --xrootd
This is the syntax done with condor and Minnesota, but running with crab over the MuOniaParked2012B dataset should be equivalent.
The zjpsi_DoubleMuonFile.root and zjpsi_DoubleElectronFile.root come from ZFinder/Event/zfinder_cfg.py run on the DoubleMuon and DoubleElectron datasets respectively.

If running the analysis more than once, it is useful to skim the double electron and double muon datasets, which can be done with:

JPsiFilter/JPsiFilter/python/jpsi_z_to_ee_filter_cfg.py
JPsiFilter/JPsiFilter/python/jpsi_z_to_mumu_filter_cfg.py
And running via crab (or condor with --xrootd at Minnesota). This will filter out most events, keeping events that are able to have both a Z and a JPsi to mumu.
