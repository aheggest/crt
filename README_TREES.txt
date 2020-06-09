Last modified: June 9th, 2020
Modified by: Chris.Hilgenberg@colostate.edu

Contents of trees (n-tuples) produced by analyzer module, icaruscode/CRT/CRTSimAnalysis_module.cc 
Note the coordinate system has origin at LAr center (between 2 cold vessels/cryostats) with
  x-> along drift direction (east to west)
  y-> vertical (bottom to top)
  z-> along BNB direction (south to north)

Region codes used below:
  - LAr regions
    + 4  = Cryostat 0, TPC 0
    + 5  = Cryostat 0, TPC 1
    + 10 = Cryostat 0, Inactive
    + 7  = Cryostat 1, TPC 0
    + 8  = Cryostat 1, TPC 1
    + 12 = Cryostat 1, Inactive
  - CRT regions
    + 30 = Top, Roof
    + 31 = Top, Rim-West
    + 32 = Top, Rim-East
    + 33 = Top, Rim-South
    + 34 = Top, Rim-North
    + 40 = Side, West-South
    + 41 = Side, West-Center
    + 42 = Side, West-North
    + 43 = Side, East-South
    + 44 = Side, East-Center
    + 45 = Side, East-North
    + 46 = Side, South
    + 47 = Side, North
    + 50 = Bottom

// trees with branch discription and data type
GenTree -> generator level info with 1 entry per particle
  - event:     event ID (int)
  - run:       run ID   (int)
  - subRun:    subRun ID (int)
  - nGen:      number of particles generated (int)
  - trackID:   generator trackID (int)
  - pdg:       particle PDG code (e.g. 13 = muon, 11 = electron, 22 = photon)
  - startXYZT: particle start position [cm] in world coords, start time [ns] (vector<double> = {x,y,z,t})
  - endXYZT:   particle end position [cm] in world coords, end time [ns] (vector<double> = {x,y,z,t})
  - startPE:   particle start momentum components [GeV/c] in world coords, start total energy [GeV] (vector<double> = {Px,Py,Pz,E})
  - endPE:     particle end momentum components [GeV/c] in world coords, end total energy [GeV] (vector<double> = {Px,Py,Pz,E})


SimTree -> G4 level info with 1 entry per particle
  - event:     event ID (int)
  - run:       run ID   (int)
  - subRun:    subRun ID (int)
  - nPoints:   number of trajectory points (int)
  - trackID:   Geant4 trackID (int)
  - pdg:       particle PDG code (e.g. 13 = muon, 11 = electron, 22 = photon)
  - trackLength: distance between start and end positions [cm] (float)
  - process:   code for process that produced the particle (int)
  - endProcess: code for process that killed the particle (int)
  - mother:    mother particle G4 trackID (int)
  - parentPDG: PDG code of parent particle (int)
  - parentE:   total energy of the parent particle [GeV] (float)
  - progenitor: particle's track ID from the generator that began the lineage (int)
  - nDaught:   number of daughters this particle produces (int)
  - startXYZT: particle start position [cm] in world coords, start time [ns] (vector<double> = {x,y,z,t})
  - endXYZT:   particle end position [cm] in world coords, end time [ns] (vector<double> = {x,y,z,t})
  - startPE:   particle start momentum components [GeV/c] in world coords, start total energy [GeV] (vector<double> = {Px,Py,Pz,E})
  - endPE:     particle end momentum components [GeV/c] in world coords, end total energy [GeV] (vector<double> = {Px,Py,Pz,E})

  // one entry per CRT strip hit by the particle
  - nChan:             number of CRT strips hit by this particle (int)
  - auxDetSensitiveID: id of the CRT strip (unique only within a CRT module) hit by the particle (vector<int>)
  - auxDetID:          id of the CRT module (globaly unique) hit by the particle (vector<int>)
  - auxDetEDep:        energy deposited in the scintillator strip [GeV] (vector<float>)
  - auxDetdEdx         <dE/dx> [GeV/cm] within the strip (vector<float>)
  - auxDetTrackLength: track length within the strip [cm] (vector<float>)
  - auxDetEnterXYZT:   track entry point/time (world coords) into the strip [cm, ns] (vector<vector<float>> = vector<{x,y,z,t}>) 
  - auxDetExitXYZT:    track exit point/time (world coords) from the strip [cm, ns] (vector<vector<float>> = vector<{x,y,z,t}>)
  - auxDetEnterPE:     track entry momentum components/energy (world coord) entering the strip [GeV/c, GeV] (vector<vector<float>> = vector<{Px,Py,Pz,E}>)
  - auxDetExitPE:      track entry momentum components/energy (world coord) entering the strip [GeV/c, GeV] (vector<vector<float>> = vector<{Px,Py,Pz,E}>)
  - auxDetRegion:      region code for the CRT strip (vector<int>)
  - mac5:              front-end board ID corresponding to the CRT module (vector<int>)
  - adType:            CRT module type with 0=CERN type [top], 1=MINOS type [side], 2=Double Chooz type [bottom] (vector<int>)

RegTree -> G4 level, particles entering/exiting different volumes with 1 entry per particle
  - event:     event ID (int)
  - run:       run ID   (int)
  - subRun:    subRun ID (int)
  - nReg:      number of regions crossed by particle (int)
  - fiducial:  number of fiducial LAr volumes crossed (int)
  - active:    number of active LAr volumes crossed (int)
  - inactive:  number of inactive LAr volumes crossed (int)
  - crts:      number of CRT regions crossed (int)
  - regions:   region codes crossed by particle, time ordered (vector<int>)
  - pdg:       PDG code of particle (int)
  - trackID:   Geant4 track ID of particle (int)
  - eDep:      energy deposited by particle in each region [GeV] (float)
  - dL:        path length within each module [cm] (float)
  - opDetID:   if region is cryostat, (globally unique) id of closest PMT to entry point in the volume (vector<int>)
  - distToOpDet: if region is cryostat, distance from entry point into the region to the closest PMT [cm] (vector<int>)
  - opDetXYZT:  if region is cryostat, position [cm] of PMT and (geometric) photon time of arrival [ns] at PMT (vector<vector<double>> = vector<{x,y,z,t}>)
  - entryPE:    track entry momentum components/energy (world coord) entering the region [GeV/c, GeV] (vector<vector<double>> = vector<{Px,Py,Pz,E}>)
  - exitPE:     track exit momentum components/energy (world coord) exiting the region [GeV/c, GeV] (vector<vector<double>> = vector<{Px,Py,Pz,E}>)
  - entryXYZT:  track entry point/time (world coords) into the region [cm, ns] (vector<vector<double>> = vector<{x,y,z,t}>)
  - exitXYZT:   track exit point/time (world coords) from the region [cm, ns] (vector<vector<double>> = vector<{x,y,z,t}>)
  - entrySlope: track entry slope (direction cosines) into the region  (vector<vector<double>> = vector<{Cosx,Cosy,Cosz}>)
  - exitSlope:  track exit slope (direction cosines) from the region  (vector<vector<double>> = vector<{Cosx,Cosy,Cosz}>)

CosmicDisplay -> G4 level, Trajectory points for event display
TO BE FILLED (FIX ME!)

DetTree -> CRT detector simulation info including FEB ID, charge, and timing
TO BE FILLED (FIX ME!)

HitTree -> CRT simulated hit reconstructon giving spatial, time, and CRT region
TO BE FILLED (FIX ME!)

CRTTrueHitTree -> truth level version of HitTree
TO BE FILLED (FIX ME!)

