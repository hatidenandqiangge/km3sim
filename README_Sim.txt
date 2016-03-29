The GEANT4 simulation part runs as:
./KM3Sim novrml $runnumber $output_evt_file $GDML_geometry_file $INPUT_Parameters_file null null ANTARES_EVT_FORMAT $input_evt_file
The arguments (except the self explaining) are:
1st argument: novrml, vrml, vrmle
              novrml : no vrml output
	      vrml   : vrml output for only the geometry, but for each event (run only for input file containing one event)
	      vrmle  : vrml output for the geometry, the hits and the muon only trajectories
              ATTENTION 1)RUN WITH VRML OUTPUT ONLY INTERACTIVELY, BECAUSE YOU MUST SUPPLY SCREEN INPUT WITH COLORS OF EACH VOLUME
                        2)VRML FILES TEND TO BE VERY LARGE, ESPECIALLY IF THE OMVOLUME AND/OR PMT VOLUME IS CHOSEN TO BE SHOWN
			3)the best is to choose only storey volume to be shown

2nd argument: random seed
4th argument: GDML input file
5th argument: file containing information, like QE, absorption and scattering lengths, .... vs photon wavelength
6th argument: EM parametrization tables (null here, because in full simulation are not used)
7th argument: HA parametrization tables (null here, because in full simulation are not used)
8th argument: always ANTARES_EVT_FORMAT to produce evt output
