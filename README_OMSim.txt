The OM simulation part runs as:

./OmSim $inputfile null null k40_multiple_hits PmtPositionsAndDirections_info $k40 $outputfile $runnumber hit_raw
Arguments are:
2nd, 3rd: input files for parametrized single pe pulse shape and response function (are not used usually)
4th: file containing k40 multiple corellated info
5th: info of pmt positions and directions (see README_Geom.txt)
6th: == background , to include k40 background , or == nobackground
9th: == hit_raw , to output hit info as in evt format , or == hit_rawOM to output OMhits info (hits of all pmts in an OM 
                                                              within a time window are merged to produce a single OMhit, with
							      multiplicity information and OM direction information)
