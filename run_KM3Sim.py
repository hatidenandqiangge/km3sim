#! /usr/bin/env python

import os,sys
import glob
from optparse import OptionParser

#check the environment is set, if it isn't need to source setenv_KM3Sim.sh before executing
if 'KM3SIM_PATH' not in os.environ:
  print " **** environment not set. Need to 'source setenv_KM3Sim.csh' before running this ****  "
  exit(1)


#set the option parser
parser = OptionParser()


#parser.add_option("--nu-type",dest="nu_type",help="Specify neutrino type according to Apostolos convention: e, ebar, mu, mubar, tau, taubar.")
#parser.add_option("--run-number", dest="run_number", help="The run number of the file to be processed, it is also used to inizialise the random seed.")
parser.add_option("--in-file",dest="in_file",help="Full path to the input file.")
parser.add_option("--rnd-seed",dest="rnd_seed",default=1,help="The random seed. [default = %default]")
parser.add_option("--det-file", dest="det_file",help="The GDML detector geometry file. ")
parser.add_option("--properties-file", dest="properties_file",default=os.path.join(os.environ['KM3SIM_PATH'],"INPUTParametersRun_3500_scatter_WPD3.3_p0.0075_1.0"),help="File that defines the water and other properties.")
#parser.add_option("--prod-dir",dest="prod_dir",default=os.path.join(os.environ['KM3SIM_PATH'],"bartol_genie_gene"), help="The path where the production is stored. [default = %default]")
#parser.add_option("--prod-version",dest="prod_version", default = '20', help="The version of the bartol genie production. [default = %default]")
parser.add_option("--emin",dest="emin", default = '1',  help="The minimum neutrino energy generated. [default = %default]")
parser.add_option("--emax",dest="emax", default = '100', help="The maximum neutrino energy generated. [default = %default]")
parser.add_option("-j","--job",dest="job", default=False, action="store_true", help="Submit as batch job.")
parser.add_option("--in-dir",dest="in_dir",help="It must be a full path. All files in in_dir are submitted to be processed. Accepts .gz files. Only with -j option.")
parser.add_option("--out-dir",dest="out_dir",help="Full path to the directory where output files are copied after processing. Only with -j option.")
parser.add_option("--njobs",dest="njobs", type = 'int', default = sys.maxint,  help="Number of jobs to be subitted.")

(opt, args) = parser.parse_args()
if opt.det_file==None:
  parser.print_help()
  exit(1)




#check that detector file exists
if not os.path.exists(opt.det_file):
  print "ERROR: detector file "+opt.det_file+" does not exist, check paths."
  exit(1)

#check inputs
if opt.in_dir==None:
  if opt.in_file==None:
    parser.print_help()
    exit(1)
  #check infile exists
  if not os.path.exists(opt.in_file):
    print "ERROR: input file "+opt.in_file+" does not exist, check paths."
    exit(1)
else:
  if not os.path.exists(opt.in_dir):
    print "ERROR: input directory "+opt.in_dir+" does not exist, check paths."
    exit(1)

if not opt.in_file==None:
  infile  = opt.in_file
  print 'infile: '+infile

#submit it as a job?
if opt.job:
  if opt.out_dir==None or not os.path.exists(os.path.split(opt.out_dir)[0]):
    print "Need to specify a valid directory where to copy processed files with the option --out-dir."
    parser.print_help()
    exit(0)
  #create out dir if it does not exist
  os.system("mkdir -p "+opt.out_dir)
  #create log dir
  logdir=os.path.join(opt.out_dir,'log')
  os.system("mkdir -p "+logdir)
  #submit a whole directory of files?
  fidx = 0#the file index, used to increment the random seed
  if not opt.in_dir==None:
    infiles = glob.glob(os.path.join(opt.in_dir,"*evt"))
    njobs = 0
#    for root,dirs,files in os.walk(opt.in_dir):
    for infile in infiles:
        fidx+=1
        #infile = os.path.join(root,infile)
        
        ##GET THE SEED FOR THE DENSE PRODUCTION
        #seed = infile[infile.find('seed_')+len('seed_'):infile.find('seed_')+len('seed_')+2]
        #seed = seed+"%03d"%fidx
        ########################################
        
        ##GET THE SEED FOR THE gSeaGen PRODUCTION
        seed  =  infile.split(".")[-2]
        #########################################
        
        cmd = sys.argv[0]+' --in-file '+infile+' --out-dir '+opt.out_dir+' --rnd-seed '+seed+' --det-file '+opt.det_file+' --emin '+opt.emin+' --emax '+opt.emax+\
        ' --properties-file '+opt.properties_file
        logfile = os.path.join(logdir,"KM3Sim."+os.path.split(infile)[1]+'.log')
        job_name = 'mynameistanino'+str(fidx)
        qsub = 'qsub  -P P_antares -q long -V -l sps=1,ct=30:00:00,vmem=3G,fsize=20G -o '+logfile+' -j y -N '+job_name+" "+cmd#GE
        print qsub
        os.system(qsub)
        njobs+=1
        if njobs >= opt.njobs : break
    exit(0)
    
  cmd = sys.argv[0]+' --in-file '+opt.in_file+' --out-dir '+opt.out_dir+' --rnd-seed '+str(opt.rnd_seed)+' --det-file '+opt.det_file+' --emin '+opt.emin+' --emax '+opt.emax
  logfile = os.path.join(logdir,"sim"+os.path.split(infile)[1]+'.log')
  job_name = 'mynameisgiannino'
  qsub = 'qsub  -P P_antares -q long -V -l sps=1,ct=30:00:00,vmem=3G,fsize=20G -o '+logfile+' -j y -N '+job_name+" "+cmd#GE
  print qsub
  #os.system(qsub)
  exit(0)


#copy infile to local dir
os.system("'cp' %s ./"%infile)
infile  = os.path.split(infile)[1]
#unzip it
if infile[-3:]=='.gz':
  os.system("gunzip -f %s"%infile)
  infile = infile[:-3]
#copy the executable KM3Sim, the geometry and the parameter files to local dir
#copy KM3Sim executable
cmd = "'cp' "+os.path.join(os.environ['KM3SIM_PATH'],"KM3Sim")+" ./"
print cmd
os.system(cmd)
#copy geometry file
cmd = "'cp' "+opt.det_file+" ./"
print cmd
os.system(cmd)
opt.det_file = os.path.split(opt.det_file)[1]
#copy input parameter files
cmd = "'cp' "+opt.properties_file+" ./"
opt.properties_file = os.path.split(opt.properties_file)[1]
print cmd
os.system(cmd)
cmd = "'cp' "+os.path.join(os.environ['KM3SIM_PATH'],"MiePhaseFactors.in")+" ./"
print cmd
os.system(cmd)


#put the seed in the output file name
outfile = "KM3Sim."+os.path.split(infile)[1]
#outfile = outfile[:outfile.find('seed_')+len('seed_')]+str(opt.rnd_seed)+outfile[outfile.find('seed_')+len('seed_')+2:]

print 'outfile: '+outfile

#execute
cmd = "./KM3Sim novrml %s %s %s %s null null ANTARES_EVT_FORMAT %s"%(opt.rnd_seed,outfile,opt.det_file,opt.properties_file,infile)
print cmd
os.system(cmd)

#copy outfile to out_dir
if opt.out_dir!=None:
  cmd = 'cp -v '+outfile+' '+opt.out_dir
print cmd
os.system(cmd)

