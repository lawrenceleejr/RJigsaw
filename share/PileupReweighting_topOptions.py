# =============================================================================
#  Name:        PileupReweighting_topOptions.py
# 
#  Author:      Karsten Koeneke (DESY)
#  Created:     April 2011

#  Description: job options for pileup reweighting.
# 
#  Usage: Here, all neccessary job options for the DAOD can be set.
#         To run, type:
#              athena D2PDMaker/PileupReweighting_topOptions.py 2>&1 | tee log.txt
# =============================================================================


# ==============================================================================
# Include the job property flags for the RecExCommon package
# ==============================================================================
from RecExConfig.RecFlags import rec


# ==============================================================================
# If you have your own DPD Maker scripts
# (see:
#                    share/D2PD_ExampleSimple*.py
# for examples),
# then just append your script (wherever it is) to this list:
#
#       rec.DPDMakerScripts.append("MyPackage/MyScript")
#
# The example scripts are appended below, so you can see how it works!
# ==============================================================================
rec.DPDMakerScripts.append("PileupReweighting/DAOD_PileupReweighting.py")


# ==============================================================================
# Load your input file
# ==============================================================================
from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
#athenaCommonFlags.FilesInput=["/afs/cern.ch/atlas/project/PAT/data/data.physics_Egamma.AOD.pool.root"]
athenaCommonFlags.FilesInput=["/afs/cern.ch/user/k/kkoeneke/scratch0/data/AOD/mc10_7TeV.107652.AlpgenJimmyZeeNp2_pt20.merge.AOD.e737_s933_s946_r2215_r2260_tid299146_00/AOD.299146._000023.pool.root.1"]


# ==============================================================================
# You can change the location and name of your output file with
# these three flags:
# ==============================================================================
#D2PDFlags.OutputDirectoryName = "/my/directory/with/enough/space/"





# ==============================================================================
#           ---- NO NEED TO CHANGE ANYTHING BELOW THIS LINE !!! ----
# ==============================================================================


# ==============================================================================
# Set the number of events that you want to process
# ==============================================================================
#athenaCommonFlags.SkipEvents.set_Value_and_Lock(0)
if not vars().has_key('EvtMax'): EvtMax = -1
athenaCommonFlags.EvtMax=EvtMax


# ==============================================================================
# Execute the dpd maker
# ==============================================================================
from RecExConfig.RecFlags import rec
# Turn off most of RecExCommon... (optional)
rec.doCBNT       = False
rec.doWriteESD   = False
rec.doWriteAOD   = False
rec.doAOD        = False
rec.doWriteTAG   = False 
rec.doPerfMon    = False
rec.doHist       = False
rec.doTruth      = False
#rec.LoadGeometry = False
#rec.doApplyAODFix.set_Value_and_Lock(False)

# Mandatory for ESD->DPD or AOD->DPD or DPD->DPD
rec.doDPD   = True

# Print the used configuration to the log file
if rec.OutputLevel() <= INFO:
    print rec
    pass

# =============================================================================
# Final debugging options...
# Dump list of containers on StoreGate to output log
StoreGateSvc = Service( "StoreGateSvc" )
StoreGateSvc.Dump = True
svcMgr.StoreGateSvc.Dump = True     # Prints the content of StoreGate (EventStore) 
#svcMgr.MetaDataStore.Dump = True    # Prints output metadata store content
#svcMgr.InputMetaDataStore.Dump=True # Prints input metadata store content
#theApp.Dlls  += ["TruthExamples" ]
#theApp.TopAlg += [ "DumpMC" ]    # These two lines will dump some
#DumpMC = Algorithm ( "DumpMC" )  # info about the MC truth
#DumpMCG4.McEventKey = "SpclMC"
#topSequence.InsertTruthToolLooper.TrueJetInserter.OutputLevel = DEBUG 
# =============================================================================

# The job starts here!
include ("RecExCommon/RecExCommon_topOptions.py")


