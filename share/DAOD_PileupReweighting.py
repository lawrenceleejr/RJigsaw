# -----------------------------------------------------------------------------
# Name: DAOD_PileupReweighting.py
#
# Author: Karsten Koeneke (CERN)
# Email:  karsten.koeneke@cernSPAMNOT.ch
#
# Description: This calculates the pileup reweighting factors and attaches
#              them to each event as UserData.
# -----------------------------------------------------------------------------

# Include the job property flags for this package and from RecExCommon
from D2PDMaker.D2PDFlags import D2PDFlags

# This handels multiple output streams
from OutputStreamAthenaPool.MultipleStreamManager import MSMgr 


# ====================================================================
# Check if we have Monte Carlo or real data, based on the inputFileSummary
# ====================================================================
inputIsSimulation = False
if inputFileSummary.has_key("evt_type"):
    eventTypeList = inputFileSummary.get("evt_type")
    if eventTypeList.__contains__("IS_SIMULATION") :
        print "Detected that the input file is a simulated dataset"
        inputIsSimulation = True
        pass
    pass


if inputIsSimulation:
    # ====================================================================
    # Try to copy the root files locally (DON'T use https!)
    # ====================================================================
    print "Copying root files with mu distributions locally..."
    import socket
    socket.setdefaulttimeout(180)
    remoteDataFile = "http://atlas.web.cern.ch/Atlas/GROUPS/DATAPREPARATION/InteractionsperCrossing/data11_7TeV.periodB.root"
    localDataFile  = "data11_7TeV.periodB.root"
    remoteMCFile = "http://atlas.web.cern.ch/Atlas/GROUPS/DATAPREPARATION/InteractionsperCrossing/mu_mc10a.root"
    localMCFile  = "mu_mc10a.root"
    import urllib
    urllib.urlretrieve( remoteDataFile, localDataFile )
    urllib.urlretrieve( remoteMCFile, localMCFile )
    
    # ====================================================================
    # Load the pileup reweighting algorithm
    # ====================================================================
    from PileupReweighting.PileupReweightingConf import PileupReweightingAlg
    topSequence += PileupReweightingAlg( "PileupStream_PileupReweightingAlg",
                                         OutputLevel      = DEBUG,
                                         EventInfoName    = "MyEvent",
                                         dataROOTFileName = localDataFile,
                                         dataROOTHistName = "intperbx",
                                         mcROOTFileName   = localMCFile,
                                         mcROOTHistName   = "mu_mc10a"
                                         )
    pass




# ====================================================================
# Define the test DPD output stream
# ====================================================================
# This stream HAS TO start with "StreamD2AODM_"! If the input was an (D)ESD(M), this should start with "StreamD2ESD(M)_".
# See this twiki for more information: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/DPDNamingConvention
fileName   = "DAOD_WithPileupWeight.pool.root"
streamName = "DAOD_TEST"
PielupStream = MSMgr.NewPoolStream( streamName, fileName )

# Only events that pass the filters listed below are written out
# AcceptAlgs  = logical OR of filters
# RequireAlgs = logical AND of filters
#PileupStream.AcceptAlgs( [""] )



# ---------------------------------------------------
# Add the containers to the output stream
# ---------------------------------------------------
from PrimaryDPDMaker import PrimaryDPD_OutputDefinitions as dpdOutput

# Take all items from the input, except for the ones listed in the excludeList
# If the excludeList is empty, all containers from the input file (e.g. AOD)
# are copied to the output file.
excludeList = [ "TrigMuonEFContainer#HLT_MuonEF" ]
dpdOutput.addAllItemsFromInputExceptExcludeList( streamName, excludeList )
