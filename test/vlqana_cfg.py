import sys, os
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
#### for testing######################
### btageffmap=btagEff_Ttq30p_RH_loose.root applyBTagSFs=True cleanEvents=True isData=False doPUReweightingOfficial=True jecShift=0 jerShift=1 HTMin=0 storePreselEvts=True storeLHEWts=True storeTrigBits=False
###################################
options = VarParsing('analysis')
options.register('isData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Is data?"
    )
options.register('outFileName', 'singleT',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('hltORAND', 'OR',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "OR or AND HLT paths?"
    )
options.register('cleanEvents', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Clean events using EventCleaner?"
    )
options.register('doPUReweightingOfficial', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do pileup reweighting using official recipe"
    )
options.register('doBTagSFUnc', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply b-tag SF uncertainties"
    )
options.register('jecShift', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "JEC shift"
    )
options.register('jerShift', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "JER shift"
    )
options.register('topTagtau32', 0.5,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Top-tagging tau32 cut"
    )
options.register('topTagBDisc', 0.5426,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Top-tagging b-discriminator cut"
    )
options.register('HTMin', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum HT"
    )
options.register('storePreselEvts', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Store pre-selected events after pre-selection", 
    )
options.register('applyBTagSFs', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply b-tagging SFs to the MC"
    )
options.register('btageffmap', "btagEff_Ttq30p_RH_loose.root",#until new SFs arrive
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "ROOT file with Th2D histos of b tag effs for b,c, and light flavoured jets"
    )
options.register('doPreselOnly', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Only run pre-selections"
    )
options.register('storeLHEWts', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Store LHE wts?"
    )
options.register('storeTrigBits', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Store trigger bits and names?"
    )
options.register('FileNames', 'TbtH_1200_LH',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of list of input files"
    )

options.setDefault('maxEvents', -1000)
options.parseArguments()


#hltpaths = ["HLT_PFJet320_v"]
#Run2016B-G
#hltpathsOr = [
#               "HLT_AK8PFJet360_TrimMass30_v", 
#              "HLT_AK8DiPFJet280_200_TrimMass30_v",
#              "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v",
#              "HLT_AK8PFJet450_v",
#              "HLT_PFHT800_v",
#             ]
#Run2016H
#hltpathsOr = ["HLT_AK8PFJet360_TrimMass30_v", 
#              "HLT_AK8DiPFJet300_200_TrimMass30_v",
#              "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v",
#              "HLT_AK8PFJet450_v",
#              "HLT_PFHT900_v",
#             ]

hltpathsOr = ["HLT_AK8PFJet360_TrimMass30_v", 
              "HLT_AK8DiPFJet300_200_TrimMass30_v",
              "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v",
              "HLT_AK8PFJet450_v",
              "HLT_AK8DiPFJet280_200_TrimMass30_v",
              "HLT_PFHT800_v",
              "HLT_PFHT900_v",
             ]
if options.isData:
    options.doBTagSFUnc = False 
    options.jerShift = 0 
    options.doPUReweightingOfficial=False 
    options.storeLHEWts=False
    options.applyBTagSFs = False
#HTMin=1100
if options.storePreselEvts:
  HTMin = options.HTMin

process = cms.Process("VLQAna")

from inputFiles_cfi import * 
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      FileNames[options.FileNames]
      )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
      options.outFileName+".root"
      )
    )
#dataFilePath = "$CMSSW_BASE/src/Analysis/VLQAna/data/"
#dataFilePath = '../data/'
dataFilePath = './'
process.load("Analysis.VLQAna.EventCleaner_cff") 
process.evtcleaner.hltORAND = cms.string (options.hltORAND)  
process.evtcleaner.hltPaths = cms.vstring (hltpathsOr)  
process.evtcleaner.cleanEvents = cms.bool(options.cleanEvents)
process.evtcleaner.isData = options.isData 
process.evtcleaner.DoPUReweightingOfficial = options.doPUReweightingOfficial
process.evtcleaner.storeLHEWts = options.storeLHEWts
process.evtcleaner.storeTrigBits = options.storeTrigBits
process.evtcleaner.File_PUDistData      = cms.string(os.path.join(dataFilePath,'RunII2016Rereco_25ns_PUXsec69000nb.root'))
process.evtcleaner.File_PUDistDataLow   = cms.string(os.path.join(dataFilePath,'RunII2016Rereco_25ns_PUXsec65550nb.root'))
process.evtcleaner.File_PUDistDataHigh  = cms.string(os.path.join(dataFilePath,'RunII2016Rereco_25ns_PUXsec72450nb.root'))
process.evtcleaner.File_PUDistMC        = cms.string(os.path.join(dataFilePath,'PUDistMC_Summer2016_25ns_Moriond17MC_PoissonOOTPU.root'))

process.evtcleanerBG = process.evtcleaner.clone()
process.evtcleanerBG.File_PUDistData = cms.string(os.path.join(dataFilePath, 'RunII2016Rereco_25ns_RunsBtoG_PUXsec69000nb.root'))
process.evtcleanerBG.storeLHEWts = cms.bool(False)
process.evtcleanerBG.storeTrigBits = cms.bool(False)
process.evtcleanerH = process.evtcleaner.clone()
process.evtcleanerH.File_PUDistData = cms.string(os.path.join(dataFilePath,'RunII2016Rereco_25ns_RunH_PUXsec69000nb.root'))
process.evtcleanerH.storeLHEWts = cms.bool(False)
process.evtcleanerH.storeTrigBits = cms.bool(False)

process.evtcleanerAlt = process.evtcleaner.clone()
process.evtcleanerAlt.File_PUDistData      = cms.string(os.path.join(dataFilePath,'RunII2016_25ns_PUXsec72000nb.root'))
process.evtcleanerAlt.File_PUDistDataLow   = cms.string(os.path.join(dataFilePath,'RunII2016_25ns_PUXsec68400nb.root'))
process.evtcleanerAlt.File_PUDistDataHigh  = cms.string(os.path.join(dataFilePath,'RunII2016_25ns_PUXsec75600nb.root'))
process.evtcleanerAlt.storeLHEWts = cms.bool(False)
process.evtcleanerAlt.storeTrigBits = cms.bool(False)

from Analysis.VLQAna.VLQAna_cfi import *

if options.isData == False: ### Careful, to be reset when B2GAnaFW_v80X_v2.4 MC are used
  for par in ['jetAK4selParams', 'jetAK8selParams', 'jetHTaggedselParams', 'jetAntiHTaggedselParams', 'jetTopTaggedselParams', 'jetAntiTopTaggedselParams', 'jetZTaggedselParams', 'jetAntiZTaggedselParams']:
    if 'AK4' in par:
        jetType = 'AK4PFchs'
    else:
        jetType = 'AK8PFchs'
        payLoadTypes = ['L2Relative', 'L3Absolute']
        payLoadFiles = []
        for payLoadType in payLoadTypes: payLoadFiles.append(os.path.join(dataFilePath,'Summer16_23Sep2016V4_MC_'+payLoadType+'_'+jetType+'.txt'))   
        setattr(getattr(ana, par), 'jecAK8GroomedPayloadNames', cms.vstring(payLoadFiles))

        setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jettau1Label'         ,cms.InputTag("jetsAK8CHS", "jetAK8CHStau1CHS"))
        setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jettau2Label'         ,cms.InputTag("jetsAK8CHS", "jetAK8CHStau2CHS"))
        setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jettau3Label'         ,cms.InputTag("jetsAK8CHS", "jetAK8CHStau3CHS"))
        setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jetPrunedMassLabel'   ,cms.InputTag("jetsAK8CHS", "jetAK8CHSprunedMassCHS"))
        setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jetTrimmedMassLabel'  ,cms.InputTag("jetsAK8CHS", "jetAK8CHStrimmedMassCHS"))
        setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jetFilteredMassLabel' ,cms.InputTag("jetsAK8CHS", "jetAK8CHSfilteredMassCHS"))
        setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jetSoftDropMassLabel' ,cms.InputTag("jetsAK8CHS", "jetAK8CHSsoftDropMassCHS"))
    setattr(getattr(ana,par), 'jecUncPayloadName', cms.string(os.path.join(dataFilePath, 'Summer16_23Sep2016V4_MC_Uncertainty_'+jetType+'.txt')))

 
process.ana = ana.clone()
process.ana.doBTagSFUnc = options.doBTagSFUnc
process.ana.jetAK4selParams.jecShift = options.jecShift 
process.ana.jetAK4selParams.jerShift = options.jerShift 
process.ana.jetAK8selParams.jecShift = options.jecShift 
process.ana.jetAK8selParams.jerShift = options.jerShift 
process.ana.jetHTaggedselParams.jecShift = options.jecShift 
process.ana.jetHTaggedselParams.jerShift = options.jerShift 
process.ana.jetTopTaggedselParams.jecShift = options.jecShift 
process.ana.jetTopTaggedselParams.jerShift = options.jerShift 
process.ana.jetAntiHTaggedselParams.jecShift = options.jecShift 
process.ana.jetAntiHTaggedselParams.jerShift = options.jerShift 
process.ana.jetAntiTopTaggedselParams.jecShift = options.jecShift 
process.ana.jetAntiTopTaggedselParams.jerShift = options.jerShift 
process.ana.jetZTaggedselParams.jecShift = options.jecShift 
process.ana.jetZTaggedselParams.jerShift = options.jerShift 
process.ana.jetAntiZTaggedselParams.jecShift = options.jecShift 
process.ana.jetAntiZTaggedselParams.jerShift = options.jerShift 
process.ana.storePreselEvts = options.storePreselEvts
process.ana.doPreselOnly = options.doPreselOnly
process.ana.HTMin = HTMin
process.ana.applyBTagSFs = options.applyBTagSFs
process.ana.btageffmap = cms.string(os.path.join(dataFilePath,options.btageffmap))
process.ana.btageffmapMed = cms.string(os.path.join(dataFilePath,options.btageffmap.replace('loose','medium'))) 
process.ana.sjbtagSFcsv = cms.string(os.path.join(dataFilePath,"subjet_CSVv2_Moriond17_B_H.csv")) 



from Analysis.VLQAna.JetMaker_cfi import *
            

process.anaCHS = process.ana.clone(
    jetTopTaggedselParams = defaultCHSTJetSelectionParameters.clone(
       subjetHighestCSVMin = cms.double(CSVv2M),
       jettau3Bytau2Max = cms.double(0.57),
       ),
    jetAntiTopTaggedselParams = defaultCHSTJetSelectionParameters.clone(
       subjetHighestCSVMin = cms.double(-1000000),
       subjetHighestCSVMax = cms.double(CSVv2M),
      jettau3Bytau2Max = cms.double(0.57),
       ),
)

process.ana0p3 = process.ana.clone(
    jetTopTaggedselParams = defaultCHSTJetSelectionParameters.clone(
       jettau3Bytau2Max = cms.double(0.57),
    ),
    jetAntiTopTaggedselParams = defaultCHSTJetSelectionParameters.clone(
      subjetHighestCSVMin = cms.double(-1000000),
      subjetHighestCSVMax = defaultCHSTJetSelectionParameters.subjetHighestCSVMin, 
      jettau3Bytau2Max = cms.double(0.57),
    ),
)     
process.ana0p3Med = process.ana.clone(
    jetTopTaggedselParams = defaultCHSTJetSelectionParameters.clone(
       jettau3Bytau2Max = cms.double(0.57),
       subjetHighestCSVMin = cms.double(CSVv2M),
    ),
    jetAntiTopTaggedselParams = defaultCHSTJetSelectionParameters.clone(
      subjetHighestCSVMin = cms.double(-1000000),
      subjetHighestCSVMax = cms.double(CSVv2M), 
      jettau3Bytau2Max = cms.double(0.57),
    ),
)     
if options.isData:
    process.anaJERUp = process.ana.clone()
    process.anaJERDown = process.ana.clone()
    process.anaJESUp = process.ana.clone()
    process.anaJESDown = process.ana.clone()
    process.evtcleanerJERUp = process.evtcleaner.clone()
    process.evtcleanerJERUp.storeLHEWts = False
    process.evtcleanerJERUp.storeTrigBits = False
    process.evtcleanerJERDown = process.evtcleaner.clone()
    process.evtcleanerJERDown.storeLHEWts = False
    process.evtcleanerJERDown.storeTrigBits = False
    process.evtcleanerJESUp = process.evtcleaner.clone()
    process.evtcleanerJESUp.storeLHEWts = False
    process.evtcleanerJESUp.storeTrigBits = False
    process.evtcleanerJESDown = process.evtcleaner.clone()
    process.evtcleanerJESDown.storeLHEWts = False
    process.evtcleanerJESDown.storeTrigBits = False
else:    
    process.anaJERUp = process.ana.clone()
    process.anaJERUp.jetAK4selParams.jecShift = cms.double(0) 
    process.anaJERUp.jetAK4selParams.jerShift = cms.int32(2) 
    process.anaJERUp.jetAK8selParams.jecShift = cms.double(0) 
    process.anaJERUp.jetAK8selParams.jerShift = cms.int32(2) 
    process.anaJERUp.jetHTaggedselParams.jecShift = cms.double(0) 
    process.anaJERUp.jetHTaggedselParams.jerShift = cms.int32(2) 
    process.anaJERUp.jetTopTaggedselParams.jecShift = cms.double(0) 
    process.anaJERUp.jetTopTaggedselParams.jerShift = cms.int32(2) 
    process.anaJERUp.jetAntiHTaggedselParams.jecShift = cms.double(0) 
    process.anaJERUp.jetAntiHTaggedselParams.jerShift = cms.int32(2) 
    process.anaJERUp.jetAntiTopTaggedselParams.jecShift = cms.double(0) 
    process.anaJERUp.jetAntiTopTaggedselParams.jerShift = cms.int32(2) 
    process.anaJERUp.jetZTaggedselParams.jecShift = cms.double(0) 
    process.anaJERUp.jetZTaggedselParams.jerShift = cms.int32(2) 
    process.anaJERUp.jetAntiZTaggedselParams.jecShift = cms.double(0) 
    process.anaJERUp.jetAntiZTaggedselParams.jerShift = cms.int32(2) 
    
    process.evtcleanerJERUp = process.evtcleaner.clone()
    process.evtcleanerJERUp.storeLHEWts = False
    process.evtcleanerJERUp.storeTrigBits = False
    
     
    process.anaJERDown = process.ana.clone()
    process.anaJERDown.jetAK4selParams.jecShift = cms.double(0) 
    process.anaJERDown.jetAK4selParams.jerShift = cms.int32(-1) 
    process.anaJERDown.jetAK8selParams.jecShift = cms.double(0) 
    process.anaJERDown.jetAK8selParams.jerShift = cms.int32(-1) 
    process.anaJERDown.jetHTaggedselParams.jecShift = cms.double(0) 
    process.anaJERDown.jetHTaggedselParams.jerShift = cms.int32(-1) 
    process.anaJERDown.jetTopTaggedselParams.jecShift = cms.double(0) 
    process.anaJERDown.jetTopTaggedselParams.jerShift = cms.int32(-1) 
    process.anaJERDown.jetAntiHTaggedselParams.jecShift = cms.double(0) 
    process.anaJERDown.jetAntiHTaggedselParams.jerShift = cms.int32(-1) 
    process.anaJERDown.jetAntiTopTaggedselParams.jecShift = cms.double(0) 
    process.anaJERDown.jetAntiTopTaggedselParams.jerShift = cms.int32(-1) 
    process.anaJERDown.jetZTaggedselParams.jecShift = cms.double(0) 
    process.anaJERDown.jetZTaggedselParams.jerShift = cms.int32(-1) 
    process.anaJERDown.jetAntiZTaggedselParams.jecShift = cms.double(0) 
    process.anaJERDown.jetAntiZTaggedselParams.jerShift = cms.int32(-1)
    
    process.evtcleanerJERDown = process.evtcleaner.clone()
    process.evtcleanerJERDown.storeLHEWts = False
    process.evtcleanerJERDown.storeTrigBits = False
    
    process.anaJESUp = process.ana.clone()
    process.anaJESUp.jetAK4selParams.jecShift = cms.double(1) 
    process.anaJESUp.jetAK4selParams.jerShift = cms.int32(1) 
    process.anaJESUp.jetAK8selParams.jecShift = cms.double(1) 
    process.anaJESUp.jetAK8selParams.jerShift = cms.int32(1) 
    process.anaJESUp.jetHTaggedselParams.jecShift = cms.double(1) 
    process.anaJESUp.jetHTaggedselParams.jerShift = cms.int32(1) 
    process.anaJESUp.jetTopTaggedselParams.jecShift = cms.double(1) 
    process.anaJESUp.jetTopTaggedselParams.jerShift = cms.int32(1) 
    process.anaJESUp.jetAntiHTaggedselParams.jecShift = cms.double(1) 
    process.anaJESUp.jetAntiHTaggedselParams.jerShift = cms.int32(1) 
    process.anaJESUp.jetAntiTopTaggedselParams.jecShift = cms.double(1) 
    process.anaJESUp.jetAntiTopTaggedselParams.jerShift = cms.int32(1) 
    process.anaJESUp.jetZTaggedselParams.jecShift = cms.double(1) 
    process.anaJESUp.jetZTaggedselParams.jerShift = cms.int32(1) 
    process.anaJESUp.jetAntiZTaggedselParams.jecShift = cms.double(1) 
    process.anaJESUp.jetAntiZTaggedselParams.jerShift = cms.int32(1)
    
    process.evtcleanerJESUp = process.evtcleaner.clone()
    process.evtcleanerJESUp.storeLHEWts = False
    process.evtcleanerJESUp.storeTrigBits = False
    
    process.anaJESDown = process.ana.clone()
    process.anaJESDown.jetAK4selParams.jecShift = cms.double(-1) 
    process.anaJESDown.jetAK4selParams.jerShift = cms.int32(1) 
    process.anaJESDown.jetAK8selParams.jecShift = cms.double(-1) 
    process.anaJESDown.jetAK8selParams.jerShift = cms.int32(1) 
    process.anaJESDown.jetHTaggedselParams.jecShift = cms.double(-1) 
    process.anaJESDown.jetHTaggedselParams.jerShift = cms.int32(1) 
    process.anaJESDown.jetTopTaggedselParams.jecShift = cms.double(-1) 
    process.anaJESDown.jetTopTaggedselParams.jerShift = cms.int32(1) 
    process.anaJESDown.jetAntiHTaggedselParams.jecShift = cms.double(-1) 
    process.anaJESDown.jetAntiHTaggedselParams.jerShift = cms.int32(1) 
    process.anaJESDown.jetAntiTopTaggedselParams.jecShift = cms.double(-1) 
    process.anaJESDown.jetAntiTopTaggedselParams.jerShift = cms.int32(1) 
    process.anaJESDown.jetZTaggedselParams.jecShift = cms.double(-1) 
    process.anaJESDown.jetZTaggedselParams.jerShift = cms.int32(1) 
    process.anaJESDown.jetAntiZTaggedselParams.jecShift = cms.double(-1) 
    process.anaJESDown.jetAntiZTaggedselParams.jerShift = cms.int32(1)
     
    process.evtcleanerJESDown = process.evtcleaner.clone()
    process.evtcleanerJESDown.storeLHEWts = False
    process.evtcleanerJESDown.storeTrigBits = False

from Analysis.EventCounter.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone(isData=options.isData)
process.cleanedEvents = eventCounter.clone(isData=options.isData)
process.cleanedEventsJERUp = eventCounter.clone(isData=options.isData)
process.cleanedEventsJERDown = eventCounter.clone(isData=options.isData)
process.cleanedEventsJESUp = eventCounter.clone(isData=options.isData)
process.cleanedEventsJESDown = eventCounter.clone(isData=options.isData)
process.finalEvents = eventCounter.clone(isData=options.isData)
process.finalEvents0p3 = eventCounter.clone(isData=options.isData)
process.finalEvents0p3Med = eventCounter.clone(isData=options.isData)
process.finalEventsJERUp = eventCounter.clone(isData=options.isData)
process.finalEventsJERDown = eventCounter.clone(isData=options.isData)
process.finalEventsJESUp = eventCounter.clone(isData=options.isData)
process.finalEventsJESDown = eventCounter.clone(isData=options.isData)

process.p = cms.Path(
    process.allEvents
    *process.evtcleaner
    *process.evtcleanerBG
    *process.evtcleanerH
    *process.evtcleanerAlt
    *process.cleanedEvents
    *process.anaCHS 
    *process.finalEvents
    )
process.p0p3 = cms.Path(
    process.allEvents
    *process.evtcleaner
    *process.evtcleanerBG
    *process.evtcleanerH
    *process.cleanedEvents
    *process.ana0p3
    *process.finalEvents0p3
    )
process.p0p3Med = cms.Path(
    process.allEvents
    *process.evtcleaner
    *process.evtcleanerBG
    *process.evtcleanerH
    *process.cleanedEvents
    *process.ana0p3Med
    *process.finalEvents0p3Med
    )

process.pJERUp = cms.Path(
    process.allEvents
    *process.evtcleanerJERUp
    *process.cleanedEventsJERUp
    *process.anaJERUp 
    *process.finalEventsJERUp
    )
process.pJERDown = cms.Path(
    process.allEvents
    *process.evtcleanerJERDown
    *process.cleanedEventsJERDown
    *process.anaJERDown
    *process.finalEventsJERDown
    )
process.pJESUp = cms.Path(
    process.allEvents
    *process.evtcleanerJESUp
    *process.cleanedEventsJESUp
    *process.anaJESUp
    *process.finalEventsJESUp
    )
process.pJESDown = cms.Path(
    process.allEvents
    *process.evtcleanerJESDown
    *process.cleanedEventsJESDown
    *process.anaJESDown
    *process.finalEventsJESDown
    )
#process.out = cms.OutputModule("PoolOutputModule",
#        SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('evtcleaner')),
#            )

process.schedule = cms.Schedule(process.p, process.pJERUp, process.pJERDown, process.pJESUp, process.pJESDown)

#process.outpath = cms.EndPath(process.out)

open('dump.py','w').write(process.dumpPython())
