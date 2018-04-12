#!/usr/bin/env python
"""
This is a small script that submits a config over many datasets
"""
import os
from optparse import OptionParser
from effMaps_cfg import *
def getOptions() :
    """
    Parse and return the arguments provided by the user.
    """
    usage = ('usage: python submit_all.py -c CONFIG -d DIR ')

    parser = OptionParser(usage=usage)    
    #parser.add_option("-c", "--config", dest="config",
    #    help=("The crab script you want to submit "),
    #    metavar="CONFIG")
    #parser.add_option("-d", "--dir", dest="dir",
    #    help=("The crab directory you want to use "),
    #    metavar="DIR")
    parser.add_option("-f", "--datasets", dest="datasets",
        help=("File listing datasets to run over"),
        metavar="FILE")
    (options, args) = parser.parse_args()


    #if options.config == None or options.dir == None:
     #   parser.error(usage)
    
    return options
    

def main():
    options = getOptions()

    from WMCore.Configuration import Configuration
    config = Configuration()

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.section_("General")
    config.General.workArea = '80X_trees_3p2_10Apr18'
    config.General.transferLogs = True

    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = 'vlqana_cfg.py' 
    config.JobType.pyCfgParams = None 
    config.JobType.inputFiles = [
        '../data/PUDistMC_Summer2016_25ns_Moriond17MC_PoissonOOTPU.root'
        ,'../data/RunII2016Rereco_25ns_PUXsec69000nb.root'
        ,'../data/RunII2016Rereco_25ns_PUXsec65550nb.root'
        ,'../data/RunII2016Rereco_25ns_PUXsec72450nb.root'
        ,'../data/RunII2016_25ns_PUXsec72000nb.root'
        ,'../data/RunII2016_25ns_PUXsec68400nb.root'
        ,'../data/RunII2016_25ns_PUXsec75600nb.root'
        ,'../data/RunII2016Rereco_25ns_RunsBtoG_PUXsec69000nb.root'
        ,'../data/RunII2016Rereco_25ns_RunH_PUXsec69000nb.root'
        ,'../data/Summer16_23Sep2016V4_MC_L2Relative_AK8PFPuppi.txt'
        ,'../data/Summer16_23Sep2016V4_MC_L3Absolute_AK8PFPuppi.txt'
        ,'../data/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt'
        ,'../data/Summer16_23Sep2016V4_MC_L2Relative_AK8PFchs.txt'
        ,'../data/Summer16_23Sep2016V4_MC_L3Absolute_AK8PFchs.txt'
        ,'../data/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt'
        ,'../data/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFPuppi.txt'
        #,'../data/Summer16_23Sep2016V4_DATA_L2L3Residual_AK4PFchs.txt' 
        #,'../data/Summer16_23Sep2016V4_DATA_L2L3Residual_AK8PFchs.txt' 
        #,'../data/Summer16_23Sep2016V4_DATA_L2Relative_AK8PFchs.txt' 
        #,'../data/Summer16_23Sep2016V4_DATA_L3Absolute_AK8PFchs.txt' 
        #,'../data/Summer16_23Sep2016V4_DATA_Uncertainty_AK4PFchs.txt' 
        #,'../data/Summer16_23Sep2016V4_DATA_Uncertainty_AK8PFchs.txt' 
        ,'../data/btagEff_TTJets_loose.root'
        ,'../data/btagEff_WJets_loose.root'
        ,'../data/btagEff_ST_loose.root'
        ,'../data/btagEff_ZJets_loose.root'
        ,'../data/btagEff_ttWJets_loose.root'
        ,'../data/btagEff_ttHJets_loose.root'
        ,'../data/btagEff_ttZJets_loose.root'
        ,'../data/btagEff_WW_loose.root'
        ,'../data/btagEff_WZ_loose.root'
        ,'../data/btagEff_ZZ_loose.root'
        ,'../data/btagEff_tHQ_loose.root'
        ,'../data/btagEff_Tbq_LH_loose.root'
        ,'../data/btagEff_Tbq10p_LH_loose.root'
        ,'../data/btagEff_Tbq20p_LH_loose.root'
        ,'../data/btagEff_Tbq30p_LH_loose.root'
        ,'../data/btagEff_Tbq_tZ_RH_loose.root'
        ,'../data/btagEff_Ttq_tZ_RH_loose.root'
        ,'../data/btagEff_Tbq_tZ_LH_loose.root'
        ,'../data/btagEff_Ttq_LH_loose.root'
        ,'../data/btagEff_Tbq_RH_loose.root'
        ,'../data/btagEff_Ttq_RH_loose.root'
        ,'../data/btagEff_Ttq10p_RH_loose.root'
        ,'../data/btagEff_Ttq20p_RH_loose.root'
        ,'../data/btagEff_Ttq30p_RH_loose.root'
        ,'../data/btagEff_QCDPt_loose.root'
        ,'../data/btagEff_QCDHT_loose.root'
        ,'../data/btagEff_TTJets_medium.root'
        ,'../data/btagEff_WJets_medium.root'
        ,'../data/btagEff_ST_medium.root'
        ,'../data/btagEff_ZJets_medium.root'
        ,'../data/btagEff_ttWJets_medium.root'
        ,'../data/btagEff_ttHJets_medium.root'
        ,'../data/btagEff_ttZJets_medium.root'
        ,'../data/btagEff_WW_medium.root'
        ,'../data/btagEff_WZ_medium.root'
        ,'../data/btagEff_ZZ_medium.root'
        ,'../data/btagEff_tHQ_medium.root'
        ,'../data/btagEff_Tbq_LH_medium.root'
        ,'../data/btagEff_Tbq10p_LH_medium.root'
        ,'../data/btagEff_Tbq20p_LH_medium.root'
        ,'../data/btagEff_Tbq30p_LH_medium.root'
        ,'../data/btagEff_Tbq_tZ_RH_medium.root'
        ,'../data/btagEff_Ttq_tZ_RH_medium.root'
        ,'../data/btagEff_Tbq_tZ_LH_medium.root'
        ,'../data/btagEff_Ttq_LH_medium.root'
        ,'../data/btagEff_Tbq_RH_medium.root'
        ,'../data/btagEff_Ttq_RH_medium.root'
        ,'../data/btagEff_Ttq10p_RH_medium.root'
        ,'../data/btagEff_Ttq20p_RH_medium.root'
        ,'../data/btagEff_Ttq30p_RH_medium.root'
        ,'../data/btagEff_QCDPt_medium.root'
        ,'../data/btagEff_QCDHT_medium.root'
        ,'../data/subjet_CSVv2_Moriond17_B_H.csv'
        ]
    config.JobType.maxJobRuntimeMin = 2000
    config.JobType.maxMemoryMB = 2500

    config.section_("Data")
    config.Data.inputDataset = None
    config.Data.inputDBS = 'phys03'
    #config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
    config.Data.splitting = 'FileBased'
    #config.Data.splitting = 'LumiBased'
    config.Data.unitsPerJob = 1
    config.Data.ignoreLocality = False
    config.Data.publication = False     
    config.Data.outLFNDirBase = '/store/user/eschmitz/B2G/80X_VLQAna_B2G3p2_10Apr18'
    
    config.section_("Site")
    config.Site.storageSite = 'T2_US_Nebraska'

    #print 'Using config ' + options.config
    #print 'Writing to directory ' + options.dir
    
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print 'Cannot execute commend'
            print hte.headers

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    datasetsFile = open( options.datasets )
    jobsLines = datasetsFile.readlines()
    jobs = []
    for ijob in jobsLines :
        s = ijob.rstrip()
        if (len(s)==0 or s[0][0]=='#'): continue
        jobs.append( s )
        print '  --> added ' + s
        
    for ijob, job in enumerate(jobs) :

        print "-->  ", job
        pd = job.split('/')[1] + job.split('/')[2].split('-')[0]
        processing = (job.split('/')[2]).split('-')[0] + (job.split('/')[2]).split('-')[1] + (job.split('/')[2]).split('-')[2] + (job.split('/')[2]).split('-')[0] #for data
        if (len(pd + '_' + processing)<=100): 
          config.General.requestName = pd + '_' + processing
        else:
          config.General.requestName = pd
        pyCfgParams = ['btageffmap='+effMap[job], 'applyBTagSFs=True','cleanEvents=False','isData=False','doPUReweightingOfficial=True','jecShift=0','jerShift=1', 'HTMin=0', 'storePreselEvts=True', 'storeLHEWts=False', 'storeTrigBits=False'] 
        #pyCfgParams = ['btageffmap='+effMap[job],'applyBTagSFs=True','cleanEvents=True','isData=False','doPUReweightingOfficial=True','jecShift=0','jerShift=1', 'HTMin=0', 'storePreselEvts=True', 'storeLHEWts=True', 'storeTrigBits=False'] 
        config.JobType.pyCfgParams = pyCfgParams
        config.Data.inputDataset = job
        print 'Submitting ' + config.General.requestName + ', dataset = ' + job
        print 'Configuration :'
        print config
        try :
            from multiprocessing import Process
            p = Process(target=submit, args=(config,))
            p.start()
            p.join()
            #submit(config)
        except :
            print 'Not submitted.'

if __name__ == '__main__':
    main()            
