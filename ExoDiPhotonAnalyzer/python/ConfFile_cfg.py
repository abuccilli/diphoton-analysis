import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'
        #'root://cmsxrootd.fnal.gov///store/mc/RunIIFall15MiniAODv1/RSGravToGG_kMpl-001_M-1000_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/20000/784A23B1-2DAE-E511-83F0-00266CF89498.root'
        'root://cmsxrootd.fnal.gov//store/mc/RunIIFall15MiniAODv2/RSGravToGG_kMpl-001_M-1500_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/34B117D8-C1B8-E511-B935-782BCB20EDD2.root'
    )
)

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'

process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi");

process.demo = cms.EDAnalyzer(
'ExoDiPhotonAnalyzer',
photonsMiniAOD = cms.InputTag("slimmedPhotons"),
genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
rho = cms.InputTag("fixedGridRhoFastjetAll")
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("ExoDiphotonAnalyzer.root")
                                   )

process.p = cms.Path(process.demo)
