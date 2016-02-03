#include <iostream>
#include <memory>
#include <vector>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h" 

#include "Analysis/VLQAna/interface/PickGenPart.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include <TFile.h>
#include <TH1D.h>

class EventCleaner : public edm::EDFilter {
  public:
    explicit EventCleaner(const edm::ParameterSet&) ;
    ~EventCleaner() ; 
  private:
    virtual void beginJob() override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    double GetLumiWeightsPVBased (const std::string file, const std::string hist, const unsigned npv) ; 

    edm::LumiReWeighting LumiWeights_;

    edm::InputTag l_trigName                               ; 
    edm::InputTag l_trigBit                                ; 
    edm::InputTag l_metFiltersName                         ; 
    edm::InputTag l_metFiltersBit                          ; 
    edm::InputTag l_hbheNoiseFilter                        ; 
    edm::InputTag l_lheProd                                ; 
    std::string   l_genEvtInfoProd                         ; 
    edm::InputTag l_vtxRho                                 ; 
    edm::InputTag l_vtxZ                                   ; 
    edm::InputTag l_vtxNdf                                 ; 
    edm::InputTag l_npv                                    ;
    edm::InputTag l_puNtrueInt                             ;
    std::vector<std::string> hltPaths_                     ; 
    std::vector<std::string> metFilters_                   ; 
    const bool isData_                                     ; 
    const bool doPUReweightingOfficial_                    ;
    const bool doPUReweightingNPV_                         ;
    const std::string file_PVWt_                           ; 
    const std::string file_PUDistData_                     ;
    const std::string file_PUDistMC_                       ;
    const std::string hist_PVWt_                           ; 
    const std::string hist_PUDistData_                     ;
    const std::string hist_PUDistMC_                       ;
    edm::EDGetTokenT<LHEEventProduct> t_lheProd            ; 
    edm::EDGetTokenT<GenEventInfoProduct> t_genEvtInfoProd ; 

    edm::ParameterSet TtZParams_                           ;
    edm::ParameterSet TtHParams_                           ;
    edm::ParameterSet TbWParams_                           ;
    edm::ParameterSet BbZParams_                           ;
    edm::ParameterSet BbHParams_                           ;
    edm::ParameterSet BtWParams_                           ;

};

EventCleaner::EventCleaner(const edm::ParameterSet& iConfig) :
  l_trigName              (iConfig.getParameter<edm::InputTag>            ("trigNameLabel")),
  l_trigBit               (iConfig.getParameter<edm::InputTag>            ("trigBitLabel")),
  l_metFiltersName        (iConfig.getParameter<edm::InputTag>            ("metFiltersNameLabel")),
  l_metFiltersBit         (iConfig.getParameter<edm::InputTag>            ("metFiltersBitLabel")),
  l_hbheNoiseFilter       (iConfig.getParameter<edm::InputTag>            ("hbheNoiseFilterLabel")),
  l_lheProd               (iConfig.getParameter<std::string>              ("lheProdName")),
  l_genEvtInfoProd        (iConfig.getParameter<std::string>              ("genEvtInfoProdName")),
  l_vtxRho                (iConfig.getParameter<edm::InputTag>            ("vtxRhoLabel")),  
  l_vtxZ                  (iConfig.getParameter<edm::InputTag>            ("vtxZLabel")),  
  l_vtxNdf                (iConfig.getParameter<edm::InputTag>            ("vtxNdfLabel")),  
  l_npv                   (iConfig.getParameter<edm::InputTag>            ("npvLabel")),  
  l_puNtrueInt            (iConfig.getParameter<edm::InputTag>            ("puNtrueIntLabel")),  
  hltPaths_               (iConfig.getParameter<std::vector<std::string>> ("hltPaths")), 
  metFilters_             (iConfig.getParameter<std::vector<std::string>> ("metFilters")), 
  isData_                 (iConfig.getParameter<bool>                     ("isData")),
  doPUReweightingOfficial_(iConfig.getParameter<bool>                     ("DoPUReweightingOfficial")),
  doPUReweightingNPV_     (iConfig.getParameter<bool>                     ("DoPUReweightingNPV")),
  file_PVWt_              (iConfig.getParameter<std::string>              ("File_PVWt")),
  file_PUDistData_        (iConfig.getParameter<std::string>              ("File_PUDistData")),
  file_PUDistMC_          (iConfig.getParameter<std::string>              ("File_PUDistMC")),
  hist_PVWt_              (iConfig.getParameter<std::string>              ("Hist_PVWt")),
  hist_PUDistData_        (iConfig.getParameter<std::string>              ("Hist_PUDistData")),
  hist_PUDistMC_          (iConfig.getParameter<std::string>              ("Hist_PUDistMC")),
  TtZParams_              (iConfig.getParameter<edm::ParameterSet>        ("TtZParams")),  
  TtHParams_              (iConfig.getParameter<edm::ParameterSet>        ("TtHParams")),  
  TbWParams_              (iConfig.getParameter<edm::ParameterSet>        ("TbWParams")),  
  BbZParams_              (iConfig.getParameter<edm::ParameterSet>        ("BbZParams")),  
  BbHParams_              (iConfig.getParameter<edm::ParameterSet>        ("BbHParams")),  
  BtWParams_              (iConfig.getParameter<edm::ParameterSet>        ("BtWParams")) 
{
  if (doPUReweightingOfficial_) LumiWeights_ = edm::LumiReWeighting(file_PUDistMC_, file_PUDistData_, hist_PUDistMC_, hist_PUDistData_) ;
  t_lheProd = consumes<LHEEventProduct>(l_lheProd) ; 
  t_genEvtInfoProd = consumes<GenEventInfoProduct>(l_genEvtInfoProd) ; 
  produces<std::string>("evttype"); 
  produces<double>("evtwtGen");
  produces<double>("evtwtPV");
  produces<unsigned>("npv");
  produces<double>("partonBin");
}

EventCleaner::~EventCleaner() {}

bool EventCleaner::filter(edm::Event& evt, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;

  typedef Handle <vector<string>> hstring ; 
  typedef Handle <vector<float>> hfloat ; 
  typedef Handle <vector<int>> hint ; 

  //hstring h_generator             ; evt.getByLabel (l_generator              , h_generator            );
  hstring h_trigName              ; evt.getByLabel (l_trigName               , h_trigName             );
  hfloat  h_trigBit               ; evt.getByLabel (l_trigBit                , h_trigBit              ); 
  hstring h_metFiltersName        ; evt.getByLabel (l_metFiltersName         , h_metFiltersName       );
  hfloat  h_metFiltersBit         ; evt.getByLabel (l_metFiltersBit          , h_metFiltersBit        ); 
  hfloat  h_vtxRho                ; evt.getByLabel (l_vtxRho                 , h_vtxRho               );
  hfloat  h_vtxZ                  ; evt.getByLabel (l_vtxZ                   , h_vtxZ                 );
  hint    h_vtxNdf                ; evt.getByLabel (l_vtxNdf                 , h_vtxNdf               );
  Handle<int>   h_npv             ; evt.getByLabel (l_npv                    , h_npv                  );
  Handle<int>   h_puNtrueInt      ; evt.getByLabel (l_puNtrueInt             , h_puNtrueInt           );
  Handle <bool> h_hbheNoiseFilter ; evt.getByLabel (l_hbheNoiseFilter        , h_hbheNoiseFilter      );

  unsigned int hltdecisions(0) ; 
  for ( const string& myhltpath : hltPaths_ ) {
    vector<string>::const_iterator it ;
    for (it = (h_trigName.product())->begin(); it != (h_trigName.product())->end(); ++it ) {
      if ( it->find(myhltpath) < std::string::npos) {
        hltdecisions |= int((h_trigBit.product())->at( it - (h_trigName.product())->begin() )) << ( it - (h_trigName.product())->begin() ) ;  
      }
    }
  }
  if ( hltPaths_.size() > 0 && !hltdecisions) return false ; 

  bool hbheNoiseFilter(h_hbheNoiseFilter.product()) ; 
  if ( !hbheNoiseFilter ) return false ; 
  if ( isData_ ) {
    bool metfilterdecision(1) ; 
    for ( const string& metfilter : metFilters_ ) {
      vector<string>::const_iterator it ; 
      for (it = (h_metFiltersName.product())->begin(); it != (h_metFiltersName.product())->end(); ++it) {
        if ( it->find(metfilter) < std::string::npos) {
          metfilterdecision *= (h_metFiltersBit.product())->at( it - (h_metFiltersName.product())->begin() ) ; 
        }
      }
    }
    if ( !metfilterdecision ) return false ; 
  }

  const int npv(*h_npv) ; 
  //for ( unsigned ipv = 0; ipv < (h_vtxRho.product())->size(); ++ipv) {
  //  double vtxRho = (h_vtxho.product())->at(ipv) ; 
  //  double vtxZ = (h_vtxZ.product())->at(ipv) ; 
  //  double vtxNdf = (h_vtxNdf.product())->at(ipv) ; 
  //  if ( abs(vtxRho) < 2. && abs(vtxZ) <= 24. && vtxNdf > 4 ) ++npv ; 
  //}
  if ( npv < 1 ) return false ; 

  double evtwtGen(1.0) ; 
  double partonBin(0.0);
  if ( !isData_ ) {
    Handle<GenEventInfoProduct> h_genEvtInfoProd; 
    evt.getByToken(t_genEvtInfoProd, h_genEvtInfoProd);
    evtwtGen = h_genEvtInfoProd->weight() ; 
    evtwtGen /= abs(evtwtGen) ; 
    if (h_genEvtInfoProd->binningValues().size()>0) partonBin = h_genEvtInfoProd->binningValues()[0];
  }

  double evtwtPV(1.0) ;
  if ( !isData_ && doPUReweightingNPV_ ) evtwtPV *= GetLumiWeightsPVBased(file_PVWt_, hist_PVWt_, npv) ; 
  else if ( !isData_ && doPUReweightingOfficial_ ) { 
    int npuTrue(*h_puNtrueInt) ; 
    evtwtPV *= LumiWeights_.weight(npuTrue) ; 
  }

  string evttype(isData_ ? "EvtType_Data" : "EvtType_MC");

  if ( !isData_ ) {
    PickGenPart pickTtZ(TtZParams_) ; 
    PickGenPart pickTtH(TtHParams_) ; 
    PickGenPart pickTbW(TbWParams_) ; 
    PickGenPart pickBbZ(BbZParams_) ; 
    PickGenPart pickBbH(BbHParams_) ; 
    PickGenPart pickBtW(BtWParams_) ; 

    vlq::GenParticleCollection vlqTtZ = pickTtZ(evt) ;  
    vlq::GenParticleCollection vlqTtH = pickTtH(evt) ;  
    vlq::GenParticleCollection vlqTbW = pickTbW(evt) ;  
    vlq::GenParticleCollection vlqBbZ = pickBbZ(evt) ;  
    vlq::GenParticleCollection vlqBbH = pickBbH(evt) ;  
    vlq::GenParticleCollection vlqBtW = pickBtW(evt) ;  

    if ( vlqTtZ.size() == 2 ) evttype = "EvtType_MC_tZtZ" ; 
    if ( vlqTtH.size() == 2 ) evttype = "EvtType_MC_tHtH" ; 
    if ( vlqTbW.size() == 2 ) evttype = "EvtType_MC_bWbW" ; 

    if ( vlqTtZ.size() == 1 && vlqTtH.size() == 1 ) evttype = "EvtType_MC_tZtH" ; 
    if ( vlqTtZ.size() == 1 && vlqTbW.size() == 1 ) evttype = "EvtType_MC_tZbW" ; 
    if ( vlqTtH.size() == 1 && vlqTbW.size() == 1 ) evttype = "EvtType_MC_tHbW" ; 

    if ( vlqBbZ.size() == 2 ) evttype = "EvtType_MC_bZbZ" ; 
    if ( vlqBbH.size() == 2 ) evttype = "EvtType_MC_bHbH" ; 
    if ( vlqBtW.size() == 2 ) evttype = "EvtType_MC_tWtW" ; 

    if ( vlqBbZ.size() == 1 && vlqBbH.size() == 1 ) evttype = "EvtType_MC_bZbH" ; 
    if ( vlqBbZ.size() == 1 && vlqBtW.size() == 1 ) evttype = "EvtType_MC_bZtW" ; 
    if ( vlqBbH.size() == 1 && vlqBtW.size() == 1 ) evttype = "EvtType_MC_bHtW" ; 
  }

  auto_ptr<string>ptr_evttype(new string(evttype)); 
  auto_ptr<double>ptr_evtwtGen(new double(evtwtGen)); 
  auto_ptr<double>ptr_evtwtPV(new double(evtwtPV)); 
  auto_ptr<unsigned>ptr_npv(new unsigned(npv)); 
  auto_ptr<double>ptr_partonBin(new double(partonBin)); 

  evt.put(ptr_evttype, "evttype");
  evt.put(ptr_evtwtGen, "evtwtGen");
  evt.put(ptr_evtwtPV, "evtwtPV");
  evt.put(ptr_npv, "npv");
  evt.put(ptr_partonBin, "partonBin");

  return true ; 

}

void EventCleaner::beginJob() {
}

void EventCleaner::endJob() {
}

double EventCleaner::GetLumiWeightsPVBased (const std::string file, const std::string hist, const unsigned npv) { 
  double wtPU(1) ;
  TFile f(file.c_str()) ;
  TH1D* hwt = dynamic_cast<TH1D*>( f.Get( hist.c_str() ) ) ; 
  wtPU = npv > 0 && npv <= 60 ? hwt->GetBinContent(npv) : 1.; 
  delete hwt ; 
  f.Close() ; 
  return wtPU ;
}

DEFINE_FWK_MODULE(EventCleaner);
