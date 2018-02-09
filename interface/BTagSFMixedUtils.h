#ifndef ANALYSIS_INTERFACE_BTAGSFMIXEDUTILS_H
#define ANALYSIS_INTERFACE_BTAGSFMIXEDUTILS_H
#include <TF1.h>
#include <TFile.h>
#include <TH2D.h>
#include <map>

#include "Analysis/VLQAna/interface/Utilities.h"
#include "Analysis/VLQAna/interface/BTagCalibrationStandalone.h"
#include <boost/math/special_functions/fpclassify.hpp>

class BTagSFMixedUtils {

  private:

    std::unique_ptr<BTagCalibration>        calib_      ; 
    std::unique_ptr<BTagCalibrationReader>  reader_     ; 
    std::unique_ptr<BTagCalibrationReader>  readerUp_   ; 
    std::unique_ptr<BTagCalibrationReader>  readerDown_ ; 
    std::unique_ptr<BTagCalibrationReader>  reader2_     ; 
    std::unique_ptr<BTagCalibrationReader>  readerUp2_   ; 
    std::unique_ptr<BTagCalibrationReader>  readerDown2_ ; 

    BTagEntry::OperatingPoint op_ ; 
    BTagEntry::OperatingPoint op2_ ; 
    const double bfl_ptMin_ ;
    const double bfl_ptMax_ ;
    const double cfl_ptMin_ ;
    const double cfl_ptMax_ ;
    const double lfl_ptMin_ ;
    const double lfl_ptMax_ ;

    const std::string btageffmap_ ; 
    const std::string btageffmap2nd_ ; 

  public: 

    BTagSFMixedUtils (const std::string reader, 
        const BTagEntry::OperatingPoint op,
        const BTagEntry::OperatingPoint op2,
        const double bfl_ptMin, const double bfl_ptMax,
        const double cfl_ptMin, const double cfl_ptMax,
        const double lfl_ptMin, const double lfl_ptMax, 
        const std::string btageffmap,
        const std::string btageffmap2nd
        ) : 
      calib_          (new BTagCalibration("CSVv2",reader)), 
      reader_         (new BTagCalibrationReader(op,"central")), 
      readerUp_       (new BTagCalibrationReader(op,"up")), 
      readerDown_     (new BTagCalibrationReader(op,"down")),  
      reader2_         (new BTagCalibrationReader(op2,"central")), 
      readerUp2_       (new BTagCalibrationReader(op2,"up")), 
      readerDown2_     (new BTagCalibrationReader(op2,"down")),  
      op_(op),
      op2_(op2),
      bfl_ptMin_ (bfl_ptMin),
      bfl_ptMax_ (bfl_ptMax),
      cfl_ptMin_ (cfl_ptMin),
      cfl_ptMax_ (cfl_ptMax),
      lfl_ptMin_ (lfl_ptMin),
      lfl_ptMax_ (lfl_ptMax), 
      btageffmap_(btageffmap),
      btageffmap2nd_(btageffmap2nd)
  {
    reader_->load(*calib_,     // calibration instance
        BTagEntry::BTagEntry::FLAV_B,     // btag flavour
        "mujets") ;            // measurement type

    reader_->load(*calib_,     // calibration instance
        BTagEntry::BTagEntry::FLAV_C,     // btag flavour
        "mujets") ;            // measurement type

    reader_->load(*calib_,     // calibration instance
        BTagEntry::BTagEntry::FLAV_UDSG,  // btag flavour
        "incl") ;              // measurement type

    readerUp_->load(*calib_,   // calibration instance
        BTagEntry::BTagEntry::FLAV_B,     // btag flavour
        "mujets") ;            // measurement type

    readerUp_->load(*calib_,   // calibration instance
        BTagEntry::BTagEntry::FLAV_C,     // btag flavour
        "mujets") ;            // measurement type

    readerUp_->load(*calib_,   // calibration instance
        BTagEntry::BTagEntry::FLAV_UDSG,  // btag flavour
        "incl") ;              // measurement type

    readerDown_->load(*calib_, // calibration instance
        BTagEntry::BTagEntry::FLAV_B,     // btag flavour
        "mujets") ;            // measurement type

    readerDown_->load(*calib_, // calibration instance
        BTagEntry::BTagEntry::FLAV_C,     // btag flavour
        "mujets") ;            // measurement type

    readerDown_->load(*calib_, // calibration instance
        BTagEntry::BTagEntry::FLAV_UDSG,  // btag flavour
        "incl") ;              // measurement type

    reader2_->load(*calib_,     // calibration instance
        BTagEntry::BTagEntry::FLAV_B,     // btag flavour
        "mujets") ;            // measurement type

    reader2_->load(*calib_,     // calibration instance
        BTagEntry::BTagEntry::FLAV_C,     // btag flavour
        "mujets") ;            // measurement type

    reader2_->load(*calib_,     // calibration instance
        BTagEntry::BTagEntry::FLAV_UDSG,  // btag flavour
        "incl") ;              // measurement type

    readerUp2_->load(*calib_,   // calibration instance
        BTagEntry::BTagEntry::FLAV_B,     // btag flavour
        "mujets") ;            // measurement type

    readerUp2_->load(*calib_,   // calibration instance
        BTagEntry::BTagEntry::FLAV_C,     // btag flavour
        "mujets") ;            // measurement type

    readerUp2_->load(*calib_,   // calibration instance
        BTagEntry::BTagEntry::FLAV_UDSG,  // btag flavour
        "incl") ;              // measurement type

    readerDown2_->load(*calib_, // calibration instance
        BTagEntry::BTagEntry::FLAV_B,     // btag flavour
        "mujets") ;            // measurement type

    readerDown2_->load(*calib_, // calibration instance
        BTagEntry::BTagEntry::FLAV_C,     // btag flavour
        "mujets") ;            // measurement type

    readerDown2_->load(*calib_, // calibration instance
        BTagEntry::BTagEntry::FLAV_UDSG,  // btag flavour
        "incl") ;              // measurement type
  }

    BTagSFMixedUtils () : 
      calib_          (new BTagCalibration("CSVv2","subjet_CSVv2_Moriond17_B_H.csv")), 
      reader_         (new BTagCalibrationReader(BTagEntry::OP_LOOSE,"central")), 
      readerUp_       (new BTagCalibrationReader(BTagEntry::OP_LOOSE,"up")), 
      readerDown_     (new BTagCalibrationReader(BTagEntry::OP_LOOSE,"down")),
      reader2_         (new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central")), 
      readerUp2_       (new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"up")), 
      readerDown2_     (new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"down")),
      op_(BTagEntry::OP_LOOSE),
      op2_(BTagEntry::OP_MEDIUM),
      bfl_ptMin_ (20.),
      bfl_ptMax_ (450.),
      cfl_ptMin_ (20.),
      cfl_ptMax_ (450.),
      lfl_ptMin_ (20.),
      lfl_ptMax_ (1000.) 
  {
    reader_->load(*calib_,     // calibration instance
        BTagEntry::BTagEntry::FLAV_B,     // btag flavour
        "mujets") ;            // measurement type

    reader_->load(*calib_,     // calibration instance
        BTagEntry::BTagEntry::FLAV_C,     // btag flavour
        "mujets") ;            // measurement type

    reader_->load(*calib_,     // calibration instance
        BTagEntry::BTagEntry::FLAV_UDSG,  // btag flavour
        "incl") ;              // measurement type

    readerUp_->load(*calib_,   // calibration instance
        BTagEntry::BTagEntry::FLAV_B,     // btag flavour
        "mujets") ;            // measurement type

    readerUp_->load(*calib_,   // calibration instance
        BTagEntry::BTagEntry::FLAV_C,     // btag flavour
        "mujets") ;            // measurement type

    readerUp_->load(*calib_,   // calibration instance
        BTagEntry::BTagEntry::FLAV_UDSG,  // btag flavour
        "incl") ;              // measurement type

    readerDown_->load(*calib_, // calibration instance
        BTagEntry::BTagEntry::FLAV_B,     // btag flavour
        "mujets") ;            // measurement type

    readerDown_->load(*calib_, // calibration instance
        BTagEntry::BTagEntry::FLAV_C,     // btag flavour
        "mujets") ;            // measurement type

    readerDown_->load(*calib_, // calibration instance
        BTagEntry::BTagEntry::FLAV_UDSG,  // btag flavour
        "incl") ;              // measurement type


    reader2_->load(*calib_,     // calibration instance
        BTagEntry::BTagEntry::FLAV_B,     // btag flavour
        "mujets") ;            // measurement type

    reader2_->load(*calib_,     // calibration instance
        BTagEntry::BTagEntry::FLAV_C,     // btag flavour
        "mujets") ;            // measurement type

    reader2_->load(*calib_,     // calibration instance
        BTagEntry::BTagEntry::FLAV_UDSG,  // btag flavour
        "incl") ;              // measurement type

    readerUp2_->load(*calib_,   // calibration instance
        BTagEntry::BTagEntry::FLAV_B,     // btag flavour
        "mujets") ;            // measurement type

    readerUp2_->load(*calib_,   // calibration instance
        BTagEntry::BTagEntry::FLAV_C,     // btag flavour
        "mujets") ;            // measurement type

    readerUp2_->load(*calib_,   // calibration instance
        BTagEntry::BTagEntry::FLAV_UDSG,  // btag flavour
        "incl") ;              // measurement type

    readerDown2_->load(*calib_, // calibration instance
        BTagEntry::BTagEntry::FLAV_B,     // btag flavour
        "mujets") ;            // measurement type

    readerDown2_->load(*calib_, // calibration instance
        BTagEntry::BTagEntry::FLAV_C,     // btag flavour
        "mujets") ;            // measurement type

    readerDown2_->load(*calib_, // calibration instance
        BTagEntry::BTagEntry::FLAV_UDSG,  // btag flavour
        "incl") ;              // measurement type

  }

    void getBTagSFs(
        std::vector<double>jetcsvs,
        std::vector<double>jetpts,
        std::vector<double>jetetas,
        std::vector<int>jetflhads,
        const double csvMin, 
        const double csvMin2, 
        double& btagsf, double& btagsf_bcUp, double& btagsf_bcDown, double& btagsf_lUp, double& btagsf_lDown
        ) {

      //// Assign flavour notation as per b-tag SF .csv file 
      std::vector<BTagEntry::JetFlavor> jetfls ; 
      for ( int flhad : jetflhads ) {
        BTagEntry::JetFlavor fl ; 
        if ( abs(flhad) == 5 ) fl = BTagEntry::FLAV_B;
        else if ( abs(flhad) == 4 ) fl = BTagEntry::FLAV_C;
        else fl = BTagEntry::FLAV_UDSG;
        jetfls.push_back(fl) ; 
      }

      //// Get pT within bounds and double uncert if needed
      std::vector<int> uncscales ; 
      for ( auto idx : index(jetpts) ) {
        BTagEntry::JetFlavor fl = jetfls.at(idx.first) ; 
        double pt = jetpts.at(idx.first) ; 
        if (pt < 0.001) continue;
        double uncscale(1.) ; 
        if ( fl == 0 || fl == 1) {
          if ( pt < bfl_ptMin_ ) { pt = 20.01 ; uncscale *= 2 ; } 
          if ( pt > bfl_ptMax_ ) { pt = 449.99 ; uncscale *= 2 ; } 
          if ( fl == 1 ) uncscale *= 2 ; 
        }
        else {
          if ( pt < lfl_ptMin_ ) { pt =20.01 ; uncscale *= 2 ; }
          if ( pt > lfl_ptMax_ ) { pt = 999.99 ; uncscale *= 2; }
        }
        jetpts.at(idx.first) = pt ; 
        uncscales.push_back(uncscale) ; 
      }

      //// Get scale factors and efficiencies 
      std::vector<double> sfs, sfsUp, sfsDown, dsfsUp, dsfsDown, effs ;
      for ( auto idx : index(jetpts) ) {
        BTagEntry::JetFlavor fl = jetfls.at(idx.first) ; 
        int flhad = jetflhads.at(idx.first) ; 
        double pt = jetpts.at(idx.first) ; 
        double eta = jetetas.at(idx.first) ; 
        double sf = reader_->eval(fl,eta,pt) ; 
        double sfUp = readerUp_->eval(fl,eta,pt) ; 
        double sfDown = readerDown_->eval(fl,eta,pt) ; 
        double sf2 = reader2_->eval(fl,eta,pt) ; 
        double sfUp2 = readerUp2_->eval(fl,eta,pt) ; 
        double sfDown2 = readerDown2_->eval(fl,eta,pt) ; 
        double dsfUp = sfUp - sf ; 
        double dsfDown = sfDown - sf ;
        double dsfUp2 = sfUp2 - sf2 ; 
        double dsfDown2 = sfDown2 - sf2 ;
        double eff = 1.0;
        double eff2 = 1.0;
        if (pt < 0.001) continue;
        if (abs(eta) > 2.4) continue;
        eff = BTagSFMixedUtils::getBTagEff(pt,eta,flhad) ; 
        eff2 = BTagSFMixedUtils::getBTagEff2nd(pt,eta,flhad) ; 
        double uncscale = uncscales.at(idx.first) ; 
        sfs.push_back(sf) ; 
        sfsUp.push_back(sfUp) ; 
        sfsDown.push_back(sfDown) ; 
        dsfsUp.push_back(dsfUp) ; 
        dsfsDown.push_back(dsfDown) ; 
        effs.push_back(eff) ; 

        double sfUpabs = sf != 0. ? sf  + (uncscale*dsfUp/sf) : 0. ; 
        double sfDownabs = sf != 0. ? sf  + (uncscale*dsfDown/sf) : 0. ; 
        double sfUpabs2 = sf2 != 0. ? sf2  + (uncscale*dsfUp2/sf2) : 0. ; 
        double sfDownabs2 = sf2 != 0. ? sf2  + (uncscale*dsfDown2/sf2) : 0. ; 
        if ( sfUpabs < 0 ) sfUpabs = 0 ;  
        if ( sfDownabs < 0 ) sfDownabs = 0 ;  

        if ( boost::math::isnan(sf) || boost::math::isnan(sfUp) || boost::math::isnan(sfUpabs)) {
          std::cout << " sf = " << sf << " is nan = " << boost::math::isnan(sf) << std::endl ; 
          std::cout << " sf2 = " << sf2 << " is nan = " << boost::math::isnan(sf) << std::endl ; 
          std::cout << " pt = " << pt << " is nan = " << boost::math::isnan(sf) << std::endl ; 
          std::cout << " eta = " << eta << " is nan = " << boost::math::isnan(sf) << std::endl ; 
          std::cout << " eff = " << eff << " is nan = " << boost::math::isnan(sf) << std::endl ; 
          std::cout << " eff2 = " << eff2 << " is nan = " << boost::math::isnan(sf) << std::endl ; 
          std::cout << " sfUp = " << sfUp << " is nan = " << boost::math::isnan(sfUp) << std::endl ; 
          std::cout << " dsfUp = " << dsfUp << " is nan = " << boost::math::isnan(dsfUp) << std::endl ; 
          std::cout << " sfUpabs = " << sfUpabs << " is nan = " << boost::math::isnan(sfUpabs) << std::endl ; 
        }

        double jetcsv = jetcsvs.at(idx.first) ;
        //std::cout << "CSV is: " << jetcsv << std::endl;
        //std::cout << "pT is: " << pt << std::endl;
        //std::cout << "eta is: " << eta << std::endl;
        //std::cout << "sf is: " << sf << std::endl;
        //std::cout << "eff is: " << eff << std::endl;
        //std::cout << "sf2 is: " << sf2 << std::endl;
        //std::cout << "eff2 is: " << eff2 << std::endl;
        //std::cout << std::endl;
 
        if ( jetcsv >= csvMin2 ) { 
          btagsf *= sf2 ; 
          //std::cout << "btagsf is: " << btagsf << std::endl;
          //std::cout << std::endl;
          //// Get uncertainties 
          if ( fl == 0 || fl == 1) {
            btagsf_bcUp *= sfUpabs2 ; 
            btagsf_bcDown *= sfDownabs2 ; 
          }
          else {
            btagsf_lUp *= sfUpabs2 ; 
            btagsf_lDown *= sfDownabs2 ; 
          }

        }
        else if (jetcsv < csvMin2 && jetcsv >= csvMin) {
          btagsf *= eff < 1 ? (std::min(1.,eff*sf) - std::min(1.,eff2*sf2))/(eff - eff2) : 0 ; 
          //std::cout << "btagsf is: " << btagsf << std::endl;
          //std::cout << std::endl;
         
          //// Get uncertainties 
          if ( fl == 0 || fl == 1) {
            btagsf_bcUp *= eff < 1 ? (std::min(1.,eff*sfUpabs) - std::min(1.,eff2*sfUpabs2))/(eff - eff2) : 0 ;
            btagsf_bcDown *= eff < 1 ? (std::min(1.,eff*sfDownabs) - std::min(1.,eff2*sfDownabs2))/(eff - eff2) : 0 ;
          }
          else {
            btagsf_lUp *= eff < 1 ? (std::min(1.,eff*sfUpabs) - std::min(1.,eff2*sfUpabs2))/(eff - eff2) : 0 ;
            btagsf_lDown *= eff < 1 ? (std::min(1.,eff*sfDownabs) - std::min(1.,eff2*sfDownabs2))/(eff - eff2) : 0 ;
          }

        }
        else if (jetcsv < csvMin) {
          btagsf *= eff < 1 ? (1 - std::min(1.,eff*sf))/(1 - eff) : 0 ; 
          //std::cout << "btagsf is: " << btagsf << std::endl;
          //std::cout << std::endl;

          //// Get uncertainties 
          if ( fl == 0 || fl == 1) {
            btagsf_bcUp *= eff < 1 ? (1 - std::min(1.,eff*sfUpabs))/(1 - eff) : 0 ;
            btagsf_bcDown *= eff < 1 ? (1 - std::min(1.,eff*sfDownabs))/(1 - eff) : 0 ;
          }
          else {
            btagsf_lUp *= eff < 1 ? (1 - std::min(1.,eff*sfUpabs))/(1 - eff) : 0 ;
            btagsf_lDown *= eff < 1 ? (1 - std::min(1.,eff*sfDownabs))/(1 - eff) : 0 ;
          }

        }
      }

      return ; 
    }

    double getBTagEff (double pt, double eta, double jetFl) {

      std::unique_ptr<TFile>f_effmap = std::unique_ptr<TFile>(new TFile(btageffmap_.c_str())) ;
      std::unique_ptr<TH2>h2_btageffmap_b = std::unique_ptr<TH2>(dynamic_cast<TH2*>(f_effmap->Get("eff_b"))) ; 
      std::unique_ptr<TH2>h2_btageffmap_c = std::unique_ptr<TH2>(dynamic_cast<TH2*>(f_effmap->Get("eff_c"))) ; 
      std::unique_ptr<TH2>h2_btageffmap_l = std::unique_ptr<TH2>(dynamic_cast<TH2*>(f_effmap->Get("eff_l"))) ; 

      double eff(1) ; 
      jetFl = abs(jetFl) ; 
      eta = abs(eta);
      if (jetFl == 5 && pt >= bfl_ptMin_)  {
        int binpt = h2_btageffmap_b->GetXaxis()->FindBin(pt);
        int bineta = h2_btageffmap_b->GetYaxis()->FindBin(eta);
        eff = h2_btageffmap_b->GetBinContent(binpt, bineta) ; 
      }
      else if (jetFl == 4 && pt >= cfl_ptMin_)  {
        int binpt = h2_btageffmap_c->GetXaxis()->FindBin(pt);
        int bineta = h2_btageffmap_c->GetYaxis()->FindBin(eta);
        eff = h2_btageffmap_c->GetBinContent(binpt, bineta) ; 
      }
      else if (jetFl == 0 && pt >= lfl_ptMin_)  {
        int binpt = h2_btageffmap_l->GetXaxis()->FindBin(pt);
        int bineta = h2_btageffmap_l->GetYaxis()->FindBin(eta);
        eff = h2_btageffmap_l->GetBinContent(binpt, bineta) ; 
      }
      else{  eff = 0.;
          }

      return eff ;

    }
    double getBTagEff2nd (double pt, double eta, double jetFl) {

      std::unique_ptr<TFile>f_effmap = std::unique_ptr<TFile>(new TFile(btageffmap2nd_.c_str())) ;
      std::unique_ptr<TH2>h2_btageffmap2nd_b = std::unique_ptr<TH2>(dynamic_cast<TH2*>(f_effmap->Get("eff_b"))) ; 
      std::unique_ptr<TH2>h2_btageffmap2nd_c = std::unique_ptr<TH2>(dynamic_cast<TH2*>(f_effmap->Get("eff_c"))) ; 
      std::unique_ptr<TH2>h2_btageffmap2nd_l = std::unique_ptr<TH2>(dynamic_cast<TH2*>(f_effmap->Get("eff_l"))) ; 

      double eff(1) ; 
      jetFl = abs(jetFl) ; 
      eta = abs(eta);
      if (jetFl == 5 && pt >= bfl_ptMin_)  {
        int binpt = h2_btageffmap2nd_b->GetXaxis()->FindBin(pt);
        int bineta = h2_btageffmap2nd_b->GetYaxis()->FindBin(eta);
        eff = h2_btageffmap2nd_b->GetBinContent(binpt, bineta) ; 
      }
      else if (jetFl == 4 && pt >= cfl_ptMin_)  {
        int binpt = h2_btageffmap2nd_c->GetXaxis()->FindBin(pt);
        int bineta = h2_btageffmap2nd_c->GetYaxis()->FindBin(eta);
        eff = h2_btageffmap2nd_c->GetBinContent(binpt, bineta) ; 
      }
      else if (jetFl == 0 && pt >= lfl_ptMin_)  {
        int binpt = h2_btageffmap2nd_l->GetXaxis()->FindBin(pt);
        int bineta = h2_btageffmap2nd_l->GetYaxis()->FindBin(eta);
        eff = h2_btageffmap2nd_l->GetBinContent(binpt, bineta) ; 
      }
      else{  eff = 0.;
          }

      return eff ;

    }

    static double getBTagEff_CSVv2M (double pt, double jetFl) {

      double eff(1) ; 
      jetFl = abs(jetFl) ; 
      if ( pt > 0 ) {
        if (jetFl == 5) eff = 0.6 ; 
        else if (jetFl == 4) eff = 0.2 ; 
        else if (jetFl == 0) eff = 0.01 ; 
        else eff = 0.; 
      }
      else eff = 0.; 
      return eff ;
    }

    static double getBTagSF_CSVv2L (double jetPt, double jetFl, double errbc, double errl) {
      double btagsf(1);

      double pt(jetPt) ; 
      jetFl = abs(jetFl) ; 
      if ((jetFl == 5 || jetFl == 4)) {
        if (jetPt < 30.) {
          pt = 30; 
          errbc *= 2;
        }
        if (jetPt > 670.) {
          pt = 670. ;
          errbc *= 2;
        }
      }
      else {
        if (jetPt < 20.) {
          pt= 20; 
          errl *= 2;
        }
        if (jetPt > 1000.) {
          pt = 1000. ;
          errl *= 2; 
        }
      }

      TF1* SFb_CSVv2L = new TF1("SFb_CSVv2L", "0.908299+(2.70877e-06*(log(x+370.144)*(log(x+370.144)*(3-(-(104.614*log(x+370.144)))))))",0,2000) ;

      std::map<std::pair<double,double>,double> SFb_CSVv2L_err = { 
        {std::make_pair(30,50)  , 0.022327613085508347}, 
        {std::make_pair(50,70)  , 0.015330483205616474}, 
        {std::make_pair(70,100) , 0.024493992328643799}, 
        {std::make_pair(100,140), 0.020933238789439201}, 
        {std::make_pair(140,200), 0.029219608753919601}, 
        {std::make_pair(200,300), 0.039571482688188553}, 
        {std::make_pair(300,670), 0.047329759448766708}  
      } ;

      TF1* SFc_CSVv2L = new TF1("SFc_CSVv2L", "0.908299+(2.70877e-06*(log(x+370.144)*(log(x+370.144)*(3-(-(104.614*log(x+370.144)))))))",0,2000) ;

      std::map<std::pair<double,double>,double> SFc_CSVv2L_err = { 
        {std::make_pair(30,50)  , 0.044655226171016693}, 
        {std::make_pair(50,70)  , 0.030660966411232948}, 
        {std::make_pair(70,100) , 0.048987984657287598}, 
        {std::make_pair(100,140), 0.041866477578878403}, 
        {std::make_pair(140,200), 0.058439217507839203}, 
        {std::make_pair(200,300), 0.079142965376377106}, 
        {std::make_pair(300,670), 0.094659518897533417}  
      } ;

      TF1* SFl_CSVv2L = new TF1("SFl_CSVv2L", "((1.07278+(0.000535714*x))+(-1.14886e-06*(x*x)))+(7.0636e-10*(x*(x*x)))",0,2000) ;
      TF1* SFl_CSVv2L_errUp = new TF1("SFl_CSVv2L_errUp", "((1.12921+(0.000804962*x))+(-1.87332e-06*(x*x)))+(1.18864e-09*(x*(x*x)))",0,2000) ; 
      TF1* SFl_CSVv2L_errDown = new TF1("SFl_CSVv2L_errDown", "((1.01637+(0.000265653*x))+(-4.22531e-07*(x*x)))+(2.23396e-10*(x*(x*x)))",0,2000) ;

      double errb(0), errc(0) ; 
      for ( auto const& ent : SFb_CSVv2L_err ) if (pt >= (ent.first).first && pt < (ent.first).second ) errb = ent.second ; 
      for ( auto const& ent : SFc_CSVv2L_err ) if (pt >= (ent.first).first && pt < (ent.first).second ) errc = ent.second ; 

      if (jetFl == 5) btagsf = SFb_CSVv2L->Eval(pt) + errb*errbc; 
      else if (jetFl == 4) btagsf = SFc_CSVv2L->Eval(pt) + errc*errbc; 
      else if (jetFl == 0) btagsf = (
          SFl_CSVv2L->Eval(pt)*(1 - abs(errl)) + 
          (SFl_CSVv2L_errUp->Eval(pt)*abs(errl)*(1+errl)/2) + 
          (SFl_CSVv2L_errDown->Eval(pt)*abs(errl)*(1-errl)/2) ) ;  

      return btagsf ; 
    }

};
#endif 
