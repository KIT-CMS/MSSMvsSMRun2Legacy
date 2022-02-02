#include "CombineHarvester/CombinePdfs/interface/MorphFunctions.h"
#include "CombineHarvester/CombinePdfs/interface/CMSHistFuncFactory.h"
#include "CombineHarvester/CombineTools/interface/Algorithm.h"
#include "CombineHarvester/CombineTools/interface/AutoRebin.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"
#include "CombineHarvester/CombineTools/interface/CardWriter.h"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/MSSMvsSMRun2Legacy/interface/HttSystematics_MSSMvsSMRun2.h"
#include "CombineHarvester/MSSMvsSMRun2Legacy/interface/BinomialBinByBin.h"
#include "CombineHarvester/MSSMvsSMRun2Legacy/interface/dout_tools.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TF1.h"
#include "TH2.h"
#include "boost/algorithm/string/predicate.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/regex.hpp"
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <math.h>

using namespace std;
using boost::starts_with;
namespace po = boost::program_options;

template<typename T>
void update_vector_by_byparser(T& parameter, const T& parser, const string name="") {
    if (parser.size() > 0)
    {
        const std::string lower_str = boost::algorithm::to_lower_copy(parser[0]);
        doutnonl("WARNING: The", name, "are set manually:");

        if (parser.size() == 1 && (lower_str == "none" || lower_str == "null" || lower_str == "pass"))
          parameter = {};
        else
          parameter = parser;

        dprintVector(parameter);
    }
}


void ConvertShapesToLnN (ch::CombineHarvester& cb, string name) {
  auto cb_syst = cb.cp().syst_name({name});
  cb_syst.ForEachSyst([&](ch::Systematic *syst) {
    if (syst->type() == "shape") {
      std::cout << "Converting systematic " << syst->name() << " for process " << syst->process() << " in bin " << syst->bin() << " to lnN." <<std::endl;
      syst->set_type("lnN");
      return;
    }
  });
}

int main(int argc, char **argv) {
  typedef vector<string> VString;
  typedef vector<pair<int, string>> Categories;
  using ch::syst::bin_id;
  using ch::JoinStr;

  // Define program options
  string output_folder = "output_MSSMvsSM_Run2";
  string base_path = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/MSSMvsSMRun2Legacy/shapes/variables/";
  string category = "mt_nobtag_lowmsv";
  string variable = "mjj";

  bool real_data = true;
  bool verbose = false;
  bool use_automc = true;
  bool manual_rebinning = false;
  bool use_mc = false;
  std::string mode = "control-plots";
  std::string chan;

  int era = 2016; // 2016, 2017 or 2018
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()
      ("base-path,base_path", po::value<string>(&base_path)->default_value(base_path), "inputs, expected to contain a subdirectory <era>/<channel>")
      ("category", po::value<string>(&category)->default_value(category))
      ("variable", po::value<string>(&variable)->default_value(variable))
      ("channel", po::value<string>(&chan)->default_value(chan))
      ("real_data", po::value<bool>(&real_data)->default_value(real_data))
      ("use_automc", po::value<bool>(&use_automc)->default_value(use_automc))
      ("manual_rebinning", po::value<bool>(&manual_rebinning)->default_value(manual_rebinning))
      ("verbose", po::value<bool>(&verbose)->default_value(verbose))
      ("output_folder", po::value<string>(&output_folder)->default_value(output_folder))
      ("use_mc", po::value<bool>(&use_mc)->default_value(use_mc))
      ("era", po::value<int>(&era)->default_value(era))
      ("mode", po::value<string>(&mode)->default_value(mode))
      ("help", "produce help message");
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);
  if (vm.count("help"))
  {
      cout << config << "\n";
      return 0;
  }
  // Check if mode has a sensible value. Allow for modes control-plots, categorisation-plots
  bool found_Mode = false;
  for (std::string avail_mode: {"control-plots", "categorisation-plots"}) {
    if (mode == avail_mode) {
        found_Mode = true;
        break;
    }
  }
  if (!found_Mode) {
      cerr << "Given mode " << mode << " not known. Choose between `control-plots` and `categorisation-plots`." << std::endl;
      return 1;
  }

  // Define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
  std::string era_tag;
  if (era == 2016) era_tag = "2016";
  else if (era == 2017) era_tag = "2017";
  else if (era == 2018) era_tag = "2018";
  else std::runtime_error("Given era is not implemented.");

  if (mode != "categorisation-plots") {
      output_folder = output_folder + "_" + variable;
  }
  std::map<string, string> input_dir;
  if (base_path.back() != '/' ) base_path += "/";
  if (!boost::filesystem::exists(output_folder)) boost::filesystem::create_directories(output_folder);
  input_dir["mt"] = base_path;
  input_dir["et"] = base_path;
  input_dir["tt"] = base_path;
  input_dir["em"] = base_path;


  // Define channels
  VString chns;
  std::cout << chan << std::endl;
  if (chan.find("mt") != std::string::npos)
      chns.push_back("mt");
  if (chan.find("et") != std::string::npos)
      chns.push_back("et");
  if (chan.find("tt") != std::string::npos)
      chns.push_back("tt");
  if (chan.find("em") != std::string::npos)
      chns.push_back("em");

  // Define restriction to the channel defined by '--category' option
  if(category != "all"){
    std::vector<std::string> category_split;
    boost::split(category_split, category, boost::is_any_of("_"));
    chns = {category_split.at(0)};
  }
  doutnonl("Channels:\n\t");
  dprintVector(chns);

  // Define background and signal processes
  map<string, VString> bkg_procs;
  VString bkgs, bkgs_tt, bkgs_em, sm_signals, main_sm_signals, mssm_bbH_signals, hww_signals;

  sm_signals = {"WH125", "ZH125", "ttH125"};
  hww_signals = {"ggHWW125", "qqHWW125", "WHWW125", "ZHWW125"};
  main_sm_signals = {"ggH125", "qqH125"};

  mssm_bbH_signals = {"bbH_500", "bbH_1400"};
  auto all_signals = ch::JoinStr({main_sm_signals, sm_signals, hww_signals, mssm_bbH_signals});
  if (mode == "categorisation-plots") {
      mssm_bbH_signals = {"bbH_400"};
      all_signals = ch::JoinStr({main_sm_signals, mssm_bbH_signals});
  }

  bkgs = {"EMB", "ZL", "TTL", "VVL", "jetFakes"};
  bkgs_tt = {"EMB", "ZL", "TTL", "VVL", "jetFakes", "wFakes"};
  bkgs_em = {"EMB", "W", "QCD", "ZL", "TTL", "VVL"};

  if (use_mc) {
    std::cout << "WARNING: the EMB process is removed from backgrounds" << std::endl;
    // Remove embedded shapes.
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "EMB"), bkgs.end());
    bkgs_tt.erase(std::remove(bkgs_tt.begin(), bkgs_tt.end(), "EMB"), bkgs_tt.end());
    bkgs_em.erase(std::remove(bkgs_em.begin(), bkgs_em.end(), "EMB"), bkgs_em.end());
    VString to_add = {"ZTT", "TTT", "VVT"};
    for (auto proc: to_add) {
        bkgs.push_back(proc);
        bkgs_tt.push_back(proc);
        bkgs_em.push_back(proc);
    }
  }

  bkg_procs["et"] = bkgs;
  bkg_procs["mt"] = bkgs;
  bkg_procs["tt"] = bkgs_tt;
  bkg_procs["em"] = bkgs_em;


  // Define MSSM model-independent mass parameter MH
  RooRealVar MH("MH", "MH", 500., 90., 4000.);
  if (mode == "categorisation-plots") {
      MH.setVal(400.);
  }
  MH.setConstant(true);

  // Define categories
  map<string, Categories> cats;
  // STXS stage 0 categories (optimized on ggH and VBF)
  if (mode == "categorisation-plots") {
      cats["mt"] = {
          {  301, "mt_Nbtag0"},
          {  302, "mt_NbtagGt1"},
      };
      cats["et"] = {
          {  301, "et_Nbtag0"},
          {  302, "et_NbtagGt1"},
      };
      cats["em"] = {
          {  301, "em_Nbtag0"},
          {  302, "em_NbtagGt1"},
      };
  }
  else {
      cats["et"] = {
          {  300, "et_" + variable},
      };
      cats["mt"] = {
          {  300, "mt_" + variable},
      };
      cats["tt"] = {
          {  300, "tt_" + variable},
      };
      cats["em"] = {
          {  300, "em_" + variable},
      };
  }

  // Create combine harverster object
  ch::CombineHarvester cb;
  cb.SetFlag("workspaces-use-clone", true);



  for (auto chn : chns) {
    cb.AddObservations({"*"}, {"htt"}, {era_tag}, {chn}, cats[chn]);
    cb.AddProcesses({"*"}, {"htt"}, {era_tag}, {chn}, bkg_procs[chn], cats[chn], false);
    if (mode == "categorisation-plots") {
        cb.AddProcesses({""}, {"htt"}, {era_tag}, {chn}, main_sm_signals, cats[chn], false);
        cb.AddProcesses({""}, {"htt"}, {era_tag}, {chn}, mssm_bbH_signals, cats[chn], true);
    }
    else {
        cb.AddProcesses({""}, {"htt"}, {era_tag}, {chn}, main_sm_signals, cats[chn], true);
    }
  }

  // Add systematics
  ch::AddMSSMvsSMRun2Systematics(cb, true, !use_mc, true, false, false, era);

  // Define restriction to the desired category
  if(category != "all"){
    cb = cb.bin({category});
  }
  // cb.PrintAll();

  // Extract shapes from input ROOT files
  for (string chn : chns) {
    string input_file_base = input_dir[chn] + "htt_" + category + ".inputs-mssm-vs-sm-Run" + era_tag + ".root";
    if (mode == "categorisation-plots") {
        input_file_base = input_dir[chn] + "htt_" + category + ".inputs-mssm-vs-sm-Run" + era_tag + "-" + variable + ".root";
    }

    cb.cp().channel({chn}).backgrounds().ExtractShapes(
      input_file_base, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");
    if (mode == "categorisation-plots") {
        cb.cp().channel({chn}).process(mssm_bbH_signals).ExtractShapes(
          input_file_base, "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC");
    }
    else {
        cb.cp().channel({chn}).process(main_sm_signals).ExtractShapes(
          input_file_base, "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC");
    }
  }

  // Delete processes (other than mssm signals) with 0 yield
  cb.FilterProcs([&](ch::Process *p) {
    bool null_yield = !(p->rate() > 0.0);
    if (null_yield) {
      std::cout << "[WARNING] Removing process with null yield: \n ";
      std::cout << ch::Process::PrintHeader << *p << "\n";
      cb.FilterSysts([&](ch::Systematic *s) {
        bool remove_syst = (MatchingProcess(*p, *s));
        return remove_syst;
      });
    }
    return null_yield;
  });


  // Modify systematic variations with yield <= 0
  cb.FilterSysts([&](ch::Systematic *s) {
    // For remaining processes: Delete systematics since these result in a bogus norm error in combine for the remaining
    if (s->type() == "shape") {
      if (s->shape_u()->Integral() <= 0.0) {
        std::cout << "[WARNING] Removing systematic with null yield in up shift:" << std::endl;
        std::cout << ch::Systematic::PrintHeader << *s << "\n";
        return true;
      }
      if (s->shape_d()->Integral() <= 0.0) {
        std::cout << "[WARNING] Removing systematic with null yield in down shift:" << std::endl;
        std::cout << ch::Systematic::PrintHeader << *s << "\n";
        return true;
      }
    }
    return false;
  });

  // Convert systematics from shape to logNormal as done in the nominal analysis.
  // turn all JES+JER+met-unclustered uncertainties into lnN for MSSM categories - should check the met recoil uncertainties as well for small processes
  std::vector<std::string> jetmet_systs = {
    "CMS_scale_j_Absolute",
    "CMS_scale_j_BBEC1",
    "CMS_scale_j_EC2",
    "CMS_scale_j_FlavorQCD",
    "CMS_scale_j_HF",
    "CMS_scale_j_RelativeBal",
    "CMS_scale_j_Absolute_2016",
    "CMS_scale_j_Absolute_2017",
    "CMS_scale_j_Absolute_2018",
    "CMS_scale_j_BBEC1_2016",
    "CMS_scale_j_BBEC1_2017",
    "CMS_scale_j_BBEC1_2018",
    "CMS_scale_j_EC2_2016",
    "CMS_scale_j_EC2_2017",
    "CMS_scale_j_EC2_2018",
    "CMS_scale_j_HF_2016",
    "CMS_scale_j_HF_2017",
    "CMS_scale_j_HF_2018",
    "CMS_scale_j_RelativeSample_2016",
    "CMS_scale_j_RelativeSample_2017",
    "CMS_scale_j_RelativeSample_2018",
    "CMS_res_j_2016",
    "CMS_res_j_2017",
    "CMS_res_j_2018",
    "CMS_scale_met_unclustered_2016",
    "CMS_scale_met_unclustered_2017",
    "CMS_scale_met_unclustered_2018",
  };

  // Convert all JES ,JER, and MET uncertainties to lnN except for the ttbar uncertainties in the em, et and mt channels
  // These uncertainties affect MET for the ttbar and diboson so we need to include them as shapes (diboson is small enough to be converted to lnN, and is ttbar in the tt channel)
  // convert all processes except ttbar
  for(auto u : jetmet_systs) ConvertShapesToLnN (cb.cp().process({"TTL","TTT"},false), u);
  // also convert ttbar in the tt channel
  for(auto u : jetmet_systs) ConvertShapesToLnN (cb.cp().channel({"tt"}).process({"TTL","TTT"}), u);
  // also convert ttbar in the nobtag categories, but only when fitting 1D variables

  // some FF unc1 systematics for the tt channel only affect the normalisations so can be converted to lnN:
  for (string y : {"2016","2017","2018"}) {
    ConvertShapesToLnN(cb.cp().bin_id({32,35,135,235,335,435,132,232,332,432}), "CMS_ff_total_qcd_stat_dR_unc1_tt_"+y);
    ConvertShapesToLnN(cb.cp().bin_id({32,35,135,235,335,435,132,232,332,432}), "CMS_ff_total_qcd_stat_pt_unc1_tt_"+y);
    ConvertShapesToLnN(cb.cp().bin_id({35,135,235,335,435}),    "CMS_ff_total_qcd_stat_njet1_jet_pt_low_unc1_tt_"+y);
    ConvertShapesToLnN(cb.cp().bin_id({35,135,235,335,435}),    "CMS_ff_total_qcd_stat_njet1_jet_pt_med_unc1_tt_"+y);
    ConvertShapesToLnN(cb.cp().bin_id({35,135,235,335,435}),    "CMS_ff_total_qcd_stat_njet1_jet_pt_high_unc1_tt_"+y);
  }

  // Try to rename energy scale in embedded samples to see check decorrelation
  // for (string era : {"2016","2017","2018"}) {
  //     cb.cp().channel({"et", "em"}).era({era}).RenameSystematic(cb, "CMS_scale_e_emb", "CMS_scale_e_emb_"+era);
  // }
  // cb.cp().channel({"et", "em"}).RenameSystematic(cb,"scale_e_toMC", "scale_e_emb_toMC");
  // rename some fake factor systematics so that they are decorrelated between categories to match how closure corrections are measured
  for (string y : {"2016","2017","2018"}) {

    for (string u : {"unc1", "unc2"}) {

      cb.cp().bin_id({32,301}).channel({"tt"}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_dR_"+u+"_tt_"+y,"CMS_ff_total_qcd_stat_dR_"+u+"_tt_Nbtag0_"+y);
      cb.cp().bin_id({35,302}).channel({"tt"}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_dR_"+u+"_tt_"+y,"CMS_ff_total_qcd_stat_dR_"+u+"_tt_NbtagGt1_"+y);

      cb.cp().bin_id({32,301}).channel({"tt"}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_pt_"+u+"_tt_"+y,"CMS_ff_total_qcd_stat_pt_"+u+"_tt_Nbtag0_"+y);
      cb.cp().bin_id({35,302}).channel({"tt"}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_pt_"+u+"_tt_"+y,"CMS_ff_total_qcd_stat_pt_"+u+"_tt_NbtagGt1_"+y);

      for (string c : {"mt","et"}) {
        // TODO: What do we do here for the categorisation split?
        // cb.cp().bin_id({32,132,232,332,432}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_"+y,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_Nbtag0_MTLt40_"+y);
        // cb.cp().bin_id({33,133,233,333,433}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_"+y,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_Nbtag0_MT40To70_"+y);
        // cb.cp().bin_id({35,135,235,335,435}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_"+y,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_NbtagGt1_MTLt40_"+y);
        // cb.cp().bin_id({36,136,236,336,436}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_"+y,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_NbtagGt1_MT40To70_"+y);

        // cb.cp().bin_id({32,35,135,235,335,435,132,232,332,432}).channel({c}).RenameSystematic(cb,"CMS_ff_total_ttbar_stat_l_pt_"+u+"_"+c+"_"+y,"CMS_ff_total_ttbar_stat_l_pt_"+u+"_"+c+"_MTLt40_"+y);
        // cb.cp().bin_id({33,36,136,236,336,436,133,233,333,433}).channel({c}).RenameSystematic(cb,"CMS_ff_total_ttbar_stat_l_pt_"+u+"_"+c+"_"+y,"CMS_ff_total_ttbar_stat_l_pt_"+u+"_"+c+"_MT40To70_"+y);

        cb.cp().bin_id({32,33,301}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_os_"+u+"_"+c+"_"+y,"CMS_ff_total_qcd_stat_os_"+u+"_"+c+"_Nbtag0_"+y);
        cb.cp().bin_id({35,36,302}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_os_"+u+"_"+c+"_"+y,"CMS_ff_total_qcd_stat_os_"+u+"_"+c+"_NbtagGt1_"+y);

      }
    }
  }

  // the following code is used to decorrelate wjets and qcd systematics by category
  // the btag and nobtag are always decorrelated
  // we also decorrelate the wjets by loose and tight mT since we are extrapolating to different mT regions
  // in cases where these uncertainties were not derived seperatly for Nbjets>0 (usually due to limited stats) we double the uncertainty
  // note this doubling is already done in the FF workspaces for the wjets_syst_extrap so we don't need to do it again here
  // we decorrelate the wjets extrapolation uncertainties by category

  // TODO: What do we do here for the categorisation split?
  // for (string y : {"2016","2017","2018"}) {
  //   for (string c : {"mt","et"}) {
  //     cb.cp().bin_id({32,132,232,332,432}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_extrap_"+c+"_"+y,"CMS_ff_total_wjets_syst_extrap_"+c+"_Nbtag0_MTLt40_"+y);
  //     cb.cp().bin_id({33,133,233,333,433}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_extrap_"+c+"_"+y,"CMS_ff_total_wjets_syst_extrap_"+c+"_Nbtag0_MT40To70_"+y);
  //     cb.cp().bin_id({35,135,235,335,435}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_extrap_"+c+"_"+y,"CMS_ff_total_wjets_syst_extrap_"+c+"_NbtagGt1_MTLt40_"+y);
  //     cb.cp().bin_id({36,136,236,336,436}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_extrap_"+c+"_"+y,"CMS_ff_total_wjets_syst_extrap_"+c+"_NbtagGt1_MT40To70_"+y);

  //     cb.cp().bin_id({32,132,232,332,432}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_"+c+"_"+y,"CMS_ff_total_wjets_syst_"+c+"_Nbtag0_MTLt40_"+y);
  //     cb.cp().bin_id({33,133,233,333,433}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_"+c+"_"+y,"CMS_ff_total_wjets_syst_"+c+"_Nbtag0_MT40To70_"+y);
  //     cb.cp().bin_id({35,135,235,335,435}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_"+c+"_"+y,"CMS_ff_total_wjets_syst_"+c+"_NbtagGt1_MTLt40_"+y);
  //     cb.cp().bin_id({36,136,236,336,436}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_"+c+"_"+y,"CMS_ff_total_wjets_syst_"+c+"_NbtagGt1_MT40To70_"+y);
  //   }
  // }
  for (string y : {"2016","2017","2018"}) {
    string c = "tt";
    cb.cp().bin_id({32,301}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_"+c+"_"+y,"CMS_ff_total_wjets_syst_"+c+"_Nbtag0_"+y);
    cb.cp().bin_id({35,302}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_"+c+"_"+y,"CMS_ff_total_wjets_syst_"+c+"_NbtagGt1_"+y);

    cb.cp().bin_id({32,301}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_syst_dr_closure_"+c+"_"+y,"CMS_ff_total_qcd_syst_dr_closure_"+c+"_Nbtag0_"+y);
    cb.cp().bin_id({35,302}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_syst_dr_closure_"+c+"_"+y,"CMS_ff_total_qcd_syst_dr_closure_"+c+"_NbtagGt1_"+y);

    // scale Wjets uncertainty by 2 for tt channels in btag category.
    cb.cp().bin_id({35,302}).channel({c}).syst_name({"CMS_ff_total_wjets_syst_"+c+"_NbtagGt1_"+y}).ForEachSyst([&](ch::Systematic *syst) {
      if (syst->type().find("shape") != std::string::npos) {
        syst->set_scale(syst->scale() * 2.);
      }
      if (syst->type().find("lnN") != std::string::npos) {
        syst->set_value_u((syst->value_u() - 1.) * 2. + 1.);
        if (syst->asymm()){
          syst->set_value_d((syst->value_d() - 1.) * 2. + 1.);
        }
      }
    });
  }

  // we decorrelate the qcd extrapolation systematics by NBtag

  for (string y : {"2016","2017","2018"}) {
    for (string c : {"mt","et","tt"}) {
      cb.cp().bin_id({32,33,301}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_syst_"+c+"_"+y,"CMS_ff_total_qcd_syst_"+c+"_Nbtag0_"+y);
      cb.cp().bin_id({35,36,302}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_syst_"+c+"_"+y,"CMS_ff_total_qcd_syst_"+c+"_NbtagGt1_"+y);
      if(c != "tt") {
        cb.cp().bin_id({32,33,301}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_syst_iso_"+c+"_"+y,"CMS_ff_total_qcd_syst_iso_"+c+"_Nbtag0_"+y);
        cb.cp().bin_id({35,36,302}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_syst_iso_"+c+"_"+y,"CMS_ff_total_qcd_syst_iso_"+c+"_NbtagGt1_"+y);
        // scale QCD extrapolation uncertainties by 2 for et and mt channels in btag category.
        cb.cp().bin_id({35,36,302}).channel({c}).syst_name({"CMS_ff_total_qcd_syst_"+c+"_NbtagGt1_"+y,"CMS_ff_total_qcd_syst_iso_"+c+"_NbtagGt1_"+y}).ForEachSyst([&](ch::Systematic *syst) {
          if (syst->type().find("shape") != std::string::npos) {
            syst->set_scale(syst->scale() * 2.);
          }
          if (syst->type().find("lnN") != std::string::npos) {
            syst->set_value_u((syst->value_u() - 1.) * 2. + 1.);
            if (syst->asymm()){
              syst->set_value_d((syst->value_d() - 1.) * 2. + 1.);
            }
          }
        });
      }
    }
  }

  //  ConvertShapesToLnN(cb.cp().channel({"em"}), "subtrMC");
  ConvertShapesToLnN(cb.cp().channel({"em"}), "htt_em_QFlip");
  //ConvertShapesToLnN(cb.cp().channel({"em"}), "htt_em_MuToEFakes_SS"); // very small so removed
  ConvertShapesToLnN(cb.cp().channel({"em"}).bin_id({33,133},false), "htt_em_MuToEFakes_OS"); // shape only in cateogory with large ZL
  ConvertShapesToLnN(cb.cp().channel({"em"}).process({"TTL","VVL","QCD"}), "htt_em_JToMuFakes");
  ConvertShapesToLnN(cb.cp().channel({"em"}).process({"ZL"},false), "htt_em_JToEFakes"); // To lnN for all processes except ZL in largest category
  ConvertShapesToLnN(cb.cp().channel({"em"}).bin_id({33,133},false), "htt_em_JToEFakes");
  ConvertShapesToLnN(cb.cp().channel({"em"}), "CMS_htt_qcd_iso");
  for (std::string y: {"2016", "2017", "2018"}) {
    for (std::string jetbin: {"0jet","1jet","2jet"}) {
      for (std::string shapebin: {"rate","shape","shape2"}) {
	ConvertShapesToLnN(cb.cp().channel({"em"}), "CMS_htt_qcd_"+jetbin+"_"+shapebin+"_"+y);
      }
    }
  }

  cb.cp().channel({"em"}).era({"2016"}).RenameSystematic(cb, "htt_em_QFlip", "htt_em_QFlip_2016"); // this uncertainty is systematics limited but correction for 2016 has different magnitude to 2017 and 2018 so decorrelate 

  for (std::string y: {"2016", "2017", "2018"}) {
      cb.cp().channel({"em"}).era({y}).RenameSystematic(cb, "CMS_htt_qcd_iso", "CMS_htt_qcd_iso_"+y);
      cb.cp().channel({"em"}).era({y}).RenameSystematic(cb, "htt_em_MuToEFakes_OS", "htt_em_MuToEFakes_OS_"+y);
      //cb.cp().channel({"em"}).era({y}).RenameSystematic(cb, "htt_em_MuToEFakes_SS", "htt_em_MuToEFakes_SS_"+y); // very small so removed
//      cb.cp().bin_id({32,33,34,132,232,332,432,133,233,333,433,134,234,334,434}).channel({"em"}).era({y}).RenameSystematic(cb, "subtrMC", "subtrMC_lowttbar_"+y);
//      cb.cp().bin_id({2,35,135,235,335,435,36,136,236,336,436,37,137,237,337,437}).channel({"em"}).era({y}).RenameSystematic(cb, "subtrMC", "subtrMC_highttbar_"+y);
  }

//  cb.cp().bin_id({32,33,34,132,232,332,432,133,233,333,433,134,234,334,434}).channel({"em"}).RenameSystematic(cb, "subtrMC", "subtrMC_lowttbar");
//  cb.cp().bin_id({2,35,135,235,335,435,36,136,236,336,436,37,137,237,337,437}).channel({"em"}).RenameSystematic(cb, "subtrMC", "subtrMC_highttbar");

  std::vector<std::string> met_uncerts = {
    "CMS_htt_boson_scale_met_2016",
    "CMS_htt_boson_res_met_2016",
    "CMS_htt_boson_scale_met_2017",
    "CMS_htt_boson_res_met_2017",
    "CMS_htt_boson_scale_met_2018",
    "CMS_htt_boson_res_met_2018",
  };

  for(auto u : met_uncerts) ConvertShapesToLnN (cb.cp().process({"ZTT"}, false), u);

  // At this point we can fix the negative bins for the remaining processes
  std::cout << "[INFO] Fixing negative bins.\n";
  cb.ForEachProc([](ch::Process *p) {
    if (ch::HasNegativeBins(p->shape())) {
      auto newhist = p->ClonedShape();
      ch::ZeroNegativeBins(newhist.get());
      p->set_shape(std::move(newhist), false);
    }
  });

  cb.ForEachSyst([](ch::Systematic *s) {
    if (s->type().find("shape") == std::string::npos)
      return;
    if (ch::HasNegativeBins(s->shape_u()) ||
        ch::HasNegativeBins(s->shape_d())) {
      auto newhist_u = s->ClonedShapeU();
      auto newhist_d = s->ClonedShapeD();
      ch::ZeroNegativeBins(newhist_u.get());
      ch::ZeroNegativeBins(newhist_d.get());
      s->set_shapes(std::move(newhist_u), std::move(newhist_d), nullptr);
    }
  });


  // Manual rebinning, if needed
  if(manual_rebinning){
    std::map<std::string, std::vector<double> > binning;
    binning["DiTauDeltaR_et_nobtag_lowmsv"] = {0.5, 2.5, 4.0};
    binning["DiTauDeltaR_tt_nobtag_lowmsv"] = {0.5, 2.5, 3.2, 4.0};
    binning["pZetaPuppiMissVis_em_inclusive"] = {-100.0,-35.0,-10.0,30.0,100.0};
    binning["jdeta_tt_nobtag_lowmsv"] = {0.0, 4.0, 8.0};
    binning["pt_tt_puppi_em_nobtag_lowmsv"] = {0.0, 10.0, 40.0, 120.0, 200.0, 300.0};
    binning["pt_tt_puppi_et_nobtag_lowmsv"] = {0.0, 120.0, 200.0, 300.0};
    binning["pt_tt_puppi_mt_nobtag_lowmsv"] = {0.0, 120.0, 200.0, 300.0};

    for(auto b : cb.cp().bin_set())
    {
      std::string var_bin = variable + "_" + b;
      if(binning.find(var_bin) != binning.end())
      {
          std::cout << "Rebinning by hand for variable, bin: " << var_bin <<  std::endl;
          cb.cp().bin({b}).VariableRebin(binning[var_bin]);
      }
    }
  }

  // This function modifies every entry to have a standardised bin name of
  // the form: {analysis}_{channel}_{bin_id}_{era}
  ch::SetStandardBinNames(cb, "$ANALYSIS_$CHANNEL_$BINID_$ERA");
  ch::CombineHarvester cb_obs = cb.deep().backgrounds();

  // Adding bin-by-bin uncertainties
  if (use_automc) {
    std::cout << "[INFO] Adding bin-by-bin uncertainties.\n";
    cb.SetAutoMCStats(cb, 0.0);
  }

  // Write out datacards. Naming convention important for rest of workflow. We
  // make one directory per chn-cat, one per chn and cmb. In this code we only
  // store the individual datacards for each directory to be combined later.
  ch::CardWriter writer(output_folder + "/" + era_tag + "/$BIN/$BIN.txt",
                        output_folder + "/" + era_tag + "/$BIN/common/$BIN_input.root");

  // We're not using mass as an identifier - which we need to tell the
  // CardWriter
  // otherwise it will see "*" as the mass value for every object and skip it
  writer.SetWildcardMasses({});

  // Set verbosity
  if (verbose)
    writer.SetVerbosity(1);

  writer.WriteCards("", cb);

  ch::CardWriter writer_restore(output_folder + "/" + era_tag + "/restore_binning/$BIN/$BIN.txt",
                                output_folder + "/" + era_tag + "/restore_binning/$BIN/common/$BIN_input.root");
  if (mode == "categorisation-plots") {
      writer_restore = ch::CardWriter(output_folder + "/restore_binning/$BIN.txt",
                                      output_folder + "/restore_binning/common/$BIN_input.root");
  }

  // We're not using mass as an identifier - which we need to tell the
  // CardWriter
  // otherwise it will see "*" as the mass value for every object and skip it
  writer_restore.SetWildcardMasses({});

  // Set verbosity
  if (verbose)
    writer_restore.SetVerbosity(1);

  writer_restore.WriteCards("", cb_obs);

  if (verbose)
    cb.PrintAll();

  std::cout << "[INFO] Done producing datacards.\n";
}
