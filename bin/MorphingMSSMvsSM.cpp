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
#include "boost/algorithm/string.hpp"
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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

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

std::vector<double> binning_from_map(std::map<unsigned int, std::vector<double>> binning_map) {
    std::vector <double> binning = {};
    if(binning_map.begin() != binning_map.end())
    {
        for (auto it = binning_map.begin(); it != binning_map.end(); ++it)
            for (double first = it->second[0]; first < it->second[1]; first+=it->second[2])
            {
                binning.push_back(first);
            }
        auto last = --binning_map.end();
        binning.push_back(last->second[1]);
    }
    return binning;
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
  using ch::syst::bin;
  using ch::syst::mass;
  using ch::JoinStr;
  using ch::syst::SystMap;

  // Define program options
  string output_folder = "output_MSSMvsSM_Run2";
  string base_path = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/MSSMvsSMRun2Legacy/shapes/";
  string sm_gg_fractions = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root";
  string sm_predictions = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_predictions_13TeV.json";
  string chan = "mt";
  string category = "mt_nobtag_lowmsv_0jet_tightmt";
  string variable = "m_sv_puppi";
  string non_morphed_mass = "700";

  bool do_morph = true;
  bool auto_rebin = false;
  bool manual_rebin = true;
  bool real_data = false;
  bool verbose = true;
  bool use_automc = true;
  bool mva(false), no_emb(false);
  bool sm = false;
  bool split_sm_signal_cat = false;
  bool rebin_sm = true;
  bool no_shape_systs = false;
  bool enable_bsm_lowmass = false;
  bool lowmass = false;
  bool prop_plot = false;

  vector<string> mass_susy_ggH({}), mass_susy_qqH({}), parser_bkgs({}), parser_bkgs_em({}), parser_sm_signals({}), parser_main_sm_signals({});

  string analysis = "bsm-model-indep"; // "sm",  "bsm-model-indep", "bsm-model-dep-full", "bsm-model-dep-additional"
  std::vector<string> analysis_choices = {"sm", "bsm-model-indep", "bsm-model-dep-full", "bsm-model-dep-additional"};
  string sub_analysis = "sm-like-light"; // for analysis = "bsm-model-dep-{full,additional}": "sm-like-light", "sm-like-heavy", "cpv"
  std::vector<string> sub_analysis_choices = {"sm-like-light", "sm-like-heavy", "cpv"};
  string hSM_treatment = "hSM-in-bg"; // for analysis = "bsm-model-indep" and = "bsm-model-dep-full" : "hSM-in-bg", "no-hSM-in-bg"; case with analysis = "bsm-model-dep-additional": "hSM-in-bg"
  std::vector<string> hSM_treatment_choices = {"hSM-in-bg", "no-hSM-in-bg"};
  string sm_like_hists = "sm125"; // used in analysis = "bsm-model-dep-full": "sm125", "bsm"
  std::vector<string> sm_like_hists_choices = {"sm125", "bsm"};
  string categorization = "classic"; // "with-sm-ml", "sm-ml-only", "classic"
  std::vector<string> categorization_choices = {"with-sm-ml", "sm-ml-only", "classic", "lowmass"};

  int era = 2016; // 2016, 2017 or 2018
  std::vector<int> era_choices = {2016, 2017, 2018};

  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()
      ("base-path,base_path", po::value<string>(&base_path)->default_value(base_path), "inputs, expected to contain a subdirectory <era>/<channel>")
      ("sm_gg_fractions", po::value<string>(&sm_gg_fractions)->default_value(sm_gg_fractions))
      ("sm_predictions", po::value<string>(&sm_predictions)->default_value(sm_predictions))
      ("channel", po::value<string>(&chan)->default_value(chan), "single channel to process")
      ("category", po::value<string>(&category)->default_value(category))
      ("variable", po::value<string>(&variable)->default_value(variable))
      ("non-morphed-mass", po::value<string>(&non_morphed_mass)->default_value(non_morphed_mass))
      ("do-morph", po::value<bool>(&do_morph)->default_value(do_morph))
      ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(auto_rebin))
      ("manual_rebin", po::value<bool>(&manual_rebin)->default_value(manual_rebin))
      ("real_data", po::value<bool>(&real_data)->default_value(real_data))
      ("use_automc", po::value<bool>(&use_automc)->default_value(use_automc))
      ("verbose", po::value<bool>(&verbose)->default_value(verbose))
      ("output_folder", po::value<string>(&output_folder)->default_value(output_folder))
      ("analysis", po::value<string>(&analysis)->default_value(analysis))
      ("sub-analysis", po::value<string>(&sub_analysis)->default_value(sub_analysis))
      ("hSM-treatment", po::value<string>(&hSM_treatment)->default_value(hSM_treatment))
      ("sm-like-hists", po::value<string>(&sm_like_hists)->default_value(sm_like_hists))
      ("categorization", po::value<string>(&categorization)->default_value(categorization))
      ("era", po::value<int>(&era)->default_value(era))
      ("no-emb,no-emb,no_emb", po::bool_switch(&no_emb), "use MC samples instead of embedding")
      ("debug,d", po::bool_switch(&debug), "debug printout")
      ("mva", po::bool_switch(&mva), "mva tau id is used")
      ("sm", po::value<bool>(&sm)->default_value(sm))
      ("split_sm_signal_cat", po::value<bool>(&split_sm_signal_cat)->default_value(split_sm_signal_cat))
      ("mass-susy-ggH,mass_susy_ggH", po::value<vector<string>>(&mass_susy_ggH)->multitoken(), "mass_susy_ggH")
      ("mass-susy-qqH,mass_susy_qqH", po::value<vector<string>>(&mass_susy_qqH)->multitoken(), "mass_susy_qqH")
      ("bkgs", po::value<vector<string>>(&parser_bkgs)->multitoken(), "backgrounds")
      ("bkgs_em", po::value<vector<string>>(&parser_bkgs_em)->multitoken(), "backgrounds-em")
      ("sm_signals", po::value<vector<string>>(&parser_sm_signals)->multitoken(), "sm_signals")
      ("main_sm_signals", po::value<vector<string>>(&parser_main_sm_signals)->multitoken(), "main_sm_signals")
      ("no_shape_systs", po::value<bool>(&no_shape_systs)->default_value(no_shape_systs))
      ("enable_bsm_lowmass", po::value<bool>(&enable_bsm_lowmass)->default_value(enable_bsm_lowmass))
      ("prop_plot", po::value<bool>(&prop_plot)->default_value(false))
      ("help", "produce help message");
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);
  if (vm.count("help"))
  {
      cout << config << "\n";
      return 0;
  }

  lowmass = (categorization == "lowmass");

  // Sanity check for options with choices
  // analysis option
  if(std::find(analysis_choices.begin(), analysis_choices.end(), analysis) == analysis_choices.end()){
    std::cout << "ERROR: wrong choice of 'analysis' option. Please choose from:\n\t";
    for(auto choice : analysis_choices){
      std::cout << choice << " ";
    }
    std::cout << std::endl;
    exit(1);
  }
  // hSM_treatment option
  if(analysis == "bsm-model-indep" || analysis == "bsm-model-dep-full"){
    if(std::find(hSM_treatment_choices.begin(), hSM_treatment_choices.end(), hSM_treatment) == hSM_treatment_choices.end()){
      std::cout << "ERROR: wrong choice of 'hSM_treatment' option. In case of analysis 'bsm-model-indep' or 'bsm-model-dep-full', please choose from:\n\t";
      for(auto choice : hSM_treatment_choices){
        std::cout << choice << " ";
      }
      std::cout << std::endl;
      exit(1);
    }
  }
  else if(analysis == "bsm-model-dep-additional"){
    if(hSM_treatment != "hSM-in-bg"){
      std::cout << "ERROR: wrong choice of 'hSM_treatment' option. In case of analysis or 'bsm-model-dep-additional', please choose from:\n\thSM-in-bg" << std::endl;
      exit(1);
    }
  }

  // sub_analysis option
  if(analysis == "bsm-model-dep-full" || analysis == "bsm-model-dep-additional"){
    if(std::find(sub_analysis_choices.begin(), sub_analysis_choices.end(), sub_analysis) == sub_analysis_choices.end()){
      std::cout << "ERROR: wrong choice of 'sub_analysis' option. In case of model-dependent analysis, please choose from:\n\t";
      for(auto choice : sub_analysis_choices){
        std::cout << choice << " ";
      }
      std::cout << std::endl;
      exit(1);
    }
  }
  // sm_like_hists option
  if(analysis == "bsm-model-dep-full"){
    if(std::find(sm_like_hists_choices.begin(), sm_like_hists_choices.end(), sm_like_hists) == sm_like_hists_choices.end()){
      std::cout << "ERROR: wrong choice of 'sm_like_hists' option. Please choose from:\n\t";
      for(auto choice : sm_like_hists_choices){
        std::cout << choice << " ";
      }
      std::cout << std::endl;
      exit(1);
    }
  }
  // categorization option
  if(std::find(categorization_choices.begin(), categorization_choices.end(), categorization) == categorization_choices.end()){
    std::cout << "ERROR: wrong choice of 'categorization' option. Please choose from:\n\t";
    for(auto choice : categorization_choices){
      std::cout << choice << " ";
    }
    std::cout << std::endl;
    exit(1);
  }
  // era option
  if(std::find(era_choices.begin(), era_choices.end(), era) == era_choices.end()){
    std::cout << "ERROR: wrong choice of 'era' option. Please choose from:\n\t";
    for(auto choice : era_choices){
      std::cout << choice << " ";
    }
    std::cout << std::endl;
    exit(1);
  }

  // Define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
  std::string era_tag;
  if (era == 2016) era_tag = "2016";
  else if (era == 2017) era_tag = "2017";
  else if (era == 2018) era_tag = "2018";
  else std::runtime_error("Given era is not implemented.");

  output_folder = output_folder + "_" + analysis;
  std::map<string, string> input_dir;
  if (base_path.back() != '/' ) base_path += "/";
  if (!boost::filesystem::exists(output_folder)) boost::filesystem::create_directories(output_folder);
  // input_dir["mt"] = base_path + "/" +era_tag + "_hig-19-010" + "/mt/";
  // input_dir["et"] = base_path + "/" +era_tag + "_hig-19-010" + "/et/";
  // input_dir["tt"] = base_path + "/" +era_tag + "_hig-19-010" + "/tt/";
  // input_dir["em"] = base_path + "/" +era_tag + "_hig-19-010" + "/em/";
  input_dir["mt"] = base_path + "/" +era_tag + "/mt/";
  input_dir["et"] = base_path + "/" +era_tag + "/et/";
  input_dir["tt"] = base_path + "/" +era_tag + "/tt/";
  input_dir["em"] = base_path + "/" +era_tag + "/em/";

  // Define channels
  VString chns;
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
  VString bkgs, bkgs_em, bkgs_tt, bkgs_HWW, sm_signals, main_sm_signals, bkgs_em_noCR;

  VString mssm_ggH_signals, mssm_ggH_signals_additional, mssm_ggH_signals_smlike, mssm_ggH_signals_scalar, mssm_ggH_signals_pseudoscalar;
  VString mssm_bbH_signals, mssm_bbH_signals_additional, mssm_bbH_signals_smlike, mssm_bbH_signals_scalar, mssm_bbH_signals_pseudoscalar;
  VString mssm_qqH_signals;
  VString ggX_signals;

  VString mssm_ggH_lowmass_signals, mssm_ggH_lowmass_signals_additional, mssm_ggH_lowmass_signals_smlike, mssm_ggH_lowmass_signals_scalar, mssm_ggH_lowmass_signals_pseudoscalar;
  VString mssm_bbH_lowmass_signals, mssm_bbH_lowmass_signals_additional, mssm_bbH_lowmass_signals_smlike, mssm_bbH_lowmass_signals_scalar, mssm_bbH_lowmass_signals_pseudoscalar;

  VString mssm_signals, mssm_lowmass_signals, qqh_bsm_signals, wh_bsm_signals, zh_bsm_signals;

  std::string smlike = "h";
  if(sub_analysis == "sm-like-light"){
    smlike = "h";
  }
  else if(sub_analysis == "sm-like-heavy"){
    smlike = "H";
  }
  else if(sub_analysis == "cpv"){
    smlike = "H1";
  }

  if (sm == true){
    sm_signals = {"WH125", "ZH125", "bbH125"};
  }
  else {
    sm_signals = {"bbH125"};
  }
  main_sm_signals = {"ggH125", "qqH125"}; // qqH125 for mt,et,tt,em contains VBF+VH
  update_vector_by_byparser(sm_signals, parser_sm_signals, "sm_signals");
  update_vector_by_byparser(main_sm_signals, parser_main_sm_signals, "main_sm_signals");

  if(analysis == "bsm-model-indep")
  {
    mssm_ggH_signals = {"ggh_t", "ggh_b", "ggh_i"};
    mssm_bbH_signals = {"bbh"};
    mssm_ggH_lowmass_signals = {"ggh_t_lowmass", "ggh_b_lowmass", "ggh_i_lowmass"};
    mssm_bbH_lowmass_signals = {"bbh_lowmass"};
  }
  else if(analysis == "bsm-model-dep-full" || analysis == "bsm-model-dep-additional")
  {
    if(sub_analysis == "sm-like-light")
    {
      if(sm_like_hists == "bsm")
      {
        mssm_ggH_signals_smlike = {"ggh_t", "ggh_b", "ggh_i"};
      }
      else if(sm_like_hists == "sm125")
      {
        mssm_ggH_signals_smlike = {"ggh"};
      }
      mssm_bbH_signals_smlike = {"bbh"};
      mssm_ggH_signals_scalar = {"ggH_t", "ggH_b", "ggH_i"};
      mssm_bbH_signals_scalar = {"bbH"};
      mssm_ggH_signals_pseudoscalar = {"ggA_t", "ggA_b", "ggA_i"};
      mssm_bbH_signals_pseudoscalar = {"bbA"};
      qqh_bsm_signals = {"qqh"};
      if (sm == true){
        wh_bsm_signals = {"Wh"};
        zh_bsm_signals = {"Zh"};
      }
    }
    else if(sub_analysis == "sm-like-heavy")
    {
      if(sm_like_hists == "bsm")
      {
        mssm_ggH_signals_smlike = {"ggH_t", "ggH_b", "ggH_i"};
      }
      else if(sm_like_hists == "sm125")
      {
        mssm_ggH_signals_smlike = {"ggH"};
      }
      mssm_bbH_signals_smlike = {"bbH"};
      mssm_ggH_signals_scalar = {"ggh_t", "ggh_b", "ggh_i"};
      mssm_bbH_signals_scalar = {"bbh"};
      mssm_ggH_signals_pseudoscalar = {"ggA_t", "ggA_b", "ggA_i"};
      mssm_bbH_signals_pseudoscalar = {"bbA"};
      qqh_bsm_signals = {"qqH"};
      if (sm == true){
        wh_bsm_signals = {"WH"};
        zh_bsm_signals = {"ZH"};
      }
    }
    else if(sub_analysis == "cpv") // caution! 'scalar' and 'pseudoscalar' are used here for lists only! physics-wise ill-defined!
    {
      if(sm_like_hists == "bsm")
      {
        mssm_ggH_signals_smlike = {"ggH1_t", "ggH1_b", "ggH1_i"};
      }
      else if(sm_like_hists == "sm125")
      {
        mssm_ggH_signals_smlike = {"ggH1"};
      }
      mssm_bbH_signals_smlike = {"bbH1"};
      mssm_ggH_signals_scalar = {"ggH2_t", "ggH2_b", "ggH2_i"};
      mssm_bbH_signals_scalar = {"bbH2"};
      mssm_ggH_signals_pseudoscalar = {"ggH3_t", "ggH3_b", "ggH3_i"};
      mssm_bbH_signals_pseudoscalar = {"bbH3"};
      qqh_bsm_signals = {"qqH1"};
      if (sm == true){
        wh_bsm_signals = {"WH1"};
        zh_bsm_signals = {"ZH1"};
      }
    }
    for(auto proc : mssm_ggH_signals_scalar){
        mssm_ggH_lowmass_signals_scalar.push_back(proc + "_lowmass");
    }
    for(auto proc : mssm_ggH_signals_pseudoscalar){
        mssm_ggH_lowmass_signals_pseudoscalar.push_back(proc + "_lowmass");
    }
    for(auto proc : mssm_ggH_signals_smlike){
        mssm_ggH_lowmass_signals_smlike.push_back(proc + "_lowmass");
    }

    for(auto proc : mssm_bbH_signals_scalar){
        mssm_bbH_lowmass_signals_scalar.push_back(proc + "_lowmass");
    }
    for(auto proc : mssm_bbH_signals_pseudoscalar){
        mssm_bbH_lowmass_signals_pseudoscalar.push_back(proc + "_lowmass");
    }
    for(auto proc : mssm_bbH_signals_smlike){
        mssm_bbH_lowmass_signals_smlike.push_back(proc + "_lowmass");
    }

    mssm_ggH_signals_additional = ch::JoinStr({mssm_ggH_signals_scalar, mssm_ggH_signals_pseudoscalar});
    mssm_bbH_signals_additional = ch::JoinStr({mssm_bbH_signals_scalar, mssm_bbH_signals_pseudoscalar});
    mssm_ggH_lowmass_signals_additional = ch::JoinStr({mssm_ggH_lowmass_signals_scalar, mssm_ggH_lowmass_signals_pseudoscalar});
    mssm_bbH_lowmass_signals_additional = ch::JoinStr({mssm_bbH_lowmass_signals_scalar, mssm_bbH_lowmass_signals_pseudoscalar});

    mssm_ggH_signals = ch::JoinStr({mssm_ggH_signals_smlike, mssm_ggH_signals_scalar, mssm_ggH_signals_pseudoscalar});
    mssm_bbH_signals = ch::JoinStr({mssm_bbH_signals_smlike, mssm_bbH_signals_scalar, mssm_bbH_signals_pseudoscalar});
    mssm_ggH_lowmass_signals = ch::JoinStr({mssm_ggH_lowmass_signals_smlike, mssm_ggH_lowmass_signals_scalar, mssm_ggH_lowmass_signals_pseudoscalar});
    mssm_bbH_lowmass_signals = ch::JoinStr({mssm_bbH_lowmass_signals_smlike, mssm_bbH_lowmass_signals_scalar, mssm_bbH_lowmass_signals_pseudoscalar});
  }
  if (variable=="m_sv_puppi" || variable=="m_sv_VS_pt_tt_splitpT" || lowmass) {
    mssm_qqH_signals = {"qqX"};
    ggX_signals = {"ggX_t","ggX_b","ggX_i"};
  }
  mssm_signals = ch::JoinStr({mssm_ggH_signals, mssm_bbH_signals, mssm_qqH_signals, ggX_signals});
  mssm_lowmass_signals = ch::JoinStr({mssm_ggH_lowmass_signals, mssm_bbH_lowmass_signals});


  std::cout << "Used BSM signals: ";
  for(auto proc : mssm_signals){
    std::cout << proc << " ";
  }
  std::cout << std::endl;

  std::cout << "Used BSM lowmass signals: ";
  for(auto proc : mssm_lowmass_signals){
    std::cout << proc << " ";
  }
  std::cout << std::endl;

  bkgs = {"EMB", "ZL", "TTL", "VVL", "jetFakes"};
  bkgs_tt = {"EMB", "ZL", "TTL", "VVL", "jetFakes", "wFakes"};
  bkgs_HWW = {"ggHWW125", "qqHWW125", "WHWW125", "ZHWW125"};
  bkgs_em = {"EMB", "W", "ZL", "TTL", "VVL"};
  bkgs_em_noCR = {"QCD"};
  if ( sm == true){
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "jetFakes"), bkgs.end());
    bkgs.push_back("jetFakesSM");

    bkgs_tt.erase(std::remove(bkgs_tt.begin(), bkgs_tt.end(), "wFakes"), bkgs_tt.end());
    bkgs_tt.erase(std::remove(bkgs_tt.begin(), bkgs_tt.end(), "jetFakes"), bkgs_tt.end());
    bkgs_tt.push_back("jetFakesSM");
  }
  update_vector_by_byparser(bkgs, parser_bkgs, "bkgs");
  update_vector_by_byparser(bkgs_tt, parser_bkgs, "bkgs_tt");
  update_vector_by_byparser(bkgs_em, parser_bkgs_em, "bkgs_em");

  if (no_emb) {
    dout("WARNING: the EMB process is removed from backgrounds and ZTT, TTT and VVT templates are added");
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "EMB"), bkgs.end());
    bkgs_em.erase(std::remove(bkgs_em.begin(), bkgs_em.end(), "EMB"), bkgs_em.end());
    bkgs_tt.erase(std::remove(bkgs_tt.begin(), bkgs_tt.end(), "EMB"), bkgs_tt.end());
    bkgs.push_back("ZTT"); bkgs.push_back("TTT"); bkgs.push_back("VVT");
    bkgs_em.push_back("ZTT"); bkgs_em.push_back("TTT"); bkgs_em.push_back("VVT");
    bkgs_tt.push_back("ZTT"); bkgs_tt.push_back("TTT"); bkgs_tt.push_back("VVT");
  }
  map<int, VString> SUSYggH_masses;
  map<int, VString> SUSYbbH_masses;
  map<int, VString> SUSYqqH_masses;
  map<int, VString> ggX_masses;

  map<int, VString> SUSYggH_lowmasses;
  map<int, VString> SUSYbbH_lowmasses;

  if(do_morph) {

    if(lowmass) {
      SUSYqqH_masses[2018] = {"95"};
      SUSYqqH_masses[2017] = {"95"};
      SUSYqqH_masses[2016] = {"95"};

      ggX_masses[2018] = {"95"};
      ggX_masses[2017] = {"95"};
      ggX_masses[2016] = {"95"};
    }

    if(variable=="m_sv_VS_pt_tt" || variable=="m_sv_puppi" || variable=="m_sv_VS_pt_tt_splitpT") {

      SUSYbbH_masses[2018] = {"60","80","100","120","125","130","140","160","180","200","250"};
      SUSYbbH_masses[2017] = SUSYbbH_masses[2018];
      SUSYbbH_masses[2016] = SUSYbbH_masses[2018]; 

      SUSYggH_masses[2018] = {"60","80","95","100","120","125","130","140","160","180","200","250"};
      SUSYggH_masses[2016] = SUSYggH_masses[2018];
      SUSYggH_masses[2017] = SUSYggH_masses[2018];

      SUSYbbH_lowmasses[2018] = SUSYbbH_masses[2018]; 
      SUSYbbH_lowmasses[2017] = SUSYbbH_masses[2018];
      SUSYbbH_lowmasses[2016] = SUSYbbH_masses[2018]; 

      SUSYggH_lowmasses[2018] = SUSYggH_masses[2018];
      SUSYggH_lowmasses[2016] = SUSYggH_masses[2018];
      SUSYggH_lowmasses[2017] = SUSYggH_masses[2018];

      SUSYqqH_masses[2018] = {"95"};
      SUSYqqH_masses[2017] = {"95"};
      SUSYqqH_masses[2016] = {"95"};

      ggX_masses[2018] = {"95"};
      ggX_masses[2017] = {"95"};
      ggX_masses[2016] = {"95"};

    } else {

      SUSYbbH_masses[2018] = {"60","80","100","120","125","130","140","160","180","200","250","300","350","400","450","500","600","700","800","900","1000","1200","1400","1600","1800","2000","2300","2600","2900","3200","3500"};
      SUSYbbH_masses[2017] = SUSYbbH_masses[2018];
      SUSYbbH_masses[2016] = {"60","80","100","120","125","130","140","160","180","200","250","350","400","450","500","600","800","900","1200","1400","1600","1800","2000","2300","2600","2900","3200","3500"};  // Missing 300,700,1000
        
      SUSYggH_masses[2018] = {"60","80","100","120","125","130","140","160","180","200","250","300","350","400","450","500","600","700","800","900","1000","1200","1400","1600","1800","2000","2300","2600","2900","3200","3500"};
      SUSYggH_masses[2016] = SUSYggH_masses[2018];
      SUSYggH_masses[2017] = SUSYggH_masses[2018];

      SUSYbbH_lowmasses[2018] = {"60","80","100","120","125","130","140","160","180","200","250","300","350","400","450","500","600","700","800"};
      SUSYbbH_lowmasses[2017] = SUSYbbH_lowmasses[2018];
      SUSYbbH_lowmasses[2016] = {"60","80","100","120","125","130","140","160","180","200","250","350","400","450","500","600","800"};  // Missing 300,700

      SUSYggH_lowmasses[2018] = {"60","80","100","120","125","130","140","160","180","200","250","300","350","400","450","500","600","700","800"};
      SUSYggH_lowmasses[2016] = SUSYggH_lowmasses[2018];
      SUSYggH_lowmasses[2017] = SUSYggH_lowmasses[2018];
    }
  } else {
    // dont use mass morphing - need to specify a mass here
    SUSYggH_masses[2016] = {non_morphed_mass};
    SUSYggH_masses[2017] = {non_morphed_mass};
    SUSYggH_masses[2018] = {non_morphed_mass};
    if (non_morphed_mass!="95") {
      SUSYbbH_masses[2016] = {non_morphed_mass};
      SUSYbbH_masses[2017] = {non_morphed_mass};
      SUSYbbH_masses[2018] = {non_morphed_mass};
    }
  }

  if(prop_plot) {
    SUSYggH_masses[2016] = {"100"};
    SUSYggH_masses[2017] = {"100"};
    SUSYggH_masses[2018] = {"100"};
    SUSYbbH_masses[2016] = {"100"};
    SUSYbbH_masses[2017] = {"100"};
    SUSYbbH_masses[2018] = {"100"};

    SUSYqqH_masses[2018] = {};
    SUSYqqH_masses[2017] = {};
    SUSYqqH_masses[2016] = {};
    
    ggX_masses[2018] = {};
    ggX_masses[2017] = {};
    ggX_masses[2016] = {};
  }

  update_vector_by_byparser(SUSYggH_masses[era], mass_susy_ggH, "SUSY ggH");
  update_vector_by_byparser(SUSYbbH_masses[era], mass_susy_qqH, "SUSY qqH");

  std::cout << "[INFO] Considering the following processes as main backgrounds:\n";

  if (chan.find("em") != std::string::npos || chan.find("all") != std::string::npos) {
    std::cout << "For em channel : \n\t";
    printVector(bkgs_em);
  }
  if (chan.find("tt") != std::string::npos || chan.find("all") != std::string::npos) {
    std::cout << "For tt channels : \n\t";
    printVector(bkgs_tt);
  }
  if (chan.find("mt") != std::string::npos || chan.find("et") != std::string::npos || chan.find("all") != std::string::npos) {
    std::cout << "For et,mt channels : \n\t";
    printVector(bkgs);
  }
  bkg_procs["et"] = bkgs;
  bkg_procs["mt"] = bkgs;
  bkg_procs["tt"] = bkgs_tt;
  bkg_procs["em"] = bkgs_em;

  if(analysis == "sm"){
    bkg_procs["tt"] = JoinStr({bkg_procs["tt"],bkgs_HWW});
    bkg_procs["mt"] = JoinStr({bkg_procs["mt"],bkgs_HWW});
    bkg_procs["et"] = JoinStr({bkg_procs["et"],bkgs_HWW});
    bkg_procs["em"] = JoinStr({bkg_procs["em"],bkgs_HWW});
  }
  else if((analysis == "bsm-model-indep" && hSM_treatment == "hSM-in-bg") || analysis == "bsm-model-dep-additional"){
    bkg_procs["tt"] = JoinStr({bkg_procs["tt"],main_sm_signals,sm_signals});
    bkg_procs["mt"] = JoinStr({bkg_procs["mt"],main_sm_signals,sm_signals});
    bkg_procs["et"] = JoinStr({bkg_procs["et"],main_sm_signals,sm_signals});
    bkg_procs["em"] = JoinStr({bkg_procs["em"],main_sm_signals,sm_signals,bkgs_HWW});
    if(category == "et_xxh" || category == "et_tt" || category == "et_zll" || category == "et_misc" || category == "et_emb" || category == "et_ff" ||
       category == "et_xxh_bin_1" || category == "et_xxh_bin_2" || category == "et_xxh_bin_3" || category == "et_xxh_bin_4" || category == "et_xxh_bin_5" || category == "et_xxh_bin_6"){
      bkg_procs["et"] = JoinStr({bkg_procs["et"],sm_signals,main_sm_signals,bkgs_HWW});
    }
    else if(category == "mt_xxh" || category == "mt_tt" || category == "mt_zll" || category == "mt_misc" || category == "mt_emb" || category == "mt_ff" ||
            category == "mt_xxh_bin_1" || category == "mt_xxh_bin_2" || category == "mt_xxh_bin_3" || category == "mt_xxh_bin_4" || category == "mt_xxh_bin_5" || category == "mt_xxh_bin_6"){
      bkg_procs["mt"] = JoinStr({bkg_procs["mt"],sm_signals,main_sm_signals,bkgs_HWW});
    }
    else if(category == "tt_xxh" || category == "tt_misc" || category == "tt_emb" || category == "tt_ff" ||
            category == "tt_xxh_bin_1" || category == "tt_xxh_bin_2" || category == "tt_xxh_bin_3" || category == "tt_xxh_bin_4" || category == "tt_xxh_bin_5" || category == "tt_xxh_bin_6"){
      bkg_procs["tt"] = JoinStr({bkg_procs["tt"],sm_signals,main_sm_signals,bkgs_HWW});
    }
    else if(category == "em_xxh" || category == "em_tt" || category == "em_ss" || category == "em_misc" || category == "em_db" || category == "em_emb" ||
            category == "em_xxh_bin_1" || category == "em_xxh_bin_2" || category == "em_xxh_bin_3" || category == "em_xxh_bin_4" || category == "em_xxh_bin_5" || category == "em_xxh_bin_6"){
      bkg_procs["em"] = JoinStr({bkg_procs["em"], sm_signals, main_sm_signals,bkgs_HWW});
    }
  }
  else if(analysis == "bsm-model-dep-full"){
    bkg_procs["em"] = JoinStr({bkg_procs["em"],bkgs_HWW});
    if(category == "et_xxh" || category == "et_tt" || category == "et_zll" || category == "et_misc" || category == "et_emb" || category == "et_ff" ||
       category == "et_xxh_bin_1" || category == "et_xxh_bin_2" || category == "et_xxh_bin_3" || category == "et_xxh_bin_4" || category == "et_xxh_bin_5" || category == "et_xxh_bin_6"){
      bkg_procs["et"] = JoinStr({bkg_procs["et"],bkgs_HWW});
    }
    else if(category == "mt_xxh" || category == "mt_tt" || category == "mt_zll" || category == "mt_misc" || category == "mt_emb" || category == "mt_ff" ||
            category == "mt_xxh_bin_1" || category == "mt_xxh_bin_2" || category == "mt_xxh_bin_3" || category == "mt_xxh_bin_4" || category == "mt_xxh_bin_5" || category == "mt_xxh_bin_6"){
      bkg_procs["mt"] = JoinStr({bkg_procs["mt"],bkgs_HWW});
    }
    else if(category == "tt_xxh" || category == "tt_misc" || category == "tt_emb" || category == "tt_ff" ||
            category == "tt_xxh_bin_1" || category == "tt_xxh_bin_2" || category == "tt_xxh_bin_3" || category == "tt_xxh_bin_4" || category == "tt_xxh_bin_5" || category == "tt_xxh_bin_6"){
      bkg_procs["tt"] = JoinStr({bkg_procs["tt"],bkgs_HWW});
    }
    else if(category == "em_xxh" || category == "em_tt" || category == "em_ss" || category == "em_misc" || category == "em_db" || category == "em_emb" ||
            category == "em_xxh_bin_1" || category == "em_xxh_bin_2" || category == "em_xxh_bin_3" || category == "em_xxh_bin_4" || category == "em_xxh_bin_5" || category == "em_xxh_bin_6"){
      bkg_procs["em"] = JoinStr({bkg_procs["em"],bkgs_HWW});
    }
  }

  std::map< int, std::map<std::string,int> > SM_thresholds_bbH{
    {2016,{{"et", 1900}, {"mt", 1900}, {"tt", 1900}, {"em", 1900}}},
    {2017,{{"et", 1900}, {"mt", 1900}, {"tt", 1900}, {"em", 1900}}},
    {2018,{{"et", 1900}, {"mt", 1900}, {"tt", 1900}, {"em", 1900}}}};

  std::map< int, std::map<std::string,int> > SM_thresholds_ggH{
    {2016,{{"et", 1900}, {"mt", 1900}, {"tt", 1900}, {"em", 1900}}},
    {2017,{{"et", 1900}, {"mt", 1900}, {"tt", 1900}, {"em", 1900}}},
    {2018,{{"et", 1900}, {"mt", 1900}, {"tt", 1900}, {"em", 1900}}}};

  // Define MSSM model-dependent mass parameters mA, mH, mh
  RooRealVar mA("mA", "mA", 125., 90., 4000.);
  RooRealVar mH("mH", "mH", 125., 90., 4000.);
  RooRealVar mh("mh", "mh", 125., 90., 4000.);

  std::string max_lowmass = "60";
  if(SUSYggH_lowmasses.size()>0) SUSYggH_lowmasses[2018].back(); // this is set the same for all years for the time-being

  TString expression = max_lowmass + "*(mA >=" + max_lowmass +") + mA*(mA < "+ max_lowmass + ")";
  RooFormulaVar mA_lowmass("mA_lowmass", "mA_lowmass", expression, mA);
  expression = max_lowmass + "*(mH >=" + max_lowmass +") + mH*(mH < "+ max_lowmass + ")";
  RooFormulaVar mH_lowmass("mH_lowmass", "mH_lowmass", expression, mH);
  expression = max_lowmass + "*(mh >=" + max_lowmass +") + mh*(mh < "+ max_lowmass + ")";
  RooFormulaVar mh_lowmass("mh_lowmass", "mh_lowmass", expression, mh);

  // mA is used as model parameter in case of sub_analysis "sm-like-light", for "sm-like-heavy" and "cpv" it is mHp
  if(sub_analysis == "sm-like-light")
  {
    mA.setConstant(true);
  }

  // Define MSSM CPV model-dependent mass parameters mH3, mH2, mH1 (sub_analysis "cpv")
  RooRealVar mH3("mH3", "mH3", 125., 90., 4000.);
  RooRealVar mH2("mH2", "mH2", 125., 90., 4000.);
  RooRealVar mH1("mH1", "mH1", 125., 90., 4000.);

  expression = max_lowmass + "*(mH3 >=" + max_lowmass +") + mH3*(mH3 < "+ max_lowmass + ")";
  RooFormulaVar mH3_lowmass("mH3_lowmass", "mH3_lowmass", expression, mH3);
  expression = max_lowmass + "*(mH2 >=" + max_lowmass +") + mH2*(mH2 < "+ max_lowmass + ")";
  RooFormulaVar mH2_lowmass("mH2_lowmass", "mH2_lowmass", expression, mH2);
  expression = max_lowmass + "*(mH1 >=" + max_lowmass +") + mH1*(mH1 < "+ max_lowmass + ")";
  RooFormulaVar mH1_lowmass("mH1_lowmass", "mH1_lowmass", expression, mH1);

  // Define MSSM model-independent mass parameter MH
  RooRealVar MH("MH", "MH", 125., 60., 4000.);
  MH.setConstant(true);

  // Define categories
  Categories sm_signal = {};
  map<string, Categories> cats;
  if(categorization == "classic"){
    cats["et"] = {
        { 32, "et_Nbtag0_MTLt40"},
        { 33, "et_Nbtag0_MT40To70"},
        { 35, "et_NbtagGt1_MTLt40"},
        { 36, "et_NbtagGt1_MT40To70"},
    };
    cats["mt"] = {
        { 32, "mt_Nbtag0_MTLt40"},
        { 33, "mt_Nbtag0_MT40To70"},
        { 35, "mt_NbtagGt1_MTLt40"},
        { 36, "mt_NbtagGt1_MT40To70"},
    };
    cats["tt"] = {
        { 32, "tt_Nbtag0"},
        { 35, "tt_NbtagGt1"},
    };
    cats["em"] = {
        {  2, "em_NbtagGt1_DZetaLtm35"},
        { 32, "em_Nbtag0_DZetaGt30"},
        { 33, "em_Nbtag0_DZetam10To30"},
        { 34, "em_Nbtag0_DZetam35Tom10"},
        { 35, "em_NbtagGt1_DZetaGt30"},
        { 36, "em_NbtagGt1_DZetam10To30"},
        { 37, "em_NbtagGt1_DZetam35Tom10"},
    };
  } else if(categorization == "lowmass"){
    cats["et"] = {
        { 32, "et_Nbtag0_MTLt40"},
        { 33, "et_Nbtag0_MT40To70"},
        { 35, "et_NbtagGt1_MTLt40"},
        { 36, "et_NbtagGt1_MT40To70"},
        { 132, "et_Nbtag0_MTLt40_pT_0To50"},
        { 133, "et_Nbtag0_MT40To70_pT_0To50"},
        { 232, "et_Nbtag0_MTLt40_pT_50To100"},
        { 233, "et_Nbtag0_MT40To70_pT_50To100"},
        { 332, "et_Nbtag0_MTLt40_pT_100To200"},
        { 333, "et_Nbtag0_MT40To70_pT_100To200"},
        { 432, "et_Nbtag0_MTLt40_pT_GT200"},
        { 433, "et_Nbtag0_MT40To70_pT_GT200"},
    };
    cats["mt"] = {
        { 32, "mt_Nbtag0_MTLt40"},
        { 33, "mt_Nbtag0_MT40To70"},
        { 35, "mt_NbtagGt1_MTLt40"},
        { 36, "mt_NbtagGt1_MT40To70"},
        { 132, "mt_Nbtag0_MTLt40_pT_0To50"},
        { 133, "mt_Nbtag0_MT40To70_pT_0To50"},
        { 232, "mt_Nbtag0_MTLt40_pT_50To100"},
        { 233, "mt_Nbtag0_MT40To70_pT_50To100"},
        { 332, "mt_Nbtag0_MTLt40_pT_100To200"},
        { 333, "mt_Nbtag0_MT40To70_pT_100To200"},
        { 432, "mt_Nbtag0_MTLt40_pT_GT200"},
        { 433, "mt_Nbtag0_MT40To70_pT_GT200"},
    };
    cats["tt"] = {
        { 32, "tt_Nbtag0"},
        { 35, "tt_NbtagGt1"},
        { 132, "tt_Nbtag0_pT_0To50"},
        { 232, "tt_Nbtag0_pT_50To100"},
        { 332, "tt_Nbtag0_pT_100To200"},
        { 432, "tt_Nbtag0_pT_GT200"},
    };
    cats["em"] = {
        {  2, "em_NbtagGt1_DZetaLtm35"},
        { 32, "em_Nbtag0_DZetaGt30"},
        { 33, "em_Nbtag0_DZetam10To30"},
        { 34, "em_Nbtag0_DZetam35Tom10"},
        { 35, "em_NbtagGt1_DZetaGt30"},
        { 36, "em_NbtagGt1_DZetam10To30"},
        { 37, "em_NbtagGt1_DZetam35Tom10"},
        { 132, "em_Nbtag0_DZetaGt30_pT_0To50"},
        { 133, "em_Nbtag0_DZetam10To30_pT_0To50"},
        { 134, "em_Nbtag0_DZetam35Tom10_pT_0To50"},
        { 232, "em_Nbtag0_DZetaGt30_pT_50To100"},
        { 233, "em_Nbtag0_DZetam10To30_pT_50To100"},
        { 234, "em_Nbtag0_DZetam35Tom10_pT_50To100"},
        { 332, "em_Nbtag0_DZetaGt30_pT_100To200"},
        { 333, "em_Nbtag0_DZetam10To30_pT_100To200"},
        { 334, "em_Nbtag0_DZetam35Tom10_pT_100To200"},
        { 432, "em_Nbtag0_DZetaGt30_pT_GT200"},
        { 433, "em_Nbtag0_DZetam10To30_pT_GT200"},
        { 434, "em_Nbtag0_DZetam35Tom10_pT_GT200"},
    };
  }
  else if(categorization == "sm-ml-only"){
    cats["et"] = {
      {13, "et_tt"},
      {15, "et_zll"},
      {16, "et_misc"},
      {20, "et_emb"},
      {21, "et_ff"}
    };
    if(split_sm_signal_cat){
      sm_signal = {
        // Split SM Signal Categories
        { 101, "et_xxh_bin_1"},
        { 102, "et_xxh_bin_2"},
        { 103, "et_xxh_bin_3"},
        { 104, "et_xxh_bin_4"},
        { 105, "et_xxh_bin_5"},
        { 106, "et_xxh_bin_6"}
      };
    }
    else{
      sm_signal = {
        // 2D SM Signal Category
        { 101, "et_xxh"}
      };
    }
    cats["et"].reserve(cats["et"].size() + std::distance(sm_signal.begin(), sm_signal.end()));
    cats["et"].insert(cats["et"].end(), sm_signal.begin(), sm_signal.end());

    cats["mt"] = {
      {13, "mt_tt"},
      {15, "mt_zll"},
      {16, "mt_misc"},
      {20, "mt_emb"},
      {21, "mt_ff"}
    };
    if(split_sm_signal_cat){
      sm_signal = {
        // Split SM Signal Categories
        { 101, "mt_xxh_bin_1"},
        { 102, "mt_xxh_bin_2"},
        { 103, "mt_xxh_bin_3"},
        { 104, "mt_xxh_bin_4"},
        { 105, "mt_xxh_bin_5"},
        { 106, "mt_xxh_bin_6"}
      };
    }
    else{
      sm_signal = {
        // 2D SM Signal Category
        { 101, "mt_xxh"}
      };
    }
    cats["mt"].reserve(cats["mt"].size() + std::distance(sm_signal.begin(), sm_signal.end()));
    cats["mt"].insert(cats["mt"].end(), sm_signal.begin(), sm_signal.end());

    cats["tt"] = {
      {16, "tt_misc"},
      {20, "tt_emb"},
      {21, "tt_ff"}
    };
    if(split_sm_signal_cat){
      sm_signal = {
        // Split SM Signal Categories
        { 101, "tt_xxh_bin_1"},
        { 102, "tt_xxh_bin_2"},
        { 103, "tt_xxh_bin_3"},
        { 104, "tt_xxh_bin_4"},
        { 105, "tt_xxh_bin_5"},
        { 106, "tt_xxh_bin_6"}
      };
    }
    else{
      sm_signal = {
        // 2D SM Signal Category
        { 101, "tt_xxh"}
      };
    }
    cats["tt"].reserve(cats["tt"].size() + std::distance(sm_signal.begin(), sm_signal.end()));
    cats["tt"].insert(cats["tt"].end(), sm_signal.begin(), sm_signal.end());

    cats["em"] = {
      {13, "em_tt"},
      {14, "em_ss"},
      {16, "em_misc"},
      {19, "em_db"},
      {20, "em_emb"}
    };
    if(split_sm_signal_cat){
      sm_signal = {
        // Split SM Signal Categories
        { 101, "em_xxh_bin_1"},
        { 102, "em_xxh_bin_2"},
        { 103, "em_xxh_bin_3"},
        { 104, "em_xxh_bin_4"},
        { 105, "em_xxh_bin_5"},
        { 106, "em_xxh_bin_6"}
      };
    }
    else{
      sm_signal = {
        // 2D SM Signal Category
        { 101, "em_xxh"}
      };
    }
    cats["em"].reserve(cats["em"].size() + std::distance(sm_signal.begin(), sm_signal.end()));
    cats["em"].insert(cats["em"].end(), sm_signal.begin(), sm_signal.end());
  }
  else if(categorization == "with-sm-ml"){
    cats["et"] = {
        {13, "et_tt"},
        {15, "et_zll"},
        {16, "et_misc"},
        {20, "et_emb"},
        {21, "et_ff"},

        {32, "et_Nbtag0_MTLt40_MHGt250"},
        {33, "et_Nbtag0_MT40To70_MHGt250"},
        {35, "et_NbtagGt1_MTLt40"},
        {36, "et_NbtagGt1_MT40To70"},
    };
    if(split_sm_signal_cat){
      sm_signal = {
        // Split SM Signal Categories
        { 101, "et_xxh_bin_1"},
        { 102, "et_xxh_bin_2"},
        { 103, "et_xxh_bin_3"},
        { 104, "et_xxh_bin_4"},
        { 105, "et_xxh_bin_5"},
        { 106, "et_xxh_bin_6"}
      };
    }
    else{
      sm_signal = {
        // 2D SM Signal Category
        { 101, "et_xxh"}
      };
    }
    cats["et"].reserve(cats["et"].size() + std::distance(sm_signal.begin(), sm_signal.end()));
    cats["et"].insert(cats["et"].end(), sm_signal.begin(), sm_signal.end());

    cats["mt"] = {
        {13, "mt_tt"},
        {15, "mt_zll"},
        {16, "mt_misc"},
        {20, "mt_emb"},
        {21, "mt_ff"},

        {32, "mt_Nbtag0_MTLt40_MHGt250"},
        {33, "mt_Nbtag0_MT40To70_MHGt250"},

        {35, "mt_NbtagGt1_MTLt40"},
        {36, "mt_NbtagGt1_MT40To70"},
    };
    if(split_sm_signal_cat){
      sm_signal = {
        // Split SM Signal Categories
        { 101, "mt_xxh_bin_1"},
        { 102, "mt_xxh_bin_2"},
        { 103, "mt_xxh_bin_3"},
        { 104, "mt_xxh_bin_4"},
        { 105, "mt_xxh_bin_5"},
        { 106, "mt_xxh_bin_6"}
      };
    }
    else{
      sm_signal = {
        // 2D SM Signal Category
        { 101, "mt_xxh"}
      };
    }
    cats["mt"].reserve(cats["mt"].size() + std::distance(sm_signal.begin(), sm_signal.end()));
    cats["mt"].insert(cats["mt"].end(), sm_signal.begin(), sm_signal.end());

    cats["tt"] = {
        {16, "tt_misc"},
        {20, "tt_emb"},
        {21, "tt_ff"},

        {32, "tt_Nbtag0_MHGt250"},

        {35, "tt_NbtagGt1"},
    };
    if(split_sm_signal_cat){
      sm_signal = {
        // Split SM Signal Categories
        { 101, "tt_xxh_bin_1"},
        { 102, "tt_xxh_bin_2"},
        { 103, "tt_xxh_bin_3"},
        { 104, "tt_xxh_bin_4"},
        { 105, "tt_xxh_bin_5"},
        { 106, "tt_xxh_bin_6"}
      };
    }
    else{
      sm_signal = {
        // 2D SM Signal Category
        { 101, "tt_xxh"}
      };
    }
    cats["tt"].reserve(cats["tt"].size() + std::distance(sm_signal.begin(), sm_signal.end()));
    cats["tt"].insert(cats["tt"].end(), sm_signal.begin(), sm_signal.end());

    cats["em"] = {
        { 2, "em_NbtagGt1_DZetaLtm35"},

        {13, "em_tt"},
        {14, "em_ss"},
        {16, "em_misc"},
        {19, "em_db"},
        {20, "em_emb"},

        {32, "em_Nbtag0_DZetaGt30_MHGt250"},
        {33, "em_Nbtag0_DZetam10To30_MHGt250"},
        {34, "em_Nbtag0_DZetam35Tom10_MHGt250"},

        {35, "em_NbtagGt1_DZetaGt30"},
        {36, "em_NbtagGt1_DZetam10To30"},
        {37, "em_NbtagGt1_DZetam35Tom10"},
    };
    if(split_sm_signal_cat){
      sm_signal = {
        // Split SM Signal Categories
        { 101, "em_xxh_bin_1"},
        { 102, "em_xxh_bin_2"},
        { 103, "em_xxh_bin_3"},
        { 104, "em_xxh_bin_4"},
        { 105, "em_xxh_bin_5"},
        { 106, "em_xxh_bin_6"}
      };
    }
    else{
      sm_signal = {
        // 2D SM Signal Category
        { 101, "em_xxh"}
      };
    }
    cats["em"].reserve(cats["em"].size() + std::distance(sm_signal.begin(), sm_signal.end()));
    cats["em"].insert(cats["em"].end(), sm_signal.begin(), sm_signal.end());

  }
  else throw std::runtime_error("Given categorization is not known.");

  // Create combine harverster object
  ch::CombineHarvester cb;
  cb.SetFlag("workspaces-use-clone", true);


  // Introduce ordering of categories for the final discriminator in MSSM
  std::vector<int> sm_categories = {13,14,15,16,19,20,21}; // Control regions from the ML SM HTT analysis
  std::vector<int> em_control_category = {2}; // Control region for em channel
  std::vector<int> mssm_btag_categories = {35,135,235,335,435,36,136,236,336,436,37,137,237,337,437}; // b-tagged MSSM-like categories with mt_tot as discriminator
  std::vector<int> mssm_nobtag_categories = {32,33,34,132,232,332,432,133,233,333,433,134,234,334,434}; // non-btagged MSSM-like categories with mt_tot as discriminator
  std::vector<int> sm_signal_category = {101,102,103,104,105,106}; // category for the SM signal


  for (auto chn : chns) {
  // build category maps used for the different analyses
    Categories sm_and_btag_cats = cats[chn]; // contain 101-106, 2, 13-21, 35-37
    Categories sm_cats = cats[chn]; // contain 101-106, 13-21
    Categories mssm_btag_cats = cats[chn]; // contain 2, 35-37
    Categories mssm_cats = cats[chn]; // contain 2, 32-37
    Categories exclude_em_control = cats[chn]; // contain all except 2
    Categories sm_signal_cat = cats[chn]; // contain 101-106

    auto catit = sm_signal_cat.begin();
    while(catit != sm_signal_cat.end())
    {
      if(std::find(sm_categories.begin(), sm_categories.end(), (*catit).first) != sm_categories.end()){
        sm_signal_cat.erase(catit);
      }
      else if(std::find(mssm_nobtag_categories.begin(), mssm_nobtag_categories.end(), (*catit).first) != mssm_nobtag_categories.end()){
        sm_signal_cat.erase(catit);
      }
      else if(std::find(mssm_btag_categories.begin(), mssm_btag_categories.end(), (*catit).first) != mssm_btag_categories.end()){
        sm_signal_cat.erase(catit);
      }
      else if(std::find(em_control_category.begin(), em_control_category.end(), (*catit).first) != em_control_category.end()){
        sm_signal_cat.erase(catit);
      }
      else
      {
        ++catit;
      }
    }

    catit = sm_cats.begin();
    while(catit != sm_cats.end())
    {
      if(std::find(mssm_nobtag_categories.begin(), mssm_nobtag_categories.end(), (*catit).first) != mssm_nobtag_categories.end()){
        sm_cats.erase(catit);
      }
      else if(std::find(mssm_btag_categories.begin(), mssm_btag_categories.end(), (*catit).first) != mssm_btag_categories.end()){
        sm_cats.erase(catit);
      }
      else if(std::find(em_control_category.begin(), em_control_category.end(), (*catit).first) != em_control_category.end()){
        sm_cats.erase(catit);
      }
      else
      {
        ++catit;
      }
    }

    catit = mssm_cats.begin();
    while(catit != mssm_cats.end())
    {
      if(std::find(sm_categories.begin(), sm_categories.end(), (*catit).first) != sm_categories.end()){
        mssm_cats.erase(catit);
      }
      else if(std::find(sm_signal_category.begin(), sm_signal_category.end(), (*catit).first) != sm_signal_category.end()){
        mssm_cats.erase(catit);
      }
      else
      {
        ++catit;
      }
    }

    Categories mssm_cats_exclude_em_control = mssm_cats;

    catit = mssm_cats_exclude_em_control.begin();
    while(catit != mssm_cats_exclude_em_control.end())
    {
      if(std::find(em_control_category.begin(), em_control_category.end(), (*catit).first) != em_control_category.end()){
        mssm_cats_exclude_em_control.erase(catit);
      }
      else
      {
        ++catit;
      }
    }

    catit = exclude_em_control.begin();
    while(catit != exclude_em_control.end())
    {
      if(std::find(em_control_category.begin(), em_control_category.end(), (*catit).first) != em_control_category.end()){
        exclude_em_control.erase(catit);
      }
      else
      {
        ++catit;
      }
    }

    catit = mssm_btag_cats.begin();
    while(catit != mssm_btag_cats.end())
    {
      if(std::find(mssm_nobtag_categories.begin(), mssm_nobtag_categories.end(), (*catit).first) != mssm_nobtag_categories.end()){
        mssm_btag_cats.erase(catit);
      }
      else if(std::find(sm_categories.begin(), sm_categories.end(), (*catit).first) != sm_categories.end()){
        mssm_btag_cats.erase(catit);
      }
      else if(std::find(sm_signal_category.begin(), sm_signal_category.end(), (*catit).first) != sm_signal_category.end()){
        mssm_btag_cats.erase(catit);
      }
      else
      {
        ++catit;
      }
    }

    catit = sm_and_btag_cats.begin();
    while(catit != sm_and_btag_cats.end())
    {
      if(std::find(mssm_nobtag_categories.begin(), mssm_nobtag_categories.end(), (*catit).first) != mssm_nobtag_categories.end()){
        sm_and_btag_cats.erase(catit);
      }
      else
      {
        ++catit;
      }
    }

    Categories sm_and_btag_cats_exclude_em_control = sm_and_btag_cats;

    catit = sm_and_btag_cats_exclude_em_control.begin();
    while(catit != sm_and_btag_cats_exclude_em_control.end())
    {
      if(std::find(em_control_category.begin(), em_control_category.end(), (*catit).first) != em_control_category.end()){
        sm_and_btag_cats_exclude_em_control.erase(catit);
      }
      else
      {
        ++catit;
      }
    }

    std::cout << "[INFO] Using the following categories:" << std::endl;
    std::cout << "   sm_and_btag_cats:" << std::endl;
    for (const auto i: sm_and_btag_cats)
      std::cout << "      " << i.first << ' ' << i.second << std::endl;
    std::cout  << std::endl;
    std::cout << "   sm_and_btag_cats_exclude_em_control:" << std::endl;
    for (const auto i: sm_and_btag_cats_exclude_em_control)
      std::cout << "      " << i.first << ' ' << i.second << std::endl;
    std::cout  << std::endl;
    std::cout << "    sm_cats:" << std::endl;
    for (const auto i: sm_cats)
      std::cout << "      " << i.first << ' ' << i.second << std::endl;
    std::cout  << std::endl;
    std::cout << "    mssm_cats:" << std::endl;
    for (const auto i: mssm_cats)
      std::cout << "      " << i.first << ' ' << i.second << std::endl;
    std::cout  << std::endl;
    std::cout << "    mssm_cats_exclude_em_control:" << std::endl;
    for (const auto i: mssm_cats_exclude_em_control)
      std::cout << "      " << i.first << ' ' << i.second << std::endl;
    std::cout  << std::endl;
    std::cout << "    mssm_btag_cats:" << std::endl;
    for (const auto i: mssm_btag_cats)
      std::cout << "      " << i.first << ' ' << i.second << std::endl;
    std::cout  << std::endl;
    std::cout << "    sm_signal_cat:" << std::endl;
    for (const auto i: sm_signal_cat)
      std::cout << "      " << i.first << ' ' << i.second << std::endl;
    std::cout  << std::endl;
    std::cout << "    exclude_em_control:" << std::endl;
    for (const auto i: exclude_em_control)
      std::cout << "      " << i.first << ' ' << i.second << std::endl;
    std::cout  << std::endl;

    cb.AddObservations({"*"}, {"htt"}, {era_tag}, {chn}, cats[chn]);
    // Adding background processes. This also contains SMH125 templates if configured in that way:
    // For bsm-model-indep analysis with Higgs boson in BG: ggH125, qqH125, bbH125, and WH125, ZH125
    // For bsm-model-dep-additional analysis: ggH125, qqH125, bbH125, and WH125, ZH125
    cb.AddProcesses({"*"}, {"htt"}, {era_tag}, {chn}, bkg_procs[chn], cats[chn], false);
    // Include QCD process in em channel for all categories except for CR
    if (chn == "em") {
        cb.AddProcesses({"*"}, {"htt"}, {era_tag}, {chn}, bkgs_em_noCR, exclude_em_control, false);
        //cb.AddProcesses({"*"}, {"htt"}, {era_tag}, {chn}, bkgs_em_noCR, cats[chn], false); // adding back QCD for tests
    }

    if(analysis == "sm"){
      cb.AddProcesses({""}, {"htt"}, {era_tag}, {chn}, ch::JoinStr({main_sm_signals, sm_signals}), cats[chn], true); // These are ggH125, qqH125, bbH125, WH125, ZH125
    }
    else if(analysis == "bsm-model-indep"){

      // Adding configured SUSY signals in all categories but the em control region 2 for bsm model-independent analyses
      cb.AddProcesses(SUSYbbH_masses[era], {"htt"}, {era_tag}, {chn}, mssm_bbH_signals, exclude_em_control, true);
      cb.AddProcesses(SUSYggH_masses[era], {"htt"}, {era_tag}, {chn}, mssm_ggH_signals, exclude_em_control, true);

      if(variable=="m_sv_puppi" || variable=="m_sv_VS_pt_tt_splitpT" || lowmass) {
        cb.AddProcesses(SUSYqqH_masses[era], {"htt"}, {era_tag}, {chn}, mssm_qqH_signals, exclude_em_control, true);
        cb.AddProcesses(ggX_masses[era], {"htt"}, {era_tag}, {chn}, ggX_signals, exclude_em_control, true);
      }
    }
    else if(analysis == "bsm-model-dep-additional" || analysis == "bsm-model-dep-full"){
      // Adding at first the additional Higgs boson signals
      Categories additional_higgses_cats;
      if(categorization == "classic" || categorization == "lowmass") // classic categories cover full mass range for additional Higgs signals
      {
         additional_higgses_cats = exclude_em_control;
      }
      else if(categorization == "with-sm-ml" || categorization == "sm-ml-only")
      {
         if(sub_analysis == "sm-like-light" || sub_analysis == "cpv") // in that case, the additional Higgs bosons are relatively heavy
         {
            if(enable_bsm_lowmass){
              additional_higgses_cats = mssm_cats_exclude_em_control;
            }
            else{
              additional_higgses_cats = exclude_em_control;
            }
         }
         else if(sub_analysis == "sm-like-heavy") // in that case, all additional Higgs bosons are relatively light (< 200 GeV) --> don't consider them in high mass no-btag categories
         {
             additional_higgses_cats = sm_and_btag_cats_exclude_em_control;
         }
      }
      // sm-like-light sub_analysis: H and A
      // sm-like-heavy sub_analysis: h and A
      // cpv sub_analysis: H2 and H3
      cb.AddProcesses(SUSYggH_masses[era], {"htt"}, {era_tag}, {chn}, mssm_ggH_signals_additional, additional_higgses_cats, true);
      cb.AddProcesses(SUSYbbH_masses[era], {"htt"}, {era_tag}, {chn}, mssm_bbH_signals_additional, additional_higgses_cats, true);

      if(enable_bsm_lowmass)
      {
        cb.AddProcesses(SUSYggH_lowmasses[era], {"htt"}, {era_tag}, {chn}, mssm_ggH_lowmass_signals_additional, sm_cats, true);
        cb.AddProcesses(SUSYbbH_lowmasses[era], {"htt"}, {era_tag}, {chn}, mssm_bbH_lowmass_signals_additional, sm_cats, true);
      }

      if(analysis == "bsm-model-dep-full") // In that case of analysis we compare full neutral BSM Higgs spectrum (h, H, A or H1, H2, H3) with SM hypothesis
      {
        // Defining categories for SM and BSM SM-like HTT signal: exclude mssm_nobtag_categories in case SM ML HTT categories are used because of m_sv >= 250 GeV cut
        Categories qq_gg_bb_phi_cats;
        if(categorization == "classic" || categorization == "lowmass")
        {
          qq_gg_bb_phi_cats = cats[chn];
        }
        else if(categorization == "with-sm-ml" || categorization == "sm-ml-only")
        {
          qq_gg_bb_phi_cats = sm_and_btag_cats;
        }

        // Adding SM Higgs processes as signal or background depending on hSM treatment for model-dependent analyses with full neutral Higgs modelling
        // (since testing then against SM Higgs + BG hypothesis)
        // These comprise: ggH125, qqH125, bbH125, and in case of sm categories WH125, ZH125
        cb.AddProcesses({""}, {"htt"}, {era_tag}, {chn}, ch::JoinStr({main_sm_signals, sm_signals}), qq_gg_bb_phi_cats, hSM_treatment == "no-hSM-in-bg");

        // Adding the qqphi process for all bsm model-dependent analyses with full neutral Higgs modelling
        // sm-like-light: phi = h
        // sm-like-heavy: phi = H
        // cpv: phi = H1
        cb.AddProcesses({""}, {"htt"}, {era_tag}, {chn}, qqh_bsm_signals, qq_gg_bb_phi_cats, true);
        // Wphi and Zphi are only added if an sm category is considered, since otherwise no templates are available for WH and ZH
        // sm-like-light: phi = h
        // sm-like-heavy: phi = H
        // cpv: phi = H1
        if(sm){
          cb.AddProcesses({""}, {"htt"}, {era_tag}, {chn}, wh_bsm_signals, qq_gg_bb_phi_cats, true);
          cb.AddProcesses({""}, {"htt"}, {era_tag}, {chn}, zh_bsm_signals, qq_gg_bb_phi_cats, true);
        }

        // Adding the SM-like processes bbphi and ggphi
        // sm-like-light: phi = h
        // sm-like-heavy: phi = H
        // cpv: phi = H1
        VString empty_masses = {""};
        VString ggH_SMlike_masses = (sm_like_hists == "sm125") ? empty_masses : SUSYggH_masses[era];
        cb.AddProcesses(ggH_SMlike_masses, {"htt"}, {era_tag}, {chn}, mssm_ggH_signals_smlike, qq_gg_bb_phi_cats, true);
        VString bbH_SMlike_masses = (sm_like_hists == "sm125") ? empty_masses : SUSYbbH_masses[era];
        cb.AddProcesses(bbH_SMlike_masses, {"htt"}, {era_tag}, {chn}, mssm_bbH_signals_smlike, qq_gg_bb_phi_cats, true);
      }
    }
  }
  dout("[INFO] Add systematics AddMSSMvsSMRun2Systematics, embedding:", ! no_emb, " sm categories:", sm);
  ch::AddMSSMvsSMRun2Systematics(cb, true, ! no_emb, true, true, true, era, mva, sm, smlike);

  if(variable=="m_sv_VS_pt_tt_splitpT") {
    //extra lnN uncertainties for FF when splitting by pt_tt bins
    cb.cp()
    .channel({"et", "mt", "tt"})
    .process({"jetFakes"})
    .AddSyst(cb, "CMS_ff_total_syst_pt_tt_bin_$BIN_$ERA", "lnN", SystMap<>::init(1.05));
  }
  if(variable=="m_sv_puppi") {
    // for btag categories we also use additional uncertainties to account for non-closures in m_sv distributions
    cb.cp()
    .channel({"et", "mt"})
    .bin_id({35})
    .process({"jetFakes"})
    .AddSyst(cb, "CMS_ff_total_ttbar_msv_shape_syst_mTtight_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
    .channel({"et", "mt"})
    .bin_id({36})
    .process({"jetFakes"})
    .AddSyst(cb, "CMS_ff_total_ttbar_msv_shape_syst_mTloose_$ERA", "shape", SystMap<>::init(1.00));

  }

  if(prop_plot){
    // shapeU seems to have issues for prop plots so change CMS_htt_ttbarShape to shape
    auto cb_syst = cb.cp().syst_name({"CMS_htt_ttbarShape"});
    cb_syst.ForEachSyst([&](ch::Systematic *syst) {
      syst->set_type("shape");
    });
  }

  dout("[INFO] Systematics added");
  // Define restriction to the desired category
  if(category != "all"){
    cb = cb.bin({category});
  }

  if(no_shape_systs){
    cb.FilterSysts([&](ch::Systematic *s){
      return s->type().find("shape") != std::string::npos;
    });
  }

  for (string chn : chns) {
    string input_file_base = input_dir[chn] + "htt_all.inputs-mssm-vs-sm-Run" + era_tag + "-" + variable + ".root";
    if (mva) input_file_base = input_dir[chn] + "htt_" + chn + ".inputs-mssm-vs-sm-" + era_tag + "-" + variable + ".root";
    dout("[INFO] Extracting shapes from ", input_file_base);

    // Adding background templates to processes. This also involves (if configured) SMH125 processes.
    // bsm-model-indep analysis with Higgs boson in BG: ggH125, qqH125, bbH125, and WH125, ZH125
    // bsm-model-dep-additional analysis: ggH125, qqH125, bbH125, and WH125, ZH125
    cb.cp().channel({chn}).backgrounds().process({"bbH125"}, false).ExtractShapes(
      input_file_base, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");
    cb.cp().channel({chn}).backgrounds().process({"bbH125"}).ExtractShapes( // "bbH125" needs special treatment because of template name spelling
      input_file_base, "$BIN/bbH_125", "$BIN/bbH_125_$SYSTEMATIC");

    if(analysis == "sm"){
      cb.cp().channel({chn}).process(ch::JoinStr({sm_signals,main_sm_signals})).process({"bbH125"}, false).ExtractShapes( // These are ggH125, qqH125, WH125, ZH125
        input_file_base, "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC");
      cb.cp().channel({chn}).process({"bbH125"}).ExtractShapes( // "bbH125" needs special treatment because of template name spelling
        input_file_base, "$BIN/bbH_125", "$BIN/bbH_125_$SYSTEMATIC");
    }
    // Adding templates for configured SUSY signals
    // Comprising BSM signal h in model-independent case
    else if(analysis == "bsm-model-indep"){
      cb.cp().channel({chn}).process(mssm_ggH_signals).ExtractShapes(
        input_file_base, "$BIN/$PROCESS_$MASS", "$BIN/$PROCESS_$MASS_$SYSTEMATIC");
      cb.cp().channel({chn}).process(mssm_bbH_signals).ExtractShapes(
        input_file_base, "$BIN/bbH_$MASS", "$BIN/bbH_$MASS_$SYSTEMATIC");
    }
    // Adding templates for configured SUSY signals of model-dependent analyses
    // sm-like-light sub_analysis: H and A
    // sm-like-heavy sub_analysis: h and A
    // cpv sub_analysis: H2 and H3
    else if(analysis == "bsm-model-dep-additional" || analysis == "bsm-model-dep-full"){
      // In case of ggPhi, it depends on the sub_analysis, how to include the templates
      // It is simpler, when the templates correspond to H, h, or A as it is in for
      // sub_analysis sm-like-light and sm-like-heavy
      if(sub_analysis == "sm-like-light" || sub_analysis == "sm-like-heavy"){
        cb.cp().channel({chn}).process(mssm_ggH_signals_additional).ExtractShapes(
          input_file_base, "$BIN/$PROCESS_$MASS", "$BIN/$PROCESS_$MASS_$SYSTEMATIC");
        if(enable_bsm_lowmass){
          for(auto ggH : mssm_ggH_lowmass_signals_additional){
            std::string template_ggH = boost::replace_all_copy(ggH, "_lowmass", "");
            cb.cp().channel({chn}).process({ggH}).ExtractShapes(
              input_file_base, "$BIN/" + template_ggH + "_$MASS", "$BIN/" + template_ggH + "_$MASS_$SYSTEMATIC");
          }
        }
      }
      // In case of sub_analysis cpv, the templates have to be included explicitly
      // for ggPhi due to different naming of the process
      else if(sub_analysis == "cpv"){
        cb.cp().channel({chn}).process({"ggH2_t"}).ExtractShapes(
          input_file_base, "$BIN/ggH_t_$MASS", "$BIN/ggH_t_$MASS_$SYSTEMATIC");
        cb.cp().channel({chn}).process({"ggH2_b"}).ExtractShapes(
          input_file_base, "$BIN/ggH_b_$MASS", "$BIN/ggH_b_$MASS_$SYSTEMATIC");
        cb.cp().channel({chn}).process({"ggH2_i"}).ExtractShapes(
          input_file_base, "$BIN/ggH_i_$MASS", "$BIN/ggH_i_$MASS_$SYSTEMATIC");

        cb.cp().channel({chn}).process({"ggH3_t"}).ExtractShapes(
          input_file_base, "$BIN/ggA_t_$MASS", "$BIN/ggA_t_$MASS_$SYSTEMATIC");
        cb.cp().channel({chn}).process({"ggH3_b"}).ExtractShapes(
          input_file_base, "$BIN/ggA_b_$MASS", "$BIN/ggA_b_$MASS_$SYSTEMATIC");
        cb.cp().channel({chn}).process({"ggH3_i"}).ExtractShapes(
          input_file_base, "$BIN/ggA_i_$MASS", "$BIN/ggA_i_$MASS_$SYSTEMATIC");

        if(enable_bsm_lowmass){
          cb.cp().channel({chn}).process({"ggH2_t_lowmass"}).ExtractShapes(
            input_file_base, "$BIN/ggH_t_$MASS", "$BIN/ggH_t_$MASS_$SYSTEMATIC");
          cb.cp().channel({chn}).process({"ggH2_b_lowmass"}).ExtractShapes(
            input_file_base, "$BIN/ggH_b_$MASS", "$BIN/ggH_b_$MASS_$SYSTEMATIC");
          cb.cp().channel({chn}).process({"ggH2_i_lowmass"}).ExtractShapes(
            input_file_base, "$BIN/ggH_i_$MASS", "$BIN/ggH_i_$MASS_$SYSTEMATIC");

          cb.cp().channel({chn}).process({"ggH3_t_lowmass"}).ExtractShapes(
            input_file_base, "$BIN/ggA_t_$MASS", "$BIN/ggA_t_$MASS_$SYSTEMATIC");
          cb.cp().channel({chn}).process({"ggH3_b_lowmass"}).ExtractShapes(
            input_file_base, "$BIN/ggA_b_$MASS", "$BIN/ggA_b_$MASS_$SYSTEMATIC");
          cb.cp().channel({chn}).process({"ggH3_i_lowmass"}).ExtractShapes(
            input_file_base, "$BIN/ggA_i_$MASS", "$BIN/ggA_i_$MASS_$SYSTEMATIC");
        }
      }
      // Inclusion of additional bbPhi is simple for the preconfigured process names
      // in mssm_bbH_signals_additional, reflecting the corresponding sub_analysis
      cb.cp().channel({chn}).process(mssm_bbH_signals_additional).ExtractShapes(
        input_file_base, "$BIN/bbH_$MASS", "$BIN/bbH_$MASS_$SYSTEMATIC");
      if(enable_bsm_lowmass){
        cb.cp().channel({chn}).process(mssm_bbH_lowmass_signals_additional).ExtractShapes(
          input_file_base, "$BIN/bbH_$MASS", "$BIN/bbH_$MASS_$SYSTEMATIC");
      }

      // In case full neutral Higgs modelling needs to be used (h, H, A or H1, H2, H3),
      // need to include additionally the SM-like Higgs boson. Configured process names:
      // sm-like-light sub_analysis: h
      // sm-like-heavy sub_analysis: H
      // cpv sub_analysis: H1
      // The way it is done for ggPhi and bbPhi depends in addition on whether sm125 or
      // bsm templates should be used for the SM-like signal.
      if(analysis == "bsm-model-dep-full"){
        // Here, the SUSY samples are taken for ggPhi and bbPhi
        if(sm_like_hists == "bsm"){
          // It stays simple for ggPhi, in case of sub_analysis sm-like-light or sm-like-heavy
          if(sub_analysis == "sm-like-light" || sub_analysis == "sm-like-heavy"){
            cb.cp().channel({chn}).process(mssm_ggH_signals_smlike).ExtractShapes(
              input_file_base, "$BIN/$PROCESS_$MASS", "$BIN/$PROCESS_$MASS_$SYSTEMATIC");
          }
          // Require explicit assignment for cpv sub_analysis for ggPhi due to different
          // template vs. process naming
          else if(sub_analysis == "cpv"){
            cb.cp().channel({chn}).process({"ggH1_t"}).ExtractShapes(
              input_file_base, "$BIN/ggh_t_$MASS", "$BIN/ggh_t_$MASS_$SYSTEMATIC");
            cb.cp().channel({chn}).process({"ggH1_b"}).ExtractShapes(
              input_file_base, "$BIN/ggh_b_$MASS", "$BIN/ggh_b_$MASS_$SYSTEMATIC");
            cb.cp().channel({chn}).process({"ggH1_i"}).ExtractShapes(
              input_file_base, "$BIN/ggh_i_$MASS", "$BIN/ggh_i_$MASS_$SYSTEMATIC");
          }
          // It stays simple for bbPhi, since using always the sample template
          cb.cp().channel({chn}).process(mssm_bbH_signals_smlike).ExtractShapes(
            input_file_base, "$BIN/bbH_$MASS", "$BIN/bbH_$MASS_$SYSTEMATIC");
        }
        // Here, the SM125 templates are used for ggPhi and bbPhi
        else if(sm_like_hists == "sm125"){
          cb.cp().channel({chn}).process(mssm_ggH_signals_smlike).ExtractShapes(
            input_file_base, "$BIN/ggH125$MASS", "$BIN/ggH125$MASS_$SYSTEMATIC");
          cb.cp().channel({chn}).process(mssm_bbH_signals_smlike).ExtractShapes(
            input_file_base, "$BIN/bbH_125$MASS", "$BIN/bbH_125$MASS_$SYSTEMATIC"); // Technically, still using SUSY sample but always the 125 GeV template
        }
        // Include qqPhi (and WPhi and ZPhi, if needed) always as 125 templates for
        // analysis bsm-model-dep-full
        cb.cp().channel({chn}).process(qqh_bsm_signals).ExtractShapes(
          input_file_base, "$BIN/qqH125$MASS", "$BIN/qqH125$MASS_$SYSTEMATIC");
        if(sm){
          cb.cp().channel({chn}).process(wh_bsm_signals).ExtractShapes(
            input_file_base, "$BIN/WH125$MASS", "$BIN/WH125$MASS_$SYSTEMATIC");
          cb.cp().channel({chn}).process(zh_bsm_signals).ExtractShapes(
            input_file_base, "$BIN/ZH125$MASS", "$BIN/ZH125$MASS_$SYSTEMATIC");
        }

        // Adding SM125 signal templates for SM hypothesis of analysis bsm-model-dep-full
        // These comprise ggH125, qqH125, bbH125, and in SM categories WH125 and ZH125
        if(hSM_treatment == "no-hSM-in-bg"){
          cb.cp().channel({chn}).process(ch::JoinStr({sm_signals, main_sm_signals})).process({"bbH125"}, false).ExtractShapes(
            input_file_base, "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC");
          cb.cp().channel({chn}).process({"bbH125"}).ExtractShapes(
            input_file_base, "$BIN/bbH_125$MASS", "$BIN/bbH_125$MASS_$SYSTEMATIC"); // "bbH125" needs special treatment because of template name spelling
        }
      }
    }
    if(variable=="m_sv_puppi" || variable=="m_sv_VS_pt_tt_splitpT" || lowmass) {
      cb.cp().channel({chn}).process(mssm_qqH_signals).ExtractShapes(
        input_file_base, "$BIN/qqH$MASS", "$BIN/qqH$MASS_$SYSTEMATIC");
      cb.cp().channel({chn}).process({"ggX_t"}).ExtractShapes(
        input_file_base, "$BIN/ggh_t_$MASS", "$BIN/ggh_t_$MASS_$SYSTEMATIC");
      cb.cp().channel({chn}).process({"ggX_b"}).ExtractShapes(
        input_file_base, "$BIN/ggh_b_$MASS", "$BIN/ggh_b_$MASS_$SYSTEMATIC");
      cb.cp().channel({chn}).process({"ggX_i"}).ExtractShapes(
        input_file_base, "$BIN/ggh_i_$MASS", "$BIN/ggh_i_$MASS_$SYSTEMATIC");

    }
  }

  // Rescale bbH125 to the right cross-section * BR (from 1pb to the value for 125.4 GeV)
  std::cout << "[INFO] Rescaling bbH125 process to correct cross-section * BR\n";
  boost::property_tree::ptree sm_preds;
  read_json(sm_predictions, sm_preds);
  std::cout << "\tCross section: " << float(sm_preds.get<float>("xs_bb_SMH125")) << std::endl;
  std::cout << "\tBranching fraction: " << float(sm_preds.get<float>("br_SMH125_tautau")) << std::endl;
  cb.cp().process({"bbH125"}).ForEachProc([&](ch::Process * proc) {
    proc->set_rate(proc->rate()*sm_preds.get<float>("xs_bb_SMH125")*sm_preds.get<float>("br_SMH125_tautau"));
    });


  // Manual rebinning for histograms
  if(manual_rebin)
  {
    std::map<std::string, std::map<unsigned int, std::map<unsigned int, std::vector<double> > > >binning_map;
    binning_map["em"] = {};
    binning_map["et"] = {};
    binning_map["mt"] = {};
    binning_map["tt"] = {};


    binning_map["em"][101] = {};
    binning_map["em"][102] = {};
    binning_map["em"][103] = {};
    binning_map["em"][104] = {};
    binning_map["em"][105] = {};
    binning_map["em"][106] = {};
    binning_map["em"][2][0] = {100., 200., 10.};
    binning_map["em"][2][1] = {200., 350., 25.};
    binning_map["em"][2][2] = {350., 500., 50.};
    binning_map["em"][2][3] = {500., 900., 100.};
    binning_map["em"][2][4] = {900., 1100., 200.};
    binning_map["em"][13] = {};
    binning_map["em"][14] = {};
    binning_map["em"][16] = {};
    binning_map["em"][20] = {};
    binning_map["em"][19] = {};

    binning_map["em"][32] = {};
    binning_map["em"][33] = {};
    binning_map["em"][34] = {};
    binning_map["em"][35] = {};
    binning_map["em"][36] = {};
    binning_map["em"][37] = {};

    binning_map["et"][101] = {};
    binning_map["et"][102] = {};
    binning_map["et"][103] = {};
    binning_map["et"][104] = {};
    binning_map["et"][105] = {};
    binning_map["et"][106] = {};
    binning_map["et"][13] = {};
    binning_map["et"][15] = {};
    binning_map["et"][16] = {};
    binning_map["et"][20] = {};
    binning_map["et"][21] = {};

    binning_map["et"][32] = {};
    binning_map["et"][33] = {};
    binning_map["et"][35] = {};
    binning_map["et"][36] = {};


    binning_map["mt"][101] = {};
    binning_map["mt"][102] = {};
    binning_map["mt"][103] = {};
    binning_map["mt"][104] = {};
    binning_map["mt"][105] = {};
    binning_map["mt"][106] = {};
    binning_map["mt"][13] = {};
    binning_map["mt"][15] = {};
    binning_map["mt"][16] = {};
    binning_map["mt"][20] = {};
    binning_map["mt"][21] = {};

    binning_map["mt"][32] = {};
    binning_map["mt"][33] = {};
    binning_map["mt"][35] = {};
    binning_map["mt"][36] = {};


    binning_map["tt"][101] = {};
    binning_map["tt"][102] = {};
    binning_map["tt"][103] = {};
    binning_map["tt"][104] = {};
    binning_map["tt"][105] = {};
    binning_map["tt"][106] = {};
    binning_map["tt"][10] = {};
    binning_map["tt"][16] = {};
    binning_map["tt"][20] = {};
    binning_map["tt"][21] = {};

    binning_map["tt"][32] = {};
    binning_map["tt"][35] = {};

    if(variable=="m_sv_VS_pt_tt_splitpT") {
      for(auto c : chns) {
        for(auto b : cb.cp().channel({c}).bin_id_set()) {
          binning_map[c][b][0]={0,60,60};
          binning_map[c][b][1]={60,200,10};
          binning_map[c][b][2]={200,260,20};
          binning_map[c][b][3]={260,300, 40};
        }
      }
    }

    if(variable=="m_sv_puppi") {
      // finer bins for nobtag categories 
      for(auto c : chns) {
        for(auto b : cb.cp().channel({c}).bin_id({mssm_nobtag_categories}).bin_id_set()){
          binning_map[c][b][0]={0,40,40};
          binning_map[c][b][1]={40,200,5};
          binning_map[c][b][2]={200,250,10};
          binning_map[c][b][3]={250,300, 25};
        }
      }
      //for btag categories use wider bins
      for(auto c : chns) {
        for(auto b : cb.cp().channel({c}).bin_id({mssm_nobtag_categories},false).bin_id_set()){
          binning_map[c][b][0]={0,60,60};
          binning_map[c][b][1]={60,80,20};
          binning_map[c][b][2]={80,120,10};
          binning_map[c][b][3]={120,200,20};
          binning_map[c][b][4]={200,240,40};
          binning_map[c][b][5]={240,300,60};
        }
      }
    }

    for(auto chn : chns)
    {
      for(auto b : cb.cp().channel({chn}).bin_id_set())
      {
        std::vector<double> binning = binning_from_map(binning_map[chn][b]);
        if(binning.size() > 0)
        {
            std::cout << "[INFO] Rebinning by hand for discriminator for bin: " << b << " in channel: " << chn << std::endl;
            cb.cp().channel({chn}).bin_id({b}).VariableRebin(binning);
        }
      }
    }
  }

  // Delete processes (other than mssm signals) with 0 yield
  std::cout << "[INFO] Filtering processes with null yield: \n";
  cb.FilterProcs([&](ch::Process *p) {
    if (std::find(mssm_signals.begin(), mssm_signals.end(), p->process()) != mssm_signals.end())
    {
      return false;
    }
    if (std::find(mssm_lowmass_signals.begin(), mssm_lowmass_signals.end(), p->process()) != mssm_lowmass_signals.end())
    {
      return false;
    }
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
    // For mssm signals: no action yet
    if (std::find(mssm_signals.begin(), mssm_signals.end(), s->process()) != mssm_signals.end())
    {
      return false;
    }
    if (std::find(mssm_lowmass_signals.begin(), mssm_lowmass_signals.end(), s->process()) != mssm_lowmass_signals.end())
    {
      return false;
    }
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

  // Special treatment for horizontally morphed mssm signals: Scale hists with negative intergral to zero, including its systematics
  // don't use this treatment for interference
  // Unless we use the prop_plot option in which case we just set the inteference to 0 if it goes negative 

  ch::CombineHarvester procs_no_i;
  if(prop_plot) procs_no_i = cb.cp();
  else procs_no_i = cb.cp().process({"ggH_i","ggh_i","ggA_i", "ggH1_i", "ggH2_i", "ggH3_i","ggH_i_lowmass","ggh_i_lowmass","ggA_i_lowmass", "ggH1_i_lowmass", "ggH2_i_lowmass", "ggH3_i_lowmass","ggX_i"}, false);

  if(prop_plot){
     // for prop_plot option if the inteference is negative we scale it positive and then add a rate parameter which will scale it negative again in the end 


     cb.cp().process({"ggH_i","ggh_i","ggA_i", "ggH1_i", "ggH2_i", "ggH3_i","ggH_i_lowmass","ggh_i_lowmass","ggA_i_lowmass", "ggH1_i_lowmass", "ggH2_i_lowmass", "ggH3_i_lowmass","ggX_i"}).ForEachProc([&](ch::Process *p) {
       if(p->rate() <= 0.0){
         std::cout << "[WARNING] Setting mssm inteference signal with negative yield to positive: \n ";
         std::cout << ch::Process::PrintHeader << *p << "\n";
         p->set_rate(p->rate()*-1.);

         cb.cp()
           .process({p->process()})
           .bin({p->bin()})
           .AddSyst(cb, "rate_minus","rateParam",SystMap<>::init(-1.0));
         cb.GetParameter("rate_minus")->set_range(-1.0,-1.0);
       }
     });
  }


  std::cout << "[INFO] Setting mssm signals with negative yield to 0 (excluding ggX interference).\n";
  procs_no_i.ForEachProc([mssm_signals,mssm_lowmass_signals](ch::Process *p) {
    if (std::find(mssm_signals.begin(), mssm_signals.end(), p->process()) != mssm_signals.end() || std::find(mssm_lowmass_signals.begin(), mssm_lowmass_signals.end(), p->process()) != mssm_lowmass_signals.end())
    {
      if(p->rate() <= 0.0){
        std::cout << "[WARNING] Setting mssm signal with negative yield to 0: \n ";
        std::cout << ch::Process::PrintHeader << *p << "\n";
        auto newhist = p->ClonedShape();
        newhist->Scale(0.0);
        p->set_shape(std::move(newhist), true);
      }
    }
  });

  procs_no_i.ForEachSyst([mssm_signals,mssm_lowmass_signals](ch::Systematic *s) {
    if (std::find(mssm_signals.begin(), mssm_signals.end(), s->process()) != mssm_signals.end() || std::find(mssm_lowmass_signals.begin(), mssm_lowmass_signals.end(), s->process()) != mssm_lowmass_signals.end())
    {
      if (s->type() == "shape") {
        if (s->shape_u()->Integral() <= 0.0 || s->shape_d()->Integral() <= 0.0) {
          std::cout << "[WARNING] Setting systematic for mssm signal with negative yield to 0: \n ";
          std::cout << ch::Systematic::PrintHeader << *s << "\n";
          auto newhist_u = s->ClonedShapeU();
          auto newhist_d = s->ClonedShapeD();
          newhist_u->Scale(0.0);
          newhist_d->Scale(0.0);
          s->set_shapes(std::move(newhist_u), std::move(newhist_d), nullptr);
        }
      }
    }
  });


  std::cout << "[Info] Summary of process yields: \n ";
  cb.cp().ForEachProc([&](ch::Process *p) {
    std::cout << ch::Process::PrintHeader << *p << "\n";
  });

  // Look for cases where a systematic changes the sign of the yield. These cases are due to statistical fluctuations so set the systematic shift to the nominal template
  // This is needed otherwise we get complaints about functions that evaluate as NaN
  cb.ForEachSyst([&](ch::Systematic *syst) {
  if (((syst->type().find("shape") != std::string::npos)
       && (syst->ClonedShapeU()->Integral()==0. || syst->ClonedShapeD()->Integral() == 0.)

       && (syst->process() == "bbH1" || syst->process() == "bbH2" || syst->process() == "bbH3" || syst->process() == "bbh" || syst->process() == "bbH" || syst->process() == "bbA" || syst->process() == "bbH1_lowmass" || syst->process() == "bbH2_lowmass" || syst->process() == "bbH3_lowmass" || syst->process() == "bbh_lowmass" || syst->process() == "bbH_lowmass" || syst->process() == "bbA_lowmass"
           || syst->process() == "ggH_i" || syst->process() == "ggh_i" || syst->process() == "ggA_i" || syst->process() == "ggH_i_lowmass" || syst->process() == "ggh_i_lowmass" || syst->process() == "ggA_i_lowmass"
           || syst->process() == "ggH1_i" || syst->process() == "ggH2_i" || syst->process() == "ggH3_i" || syst->process() == "ggH1_i_lowmass" || syst->process() == "ggH2_i_lowmass" || syst->process() == "ggH3_i_lowmass" || syst->process() == "ggX_i"))

      || ((syst->name().find("CMS_htt_boson_scale_met") != std::string::npos || syst->name().find("CMS_htt_boson_res_met") != std::string::npos
           || syst->name().find("CMS_scale_e") != std::string::npos || syst->name().find("CMS_scale_t_3prong_2018") != std::string::npos)

          && syst->ClonedShapeU()->Integral()==0 && syst->ClonedShapeD()->Integral() == 0

          && (syst->process() == "bbH1" || syst->process() == "bbH2" || syst->process() == "bbH3" || syst->process() == "bbh" || syst->process() == "bbH" || syst->process() == "bbA" || syst->process() == "bbH1_lowmass" || syst->process() == "bbH2_lowmass" || syst->process() == "bbH3_lowmass" || syst->process() == "bbh_lowmass"  || syst->process() == "bbH_lowmass" || syst->process() == "bbA_lowmass"))){

          std::cout << "Setting empty up and down templates to the nominal template \n";
          std::cout << ch::Systematic::PrintHeader << *syst << "\n";
          cb.cp().ForEachProc([&](ch::Process *proc){
          bool match_proc = (MatchingProcess(*proc,*syst));
          if(match_proc){
            if (proc->ClonedShape()->Integral() != 0){
              auto nominal = (TH1D*)proc->ClonedShape().get()->Clone();
              syst->set_value_u(1.0);
              syst->set_value_d(1.0);
              auto shape_u=(TH1D*)nominal->Clone();
              auto shape_d=(TH1D*)nominal->Clone();
              syst->set_shapes(std::unique_ptr<TH1>(static_cast<TH1*>(shape_u)),std::unique_ptr<TH1>(static_cast<TH1*>(shape_d)),nullptr);
            }
            }
          });
        }
    if (syst->type().find("lnN") != std::string::npos) {
      if(syst->value_u()<0.0) {
        std::cout << "[WARNING] Setting lnN systematic variation to the nominal as the yields would change sign otherwise \n ";
        std::cout << ch::Systematic::PrintHeader << *syst << "\n";
        syst->set_value_u(1.0);
      }
      if (syst->asymm() && syst->value_d()<0.0){
        std::cout << "[WARNING] Setting lnN systematic variation to the nominal as the yields would change sign otherwise \n ";
        std::cout << ch::Systematic::PrintHeader << *syst << "\n";
        syst->set_value_d(1.0);
      }
    }
    if (syst->type() =="shape") {
      double value_u = syst->value_u();
      double value_d = syst->value_d();

      if(value_u<=0. || value_d<=0.) {
        if(value_u<0. || value_d<0.) std::cout << "[WARNING] Setting shape systematic variation to the nominal as the yields would change sign otherwise \n ";
        else if(value_u==0. || value_d==0.) std::cout << "[WARNING] Setting shape systematic variation to the nominal as the yields would go to zero otherwise \n ";
        std::cout << ch::Systematic::PrintHeader << *syst << "\n";
        TH1D *shape_u = (TH1D*)syst->ClonedShapeU().get()->Clone();
        TH1D *shape_d = (TH1D*)syst->ClonedShapeD().get()->Clone();
        TH1D* nominal = new TH1D();
        cb.cp().ForEachProc([&](ch::Process *proc){
          bool match_proc = (MatchingProcess(*proc,*syst));
          if(match_proc) nominal = (TH1D*)proc->ClonedShape().get()->Clone();
        });
        if(value_u<=0.){
          syst->set_value_u(1.0);
          shape_u=(TH1D*)nominal->Clone();
        }
        if(value_d<=0.){
          syst->set_value_d(1.0);
          shape_d=(TH1D*)nominal->Clone();
        }
        syst->set_shapes(std::unique_ptr<TH1>(static_cast<TH1*>(shape_u)),std::unique_ptr<TH1>(static_cast<TH1*>(shape_d)),nullptr);
      }
    }
  });
  // rebinning of SM categories according to ML analysis == "  // Rebin categories to predefined binning for binning
  if (rebin_sm && sm) {
    // Rebin background categories
      for(auto b : cb.cp().bin_id_set())
      {
        TString bstr = b;
        TString cat = category;
        std::cout << "[INFO] Desciding the binning for " << b << "/" << category << std::endl;
        if (cat.Contains("xxh")){
          std::cout << "[INFO] Performing auto-rebinning for SM signal category.\n";
          auto rebin = ch::AutoRebin().SetBinThreshold(5.0).SetBinUncertFraction(0.9).SetRebinMode(1).SetPerformRebin(true).SetVerbosity(1);
          rebin.Rebin(cb, cb);
        }
        else {
          std::cout << "[INFO] Rebin background bin " << b << "\n";
          auto shape = cb.cp().bin_id({b}).backgrounds().GetShape();
          auto min = shape.GetBinLowEdge(1);
          std::vector<double> sm_binning = {min, 0.4, 0.5, 0.6, 0.7, 1.0};
          if(bstr.Contains("em") && bstr.Contains("misc")) sm_binning = {min, 0.4, 1.0};
          else if(bstr.Contains("em_emb")) sm_binning = {min, 0.4, 0.5, 0.6, 1.0};
          else if(bstr.Contains("et") && bstr.Contains("misc")) sm_binning = {min, 0.4, 0.5, 0.6, 1.0};
          else if(bstr.Contains("mt") && bstr.Contains("misc")) sm_binning = {min, 0.4, 0.5, 0.6, 1.0};
          else if(bstr.Contains("mt") && bstr.Contains("emb")) sm_binning = {min, 0.4, 0.5, 0.6, 1.0};
          else sm_binning = {min, 0.4, 0.5, 0.6, 0.7, 1.0};
          std::cout << "[INFO] Using binning: ";
          printVector(sm_binning);
          std::cout << "\n";
          cb.cp().bin_id({b}).VariableRebin(sm_binning);
        }
      }
  }

  std::vector<int> mssm_bins = {2,32,33,34,35,135,235,335,435,36,136,236,336,436,37,137,237,337,437,132,232,332,432,133,233,333,433,134,234,334,434};

  // Turn systematics into lnN
  std::cout << "[INFO] Transforming shape systematics for category " << category << std::endl;
  cb.cp().bin_id(mssm_bins, false).ForEachSyst([category, mssm_signals, mssm_lowmass_signals](ch::Systematic *s){
    TString sname = TString(s->name());
    if((s->type() == "shape") && (std::find(mssm_signals.begin(), mssm_signals.end(), s->process()) == mssm_signals.end() && std::find(mssm_lowmass_signals.begin(), mssm_lowmass_signals.end(), s->process()) == mssm_lowmass_signals.end() ))
    {
      double err_u = 0.0;
      double err_d = 0.0;
      int nbins = s->shape_u()->GetNbinsX();
      double yield_u = s->shape_u()->IntegralAndError(1,nbins,err_u);
      double yield_d = s->shape_d()->IntegralAndError(1,nbins,err_d);
      double value_u = s->value_u();
      double value_d = s->value_d();
      // std::cout << "value_u: " << value_u << " value_d: " << value_d << " unc: " << err_u/yield_u+err_d/yield_d << " diff: " << std::abs(value_u-1.0)+std::abs(value_d-1.0) << std::endl;
      if (std::abs(value_u-1.0)+std::abs(value_d-1.0)<err_u/yield_u+err_d/yield_d){
        s->set_type("lnN");
        std::cout << "\tSetting systematic " << s->name() << " of process " << s->process() << " to lnN" << std::endl;
      }
      else{
        std::cout << "\tLeaving systematic " << s->name() << " of process " << s->process() << " as shape" << std::endl;
      }
      // std::cout << ch::Systematic::PrintHeader << *s << "\n";
    }
  });

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
  for(auto u : jetmet_systs) ConvertShapesToLnN (cb.cp().bin_id(mssm_bins).process({"TTL","TTT"},false), u);
  // also convert ttbar in the tt channel
  for(auto u : jetmet_systs) ConvertShapesToLnN (cb.cp().bin_id(mssm_bins).channel({"tt"}).process({"TTL","TTT"}), u);
  // also convert ttbar in the nobtag categories, but only when fitting 1D variables

  // some FF unc1 systematics for the tt channel only affect the normalisations so can be converted to lnN:
  for (string y : {"2016","2017","2018"}) {
    ConvertShapesToLnN(cb.cp().bin_id({32,35,135,235,335,435,132,232,332,432}), "CMS_ff_total_qcd_stat_dR_unc1_tt_"+y);
    ConvertShapesToLnN(cb.cp().bin_id({32,35,135,235,335,435,132,232,332,432}), "CMS_ff_total_qcd_stat_pt_unc1_tt_"+y);
    ConvertShapesToLnN(cb.cp().bin_id({35,135,235,335,435}),    "CMS_ff_total_qcd_stat_njet1_jet_pt_low_unc1_tt_"+y);
    ConvertShapesToLnN(cb.cp().bin_id({35,135,235,335,435}),    "CMS_ff_total_qcd_stat_njet1_jet_pt_med_unc1_tt_"+y);
    ConvertShapesToLnN(cb.cp().bin_id({35,135,235,335,435}),    "CMS_ff_total_qcd_stat_njet1_jet_pt_high_unc1_tt_"+y);
  }

  // rename some fake factor systematics so that they are decorrelated between categories to match how closure corrections are measured
  for (string y : {"2016","2017","2018"}) {

    if(variable=="m_sv_puppi") {
      cb.cp().channel({"et"}).RenameSystematic(cb,"CMS_ff_total_ttbar_msv_shape_syst_mTtight_"+y,"CMS_ff_total_ttbar_msv_shape_syst_mTtight_et_"+y);
      cb.cp().channel({"et"}).RenameSystematic(cb,"CMS_ff_total_ttbar_msv_shape_syst_mTloose_"+y,"CMS_ff_total_ttbar_msv_shape_mTloose_et_"+y);
      cb.cp().channel({"mt"}).RenameSystematic(cb,"CMS_ff_total_ttbar_msv_shape_syst_mTtight_"+y,"CMS_ff_total_ttbar_msv_shape_syst_mTtight_mt_"+y);
      cb.cp().channel({"mt"}).RenameSystematic(cb,"CMS_ff_total_ttbar_msv_shape_syst_mTloose_"+y,"CMS_ff_total_ttbar_msv_shape_mTloose_mt_"+y);
    }

    for (string u : {"unc1", "unc2"}) {

      cb.cp().bin_id({32,132,232,332,432}).channel({"tt"}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_dR_"+u+"_tt_"+y,"CMS_ff_total_qcd_stat_dR_"+u+"_tt_Nbtag0_"+y);
      cb.cp().bin_id({35,135,235,335,435}).channel({"tt"}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_dR_"+u+"_tt_"+y,"CMS_ff_total_qcd_stat_dR_"+u+"_tt_NbtagGt1_"+y);

      cb.cp().bin_id({32,132,232,332,432}).channel({"tt"}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_pt_"+u+"_tt_"+y,"CMS_ff_total_qcd_stat_pt_"+u+"_tt_Nbtag0_"+y);
      cb.cp().bin_id({35,135,235,335,435}).channel({"tt"}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_pt_"+u+"_tt_"+y,"CMS_ff_total_qcd_stat_pt_"+u+"_tt_NbtagGt1_"+y);

      for (string c : {"mt","et"}) {
        cb.cp().bin_id({32,132,232,332,432}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_"+y,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_Nbtag0_MTLt40_"+y);
        cb.cp().bin_id({33,133,233,333,433}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_"+y,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_Nbtag0_MT40To70_"+y);
        cb.cp().bin_id({35,135,235,335,435}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_"+y,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_NbtagGt1_MTLt40_"+y);
        cb.cp().bin_id({36,136,236,336,436}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_"+y,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_NbtagGt1_MT40To70_"+y);

        cb.cp().bin_id({32,35,135,235,335,435,132,232,332,432}).channel({c}).RenameSystematic(cb,"CMS_ff_total_ttbar_stat_l_pt_"+u+"_"+c+"_"+y,"CMS_ff_total_ttbar_stat_l_pt_"+u+"_"+c+"_MTLt40_"+y);
        cb.cp().bin_id({33,36,136,236,336,436,133,233,333,433}).channel({c}).RenameSystematic(cb,"CMS_ff_total_ttbar_stat_l_pt_"+u+"_"+c+"_"+y,"CMS_ff_total_ttbar_stat_l_pt_"+u+"_"+c+"_MT40To70_"+y);

        cb.cp().bin_id({32,33,132,232,332,432,133,233,333,433}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_os_"+u+"_"+c+"_"+y,"CMS_ff_total_qcd_stat_os_"+u+"_"+c+"_Nbtag0_"+y);
        cb.cp().bin_id({35,135,235,335,435,36,136,236,336,436}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_os_"+u+"_"+c+"_"+y,"CMS_ff_total_qcd_stat_os_"+u+"_"+c+"_NbtagGt1_"+y);

      }
    }
  }

  // the following code is used to decorrelate wjets and qcd systematics by category
  // the btag and nobtag are always decorrelated
  // we also decorrelate the wjets by loose and tight mT since we are extrapolating to different mT regions
  // in cases where these uncertainties were not derived seperatly for Nbjets>0 (usually due to limited stats) we double the uncertainty
  // note this doubling is already done in the FF workspaces for the wjets_syst_extrap so we don't need to do it again here
  // we decorrelate the wjets extrapolation uncertainties by category

  for (string y : {"2016","2017","2018"}) {
    for (string c : {"mt","et"}) {
      cb.cp().bin_id({32,132,232,332,432}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_extrap_"+c+"_"+y,"CMS_ff_total_wjets_syst_extrap_"+c+"_Nbtag0_MTLt40_"+y);
      cb.cp().bin_id({33,133,233,333,433}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_extrap_"+c+"_"+y,"CMS_ff_total_wjets_syst_extrap_"+c+"_Nbtag0_MT40To70_"+y);
      cb.cp().bin_id({35,135,235,335,435}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_extrap_"+c+"_"+y,"CMS_ff_total_wjets_syst_extrap_"+c+"_NbtagGt1_MTLt40_"+y);
      cb.cp().bin_id({36,136,236,336,436}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_extrap_"+c+"_"+y,"CMS_ff_total_wjets_syst_extrap_"+c+"_NbtagGt1_MT40To70_"+y);

      cb.cp().bin_id({32,132,232,332,432}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_"+c+"_"+y,"CMS_ff_total_wjets_syst_"+c+"_Nbtag0_MTLt40_"+y);
      cb.cp().bin_id({33,133,233,333,433}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_"+c+"_"+y,"CMS_ff_total_wjets_syst_"+c+"_Nbtag0_MT40To70_"+y);
      cb.cp().bin_id({35,135,235,335,435}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_"+c+"_"+y,"CMS_ff_total_wjets_syst_"+c+"_NbtagGt1_MTLt40_"+y);
      cb.cp().bin_id({36,136,236,336,436}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_"+c+"_"+y,"CMS_ff_total_wjets_syst_"+c+"_NbtagGt1_MT40To70_"+y);
    }
  }
  for (string y : {"2016","2017","2018"}) {
    string c = "tt";
    cb.cp().bin_id({32,132,232,332,432}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_"+c+"_"+y,"CMS_ff_total_wjets_syst_"+c+"_Nbtag0_"+y);
    cb.cp().bin_id({35,135,235,335,435}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_syst_"+c+"_"+y,"CMS_ff_total_wjets_syst_"+c+"_NbtagGt1_"+y);

    cb.cp().bin_id({32,132,232,332,432}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_syst_dr_closure_"+c+"_"+y,"CMS_ff_total_qcd_syst_dr_closure_"+c+"_Nbtag0_"+y);
    cb.cp().bin_id({35,135,235,335,435}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_syst_dr_closure_"+c+"_"+y,"CMS_ff_total_qcd_syst_dr_closure_"+c+"_NbtagGt1_"+y);

    // scale Wjets uncertainty by 2 for tt channels in btag category.
    cb.cp().bin_id({35,135,235,335,435}).channel({c}).syst_name({"CMS_ff_total_wjets_syst_"+c+"_NbtagGt1_"+y}).ForEachSyst([&](ch::Systematic *syst) {
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
      cb.cp().bin_id({32,33,132,232,332,432,133,233,333,433}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_syst_"+c+"_"+y,"CMS_ff_total_qcd_syst_"+c+"_Nbtag0_"+y);
      cb.cp().bin_id({35,135,235,335,435,36,136,236,336,436}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_syst_"+c+"_"+y,"CMS_ff_total_qcd_syst_"+c+"_NbtagGt1_"+y);
      if(c != "tt") {
        cb.cp().bin_id({32,33,132,232,332,432,133,233,333,433}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_syst_iso_"+c+"_"+y,"CMS_ff_total_qcd_syst_iso_"+c+"_Nbtag0_"+y);
        cb.cp().bin_id({35,135,235,335,435,36,136,236,336,436}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_syst_iso_"+c+"_"+y,"CMS_ff_total_qcd_syst_iso_"+c+"_NbtagGt1_"+y);
        // scale QCD extrapolation uncertainties by 2 for et and mt channels in btag category.
        cb.cp().bin_id({35,135,235,335,435,36,136,236,336,436}).channel({c}).syst_name({"CMS_ff_total_qcd_syst_"+c+"_NbtagGt1_"+y,"CMS_ff_total_qcd_syst_iso_"+c+"_NbtagGt1_"+y}).ForEachSyst([&](ch::Systematic *syst) {
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

  for(auto u : met_uncerts) ConvertShapesToLnN (cb.cp().bin_id(mssm_bins).process({"ZTT"}, false), u);

  if(variable=="m_sv_puppi" || variable=="m_sv_VS_pt_tt_splitpT") {
    // convert ggH theory uncertainties to lnN when fitting m_sv
    std::vector<std::string> ggh_theory = {"Hdamp_ggH_t_REWEIGHT","Hdamp_ggH_b_REWEIGHT","Hdamp_ggH_i_REWEIGHT","QCDscale_ggH_REWEIGHT"}; 
    for(auto u : ggh_theory) ConvertShapesToLnN (cb.cp().bin_id(mssm_bins).process({"ZTT"}, false), u);
  }

  //// convert TER to lnN for btag category and most boosted nobtag categories
  //ConvertShapesToLnN (cb.cp().bin_id({132,232}, false), "CMS_res_t");
  //// split TER uncertainty by era:
  //for (string y : {"2016","2017","2018"}) cb.cp().era({y}).RenameSystematic(cb,"CMS_res_t","CMS_res_t_"+y);
  //for (string x : {"1prong", "1prong1pizero", "3prong", "3prong1pi0"}) {
  //  ConvertShapesToLnN (cb.cp().bin_id({132,232}, false), "CMS_res_t_"+x);
  //  for (string y : {"2016","2017","2018"}) cb.cp().era({y}).RenameSystematic(cb,"CMS_res_t_"+x,"CMS_res_t_"+x+"_"+y);
  //}

  // At this point we can fix the negative bins for the remaining processes
  // We don't want to do this for the ggH i component since this can have negative bins
  // Unless we use the prop_plot option in which case we just set the inteference to 0 if it goes negative 

  std::cout << "[INFO] Fixing negative bins.\n";
  procs_no_i.ForEachProc([](ch::Process *p) {
    if (ch::HasNegativeBins(p->shape())) {
      std::cout << "[WARNING] Fixing negative bins for process: \n ";
      std::cout << ch::Process::PrintHeader << *p << "\n";
      auto newhist = p->ClonedShape();
      ch::ZeroNegativeBins(newhist.get());
      p->set_shape(std::move(newhist), false);
    }
  });

  procs_no_i.ForEachSyst([](ch::Systematic *s) {
    if (s->type().find("shape") == std::string::npos)
      return;
    if (ch::HasNegativeBins(s->shape_u()) ||
        ch::HasNegativeBins(s->shape_d())) {
      std::cout << "[WARNING] Fixing negative bins for systematic: \n ";
      std::cout << ch::Systematic::PrintHeader << *s << "\n";
      auto newhist_u = s->ClonedShapeU();
      auto newhist_d = s->ClonedShapeD();
      ch::ZeroNegativeBins(newhist_u.get());
      ch::ZeroNegativeBins(newhist_d.get());
      s->set_shapes(std::move(newhist_u), std::move(newhist_d), nullptr);
    }
  });

  // Replacing observation with the sum of the backgrounds (Asimov data)
  // useful to be able to check this, so don't do the replacement
  // for these
  if (!real_data) {
    for (auto b : cb.cp().bin_set()) {
      std::cout << "[INFO] Replacing data with asimov in bin " << b << "\n";
      auto background_shape = cb.cp().bin({b}).backgrounds().GetShape();
      std::cout << "[INFO] Integral of background shape in bin " << b << ": " << background_shape.Integral()  << "\n";
      auto total_procs_shape = cb.cp().bin({b}).data().GetShape();
      total_procs_shape.Reset("M");
      std::cout << "[INFO] Integral of initial asimov data shape in bin " << b << ": " << total_procs_shape.Integral()  << "\n";
      // Desired Asimov model: BG + Higgs. Since ggH and qqH H->tautau treated as signal, remaining Higgs processes as BG. Need to access the signals() + bg shapes
      if(analysis == "sm"){
        auto signal_shape = cb.cp().bin({b}).signals().GetShape();
        std::cout << "[INFO] Integral of SM HTT shapes treated as signal  in bin " << b << ": " << signal_shape.Integral()  << "\n";
        bool no_signal = (signal_shape.GetNbinsX() == 1 && signal_shape.Integral() == 0.0);
        bool no_background = (background_shape.GetNbinsX() == 1 && background_shape.Integral() == 0.0);
        if(no_signal && no_background)
        {
          std::cout << "\t[WARNING] No signal and no background available in bin " << b << std::endl;
        }
        else if(no_background)
        {
          std::cout << "\t[WARNING] No background available in bin " << b << std::endl;
          total_procs_shape = total_procs_shape + signal_shape;
        }
        else if(no_signal)
        {
          std::cout << "\t[WARNING] No signal available in bin " << b << std::endl;
          total_procs_shape = total_procs_shape + background_shape;
        }
        else
        {
          std::cout << "[INFO] Setting sum of background and SM HTT signal shapes to asimov data in bin " << b << "\n";
          total_procs_shape = total_procs_shape + background_shape + signal_shape;
        }
      }
      // Desired Asimov model: BG( + Higgs). Since H->tautau treated all as background( if required), so it is sufficient to consider the bg shape
      else if(analysis == "bsm-model-indep" || analysis == "bsm-model-dep-additional" || (analysis == "bsm-model-dep-full" && hSM_treatment == "hSM-in-bg")){
        bool no_background = (background_shape.GetNbinsX() == 1 && background_shape.Integral() == 0.0);
        if(no_background)
        {
          std::cout << "\t[WARNING] No background available in bin " << b << std::endl;
        }
        else
        {
          std::cout << "[INFO] Setting background to asimov data in bin " << b << "\n";
          total_procs_shape = total_procs_shape + background_shape;
        }
      }
      // Desired Asimov model: BG + Higgs. Since H->tautau treated all as signal (together with mssm !!!), need to retrieve the SM H->tautau shapes & add it to the asimov dataset
      else if(analysis == "bsm-model-dep-full" and hSM_treatment == "no-hSM-in-bg"){
        auto sm_signal_shape = cb.cp().bin({b}).process(ch::JoinStr({sm_signals, main_sm_signals})).GetShape();
        std::cout << "[INFO] Integral of SM HTT signal shape in bin " << b << ": " << sm_signal_shape.Integral()  << "\n";
        bool no_background = (background_shape.GetNbinsX() == 1 && background_shape.Integral() == 0.0);
        bool no_sm_signal = (sm_signal_shape.GetNbinsX() == 1 && sm_signal_shape.Integral() == 0.0);
        if(no_sm_signal && no_background)
        {
          std::cout << "\t[WARNING] No signal and no background available in bin " << b << std::endl;
        }
        else if(no_background)
        {
          std::cout << "\t[WARNING] No background available in bin " << b << std::endl;
          total_procs_shape = total_procs_shape + sm_signal_shape;
        }
        else if(no_sm_signal)
        {
          std::cout << "\t[WARNING] No signal available in bin " << b << std::endl;
          total_procs_shape = total_procs_shape + background_shape;
        }
        else
        {
          std::cout << "[INFO] Setting sum of background and SM HTT signal shapes to asimov data in bin " << b << "\n";
          total_procs_shape = total_procs_shape + background_shape + sm_signal_shape;
        }
      }
      std::cout << "[INFO] Integral of final asimov data in bin " << b << ": " << total_procs_shape.Integral()  << "\n";
      cb.cp().bin({b}).ForEachObs([&](ch::Observation *obs) {
        obs->set_shape(total_procs_shape,true);
      });
    }
  }

  if (auto_rebin && !sm) {
    std::cout << "[INFO] Performing auto-rebinning.\n";
    auto rebin = ch::AutoRebin().SetBinThreshold(0.2).SetBinUncertFraction(0.9).SetRebinMode(1).SetPerformRebin(true).SetVerbosity(1);
    rebin.Rebin(cb, cb);
  }

  // This function modifies every entry to have a standardised bin name of
  // the form: {analysis}_{channel}_{bin_id}_{era}
  ch::SetStandardBinNames(cb, "$ANALYSIS_$CHANNEL_$BINID_$ERA");
  ch::CombineHarvester cb_obs = cb.deep().backgrounds();

  // Adding bin-by-bin uncertainties
  if (use_automc) {
    std::cout << "[INFO] Adding bin-by-bin uncertainties.\n";
    cb.SetAutoMCStats(cb, 0.);
  }
  // Setup morphed mssm signals for bsm analyses
  RooWorkspace ws("htt", "htt");

  std::map<std::string, RooAbsReal *> mass_var = {
    {"ggh_t", &mh}, {"ggh_b", &mh}, {"ggh_i", &mh},
    {"ggH_t", &mH}, {"ggH_b", &mH}, {"ggH_i", &mH},
    {"ggA_t", &mA}, {"ggA_b", &mA}, {"ggA_i", &mA},
    {"ggh_t_lowmass", &mh_lowmass}, {"ggh_b_lowmass", &mh_lowmass}, {"ggh_i_lowmass", &mh_lowmass},
    {"ggH_t_lowmass", &mH_lowmass}, {"ggH_b_lowmass", &mH_lowmass}, {"ggH_i_lowmass", &mH_lowmass},
    {"ggA_t_lowmass", &mA_lowmass}, {"ggA_b_lowmass", &mA_lowmass}, {"ggA_i_lowmass", &mA_lowmass},
    {"bbh", &mh}, {"bbh_lowmass", &mh_lowmass},
    {"bbH", &mH}, {"bbH_lowmass", &mH_lowmass},
    {"bbA", &mA}, {"bbA_lowmass", &mA_lowmass}
  };

  std::map<std::string, std::string> process_norm_map = {
    {"ggh_t", "prenorm"}, {"ggh_b", "prenorm"}, {"ggh_i", "prenorm"},
    {"ggH_t", "prenorm"}, {"ggH_b", "prenorm"}, {"ggH_i", "prenorm"},
    {"ggA_t", "prenorm"}, {"ggA_b", "prenorm"}, {"ggA_i", "prenorm"},
    {"ggh_t_lowmass", "prenorm"}, {"ggh_b_lowmass", "prenorm"}, {"ggh_i_lowmass", "prenorm"},
    {"ggH_t_lowmass", "prenorm"}, {"ggH_b_lowmass", "prenorm"}, {"ggH_i_lowmass", "prenorm"},
    {"ggA_t_lowmass", "prenorm"}, {"ggA_b_lowmass", "prenorm"}, {"ggA_i_lowmass", "prenorm"},
    {"bbh", "norm"}, {"bbh_lowmass", "norm"},
    {"bbH", "norm"}, {"bbH_lowmass", "norm"},
    {"bbA", "norm"}, {"bbA_lowmass", "norm"}
  };

  // Avoid morphing for the SM-like 'bbphi' process in case it is using the 125 GeV template.
  if(sm_like_hists == "sm125")
  {
    if(sub_analysis == "sm-like-light"){ mass_var.erase("bbh"); process_norm_map.erase("bbh");}
    else if(sub_analysis == "sm-like-heavy"){ mass_var.erase("bbH"); process_norm_map.erase("bbH");}
  }

  if(analysis == "bsm-model-indep")
  {
    mass_var = {
      {"ggh_t", &MH}, {"ggh_b", &MH}, {"ggh_i", &MH},
      {"bbh", &MH}
    };

    process_norm_map = {
      {"ggh_t", "prenorm"}, {"ggh_b", "prenorm"}, {"ggh_i", "prenorm"},
      {"bbh", "norm"}
    };

    if(variable=="m_sv_puppi" || variable=="m_sv_VS_pt_tt_splitpT" || lowmass) {
    process_norm_map = {
      {"ggh_t", "prenorm"}, {"ggh_b", "prenorm"}, {"ggh_i", "prenorm"},
      {"ggX_t", "prenorm"}, {"ggX_b", "prenorm"}, {"ggX_i", "prenorm"},
      {"bbh", "norm"}
    };
    }


    std::cout << "[INFO] Adding aditional terms for mssm ggh NLO reweighting.\n";
    // Assuming sm fractions of t, b and i contributions of 'ggh' in model-independent analysis
    TFile fractions_sm(sm_gg_fractions.c_str());
    std::cout << "[INFO] --> Loading WS: " << sm_gg_fractions.c_str() << std::endl;
    RooWorkspace *w_sm = (RooWorkspace*)fractions_sm.Get("w");
    if(do_morph && !prop_plot) {
      //w_sm->var("Yb_MSSM_h")->setVal(0.); // un-comment to remove bottom and top-bottom contirbutions (top only)
      //w_sm->var("Yt_MSSM_h")->setVal(0.); // un-comment to remove top and top-bottom contirbutions (bottom only)
      w_sm->var("mh")->SetName("MH");
      RooAbsReal *t_frac = w_sm->function("ggh_t_MSSM_frac");
      RooAbsReal *b_frac = w_sm->function("ggh_b_MSSM_frac");
      RooAbsReal *i_frac = w_sm->function("ggh_i_MSSM_frac");
      t_frac->SetName("ggh_t_frac");
      b_frac->SetName("ggh_b_frac");
      i_frac->SetName("ggh_i_frac");

      RooRealVar Yt("Yt_MSSM_h", "Yt_MSSM_h", 1., -1., 1.);
      Yt.setConstant(true);
      RooRealVar Yb("Yb_MSSM_h", "Yb_MSSM_h", 1., 0., 1.);
      Yb.setConstant(true);

      ws.import(MH);
      ws.import(Yt);
      ws.import(Yb);
      ws.import(*t_frac, RooFit::RecycleConflictNodes());
      ws.import(*b_frac, RooFit::RecycleConflictNodes());
      ws.import(*i_frac, RooFit::RecycleConflictNodes());
    }
    else if (prop_plot) {
      w_sm->var("mh")->setVal(100.);
      RooAbsReal *t_frac = w_sm->function("ggh_t_MSSM_frac");
      RooAbsReal *b_frac = w_sm->function("ggh_b_MSSM_frac");
      RooAbsReal *i_frac = w_sm->function("ggh_i_MSSM_frac");
      t_frac->SetName("ggh_t_frac");
      b_frac->SetName("ggh_b_frac");
      i_frac->SetName("ggh_i_frac");
      ws.import(*t_frac, RooFit::RecycleConflictNodes());
      ws.import(*b_frac, RooFit::RecycleConflictNodes());
      ws.import(*i_frac, RooFit::RecycleConflictNodes());
    }
    else {
      w_sm->var("mh")->setVal(std::stof(non_morphed_mass));
      RooAbsReal *t_frac = w_sm->function("ggh_t_MSSM_frac");
      RooAbsReal *b_frac = w_sm->function("ggh_b_MSSM_frac");
      RooAbsReal *i_frac = w_sm->function("ggh_i_MSSM_frac");
      t_frac->SetName("ggh_t_frac");
      b_frac->SetName("ggh_b_frac");
      i_frac->SetName("ggh_i_frac");
      ws.import(*t_frac, RooFit::RecycleConflictNodes());
      ws.import(*b_frac, RooFit::RecycleConflictNodes());
      ws.import(*i_frac, RooFit::RecycleConflictNodes());
    }
    fractions_sm.Close();
  }

  if(sub_analysis == "cpv")
  {
    mass_var = {
      {"ggH1_t", &mH1}, {"ggH1_b", &mH1}, {"ggH1_i", &mH1},
      {"ggH2_t", &mH2}, {"ggH2_b", &mH2}, {"ggH2_i", &mH2},
      {"ggH3_t", &mH3}, {"ggH3_b", &mH3}, {"ggH3_i", &mH3},
      {"ggH1_t_lowmass", &mH1_lowmass}, {"ggH1_b_lowmass", &mH1_lowmass}, {"ggH1_i_lowmass", &mH1_lowmass},
      {"ggH2_t_lowmass", &mH2_lowmass}, {"ggH2_b_lowmass", &mH2_lowmass}, {"ggH2_i_lowmass", &mH2_lowmass},
      {"ggH3_t_lowmass", &mH3_lowmass}, {"ggH3_b_lowmass", &mH3_lowmass}, {"ggH3_i_lowmass", &mH3_lowmass},
      {"bbH1", &mH1}, {"bbH1_lowmass", &mH1_lowmass},
      {"bbH2", &mH2}, {"bbH2_lowmass", &mH2_lowmass},
      {"bbH3", &mH3}, {"bbH3_lowmass", &mH3_lowmass}
    };

    process_norm_map = {
      {"ggH1_t", "prenorm"}, {"ggH1_b", "prenorm"}, {"ggH1_i", "prenorm"},
      {"ggH2_t", "prenorm"}, {"ggH2_b", "prenorm"}, {"ggH2_i", "prenorm"},
      {"ggH3_t", "prenorm"}, {"ggH3_b", "prenorm"}, {"ggH3_i", "prenorm"},
      {"ggH1_t_lowmass", "prenorm"}, {"ggH1_b_lowmass", "prenorm"}, {"ggH1_i_lowmass", "prenorm"},
      {"ggH2_t_lowmass", "prenorm"}, {"ggH2_b_lowmass", "prenorm"}, {"ggH2_i_lowmass", "prenorm"},
      {"ggH3_t_lowmass", "prenorm"}, {"ggH3_b_lowmass", "prenorm"}, {"ggH3_i_lowmass", "prenorm"},
      {"bbH1", "norm"}, {"bbH1_lowmass", "norm"},
      {"bbH2", "norm"}, {"bbH2_lowmass", "norm"},
      {"bbH3", "norm"}, {"bbH3_lowmass", "norm"}
    };

    // Avoid morphing for the SM-like 'bbphi' process in case it is using the 125 GeV template.
    if(sm_like_hists == "sm125")
    {
      mass_var.erase("bbH1");
      process_norm_map.erase("bbH1");
    }
  }

  dout("[INFO] Prepare demo.");
  if(do_morph && !prop_plot && (analysis == "bsm-model-indep" || analysis == "bsm-model-dep-additional" || analysis == "bsm-model-dep-full"))
  {
    //TFile morphing_demo(("htt_mssm_morphing_" + category+ "_"  + era_tag + "_" + analysis + "_demo.root").c_str(), "RECREATE");

    // Perform morphing
    std::cout << "[INFO] Performing template morphing for mssm ggh and bbh.\n";
    auto morphFactory = ch::CMSHistFuncFactory();
    morphFactory.SetHorizontalMorphingVariable(mass_var);
    morphFactory.Run(cb, ws, process_norm_map);

    if(analysis == "bsm-model-indep"){
      // Adding 'norm' terms into workspace according to desired signals
      for (auto bin : cb.cp().bin_set())
      {
        for (auto proc : ch::JoinStr({mssm_ggH_signals,ggX_signals}))
        {
          std::string prenorm_name = bin + "_" + proc + "_morph_prenorm";
          std::string norm_name = bin + "_" + proc + "_morph_norm";
          std::string proc_frac = boost::replace_all_copy(proc,"ggX","ggh");
          ws.factory(TString::Format("expr::%s('@0*@1',%s, %s_frac)", norm_name.c_str(), prenorm_name.c_str(), proc_frac.c_str()));
        }
      }
    }

    // Saving workspace with morphed signals
    //morphing_demo.cd();
    //if (verbose)
    //  ws.Print();
    //ws.Write();
    //morphing_demo.Close();
    cb.AddWorkspace(ws);
    cb.ExtractPdfs(cb, "htt", "$BIN_$PROCESS_morph");
    cb.ExtractData("htt", "$BIN_data_obs");
    std::cout << "[INFO] Finished template morphing for mssm ggh and bbh.\n";
  }
  else if(!do_morph && analysis == "bsm-model-indep" && !prop_plot){

   // TODO: for high masses, this makes only a little difference, but why required? Problem with negative value below?
   //double Tfrac = ws.function("ggh_t_frac")->getVal();
   //double Bfrac = ws.function("ggh_b_frac")->getVal();
   //double Ifrac = ws.function("ggh_i_frac")->getVal();
   double Tfrac=1., Bfrac=0., Ifrac=0.; // use t-only when no morphing option is used
   if (Ifrac<0.) {
     Ifrac=fabs(Ifrac);
     // set a constant rate parameter = -1 for the interference
     cb.cp()
      .process({"ggh_i"})
      .AddSyst(cb, "rate_minus","rateParam",SystMap<>::init(-1.0));
     cb.GetParameter("rate_minus")->set_range(-1.0,-1.0);
   }
   std::cout << "setting fractions as t,b,i = " << Tfrac << "," << Bfrac << "," << Ifrac << std::endl;

   cb.cp().process({"ggh_t"}).ForEachProc([&](ch::Process * proc) {
     proc->set_rate(proc->rate()*Tfrac);
    });

   cb.cp().process({"ggh_b"}).ForEachProc([&](ch::Process * proc) {
     proc->set_rate(proc->rate()*Bfrac);
    });

   cb.cp().process({"ggh_i"}).ForEachProc([&](ch::Process * proc) {
     proc->set_rate(proc->rate()*Ifrac);
    });
  }
  else if(prop_plot){

   // TODO: for high masses, this makes only a little difference, but why required? Problem with negative value below?
   double Tfrac = ws.function("ggh_t_frac")->getVal();
   double Bfrac = ws.function("ggh_b_frac")->getVal();
   double Ifrac = ws.function("ggh_i_frac")->getVal();
   if (Ifrac<0.) {
     Ifrac=fabs(Ifrac);
     // set a constant rate parameter = -1 for the interference
     cb.cp()
      .process({"ggh_i"})
      .AddSyst(cb, "rate_minus_overall","rateParam",SystMap<>::init(-1.0));
     cb.GetParameter("rate_minus_overall")->set_range(-1.0,-1.0);
   }
   std::cout << "setting fractions as t,b,i = " << Tfrac << "," << Bfrac << "," << Ifrac << std::endl;

   cb.cp().process({"ggh_t"}).ForEachProc([&](ch::Process * proc) {
     proc->set_rate(proc->rate()*Tfrac);
    });

   cb.cp().process({"ggh_b"}).ForEachProc([&](ch::Process * proc) {
     proc->set_rate(proc->rate()*Bfrac);
    });

   cb.cp().process({"ggh_i"}).ForEachProc([&](ch::Process * proc) {
     proc->set_rate(proc->rate()*Ifrac);
    });
  }
  ch::CombineHarvester cb_cp;
  if(prop_plot) cb_cp = cb.deep();

  std::cout << "[INFO] Writing datacards to " << output_folder << std::endl;
  // We need to do this to make sure the ttbarShape uncertainty is added properly when we use a shapeU
  if(!prop_plot){
    cb.GetParameter("CMS_htt_ttbarShape")->set_range(-1.0,1.0);
    cb.GetParameter("CMS_htt_ttbarShape")->set_err_d(-1.);
    cb.GetParameter("CMS_htt_ttbarShape")->set_err_u(1.);
  } 

  // Decide, how to write out the datacards depending on --category option
  if(category == "all") {
      // Write out datacards. Naming convention important for rest of workflow. We
      // make one directory per chn-cat, one per chn and cmb. In this code we only
      // store the individual datacards for each directory to be combined later.
      ch::CardWriter writer(output_folder + "/" + era_tag + "/$TAG/$BIN.txt",
                            output_folder + "/" + era_tag + "/$TAG/common/htt_input_" + era_tag + ".root");

      // We're not using mass as an identifier - which we need to tell the
      // CardWriter
      // otherwise it will see "*" as the mass value for every object and skip it
      writer.SetWildcardMasses({});

      // Set verbosity
      if (verbose)
        writer.SetVerbosity(1);

      // Write datacards combined and per channel
      writer.WriteCards("cmb", cb);

      for (auto chn : chns) {
        writer.WriteCards(chn, cb.cp().channel({chn}));
      }
  }
  else {
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

      ch::CardWriter writer_restore(output_folder + "/restore_binning/$BIN.txt",
                                    output_folder + "/restore_binning/common/$BIN_input.root");

      // We're not using mass as an identifier - which we need to tell the
      // CardWriter
      // otherwise it will see "*" as the mass value for every object and skip it
      writer_restore.SetWildcardMasses({});

      // Set verbosity
      if (verbose)
        writer_restore.SetVerbosity(1);

      writer_restore.WriteCards("", cb_obs);
  }

  if (verbose)
    cb.PrintAll();

  std::cout << "[INFO] Done producing datacards.\n";



  if (prop_plot) {

    auto bin_set = cb.cp().bin_id({132,232,332,432,133,233,333,433}).bin_set();

    double new_bins[19] = {0.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 220.0, 240.0, 260.0, 300.0};
    TH1F proto("proto", "proto", 18, new_bins);

    for (auto const& bin : bin_set) {
      ch::CombineHarvester cmb_bin = std::move(cb_cp.cp().bin({bin}));
      TH1F bkg = cmb_bin.cp().backgrounds().GetShape();
      TH1F sig = cmb_bin.cp().signals().process({"ggh_t","ggh_b"}).GetShape(); // just use top only for getting the signal yield otherwise they wont have all been scaled by proper fractions
      TH1F sig_i = cmb_bin.cp().signals().process({"ggh_i"}).GetShape(); // need to get inteference seperatly then scale by negative sign if we scaled this positive previously
       
      double sig_scale=5.686; // set to best fit value of signal strength for ggH
      sig.Scale(sig_scale);
      sig_i.Scale(sig_scale);

      // calculate S/(S+B) between 70-150 GeV
      double s_sb = 0.;
      int bin_lo = sig.FindBin(70.);
      int bin_hi = sig.FindBin(150.) -1;
      double i_sig_i = sig_i.Integral(bin_lo,bin_hi);
      if(cmb_bin.GetParameter("rate_minus_overall")!=nullptr) i_sig_i*=-1.;
      if(cmb_bin.GetParameter("rate_minus")!=nullptr) i_sig_i*=-1.;
      //double i_sig = sig.Integral(bin_lo,bin_hi)/1.0807788; // the factor scales down by t fraction as we base weights on t-only
      double i_sig = sig.Integral(bin_lo,bin_hi)+i_sig_i;
      double i_bkg = bkg.Integral(bin_lo,bin_hi);
      //std::cout <<"----  " <<  bin << "  " << bin_lo << "    " << bin_hi << "  " << bkg.Integral(bin_lo,bin_hi) << "    " << sig.Integral(bin_lo,bin_hi) << "  " << i_sig_i << std::endl;  
      if(i_sig>0. || i_bkg>0.) s_sb = i_sig/(i_sig+i_bkg);
      //std::cout << "weight = " << s_sb << std::endl;

      cmb_bin.ForEachObs([&](ch::Observation *e) {
        TH1 const * old_h = e->shape();
        TH1F * new_h_ = (TH1F*)old_h->Clone();
        new_h_ = (TH1F*)new_h_->Rebin(18,"",new_bins);
        double wt = s_sb;
        new_h_->Scale(wt);
        std::unique_ptr<TH1> new_h = ch::make_unique<TH1F>(*(new_h_));
        
        new_h->Scale(e->rate());
        e->set_shape(std::move(new_h), true);
      });

      cmb_bin.ForEachProc([&](ch::Process *e) {
        TH1 const * old_h = e->shape();
        TH1F * new_h_ = (TH1F*)old_h->Clone();
        new_h_ = (TH1F*)new_h_->Rebin(18,"",new_bins);
        double wt = s_sb;
        new_h_->Scale(wt);
        std::unique_ptr<TH1> new_h = ch::make_unique<TH1F>(*(new_h_));
        new_h->Scale(e->rate());
        e->set_shape(std::move(new_h), true);
      });

      cmb_bin.ForEachSyst([&](ch::Systematic *e) {
        if (e->type() != "shape") return;

        TH1D *old_hu = (TH1D*)e->ClonedShapeU().get()->Clone();
        TH1D *old_hd = (TH1D*)e->ClonedShapeD().get()->Clone();

        float val_u = e->value_u(); 
        float val_d = e->value_d();
        double wt = s_sb;

        TH1F * new_hu_ = (TH1F*)old_hu->Clone();
        new_hu_ = (TH1F*)new_hu_->Rebin(18,"",new_bins);
        new_hu_->Scale(wt*val_u);
        std::unique_ptr<TH1> new_hu = ch::make_unique<TH1F>(*(new_hu_));

        TH1F * new_hd_ = (TH1F*)old_hd->Clone();
        new_hd_ = (TH1F*)new_hd_->Rebin(18,"",new_bins);
        new_hd_->Scale(wt*val_d);
        std::unique_ptr<TH1> new_hd = ch::make_unique<TH1F>(*(new_hd_));

        e->set_shapes(std::move(new_hu), std::move(new_hd), nullptr);
      });

    }
      ch::CardWriter writer(output_folder + "/" + era_tag + "/$TAG/$BIN.txt",
                            output_folder + "/" + era_tag + "/$TAG/common/$BIN_input.root");
      writer.WriteCards("plot_cmb", cb_cp.cp().bin_id({132,232,332,432,133,233,333,433}));
      writer.WriteCards("plot_0to50", cb_cp.cp().bin_id({132,133}));
      writer.WriteCards("plot_50to100", cb_cp.cp().bin_id({232,233}));
      writer.WriteCards("plot_100to200", cb_cp.cp().bin_id({332,333}));
      writer.WriteCards("plot_GT200", cb_cp.cp().bin_id({432,433}));
  }
 

}
