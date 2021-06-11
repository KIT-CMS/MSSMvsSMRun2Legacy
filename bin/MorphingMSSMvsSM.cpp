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
  using ch::JoinStr;
  using ch::syst::SystMap;

  // Define program options
  string output_folder = "output_MSSMvsSM_Run2";
  string base_path = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/MSSMvsSMRun2Legacy/shapes/";
  string sm_gg_fractions = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root";
  string chan = "mt";
  string category = "mt_nobtag_lowmsv_0jet_tightmt";
  string variable = "m_sv_puppi";

  bool auto_rebin = false;
  bool manual_rebin = false;
  bool real_data = false;
  bool verbose = true;
  bool use_automc = true;
  bool mva(false), no_emb(false);
  bool sm = false;
  bool rebin_sm = true;
  bool no_shape_systs = false;

  vector<string> mass_susy_ggH({}), mass_susy_qqH({}), parser_bkgs({}), parser_bkgs_em({}), parser_sm_signals({}), parser_main_sm_signals({});

  string analysis = "bsm-model-indep"; // "sm",  "bsm-model-indep", "bsm-model-dep-full", "bsm-model-dep-additional"
  string sub_analysis = "hSM-in-bg"; // analysis = "bsm-model-indep": "hSM-in-bg", "no-hSM-in-bg"; case with analysis = "bsm-model-dep-{full,additional}": "sm-like-light", "sm-like-heavy", "cpv"
  string sm_like_hists = "sm125"; // used in analysis = "bsm-model-dep-full": "sm125", "bsm"
  string categorization = "classic"; // "with-sm-ml", "sm-ml-only", "classic"
  int era = 2016; // 2016, 2017 or 2018
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()
      ("base-path,base_path", po::value<string>(&base_path)->default_value(base_path), "inputs, expected to contain a subdirectory <era>/<channel>")
      ("sm_gg_fractions", po::value<string>(&sm_gg_fractions)->default_value(sm_gg_fractions))
      ("channel", po::value<string>(&chan)->default_value(chan), "single channel to process")
      ("category", po::value<string>(&category)->default_value(category))
      ("variable", po::value<string>(&variable)->default_value(variable))
      ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(auto_rebin))
      ("manual_rebin", po::value<bool>(&manual_rebin)->default_value(manual_rebin))
      ("real_data", po::value<bool>(&real_data)->default_value(real_data))
      ("use_automc", po::value<bool>(&use_automc)->default_value(use_automc))
      ("verbose", po::value<bool>(&verbose)->default_value(verbose))
      ("output_folder", po::value<string>(&output_folder)->default_value(output_folder))
      ("analysis", po::value<string>(&analysis)->default_value(analysis))
      ("sub-analysis", po::value<string>(&sub_analysis)->default_value(sub_analysis))
      ("sm-like-hists", po::value<string>(&sm_like_hists)->default_value(sm_like_hists))
      ("categorization", po::value<string>(&categorization)->default_value(categorization))
      ("era", po::value<int>(&era)->default_value(era))
      ("no-emb,no-emb,no_emb", po::bool_switch(&no_emb), "use MC samples instead of embedding")
      ("debug,d", po::bool_switch(&debug), "debug printout")
      ("mva", po::bool_switch(&mva), "mva tau id is used")
      ("sm", po::value<bool>(&sm)->default_value(sm))
      ("mass-susy-ggH,mass_susy_ggH", po::value<vector<string>>(&mass_susy_ggH)->multitoken(), "mass_susy_ggH")
      ("mass-susy-qqH,mass_susy_qqH", po::value<vector<string>>(&mass_susy_qqH)->multitoken(), "mass_susy_qqH")
      ("bkgs", po::value<vector<string>>(&parser_bkgs)->multitoken(), "backgrounds")
      ("bkgs_em", po::value<vector<string>>(&parser_bkgs_em)->multitoken(), "backgrounds-em")
      ("sm_signals", po::value<vector<string>>(&parser_sm_signals)->multitoken(), "sm_signals")
      ("main_sm_signals", po::value<vector<string>>(&parser_main_sm_signals)->multitoken(), "main_sm_signals")
      ("no_shape_systs", po::value<bool>(&no_shape_systs)->default_value(no_shape_systs))
      ("help", "produce help message");
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);
  if (vm.count("help"))
  {
      cout << config << "\n";
      return 0;
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
  VString bkgs, bkgs_em, bkgs_tt, bkgs_HWW, sm_signals, main_sm_signals;
  VString mssm_ggH_signals, mssm_ggH_signals_additional, mssm_ggH_signals_smlike, mssm_ggH_signals_scalar, mssm_ggH_signals_pseudoscalar;
  VString mssm_bbH_signals, mssm_bbH_signals_additional, mssm_bbH_signals_smlike, mssm_bbH_signals_scalar, mssm_bbH_signals_pseudoscalar;
  VString mssm_signals, qqh_bsm_signals;
  if (sm == true){
  	sm_signals = {"WH125", "ZH125", "ttH125"};
  }
  else{
  	sm_signals = {};
  }
  main_sm_signals = {"ggH125", "qqH125"}; // qqH125 for mt,et,tt contains VBF+VH
  update_vector_by_byparser(sm_signals, parser_sm_signals, "sm_signals");
  update_vector_by_byparser(main_sm_signals, parser_main_sm_signals, "main_sm_signals");

  if(analysis == "bsm-model-indep")
  {
    mssm_ggH_signals = {"ggh_t", "ggh_b", "ggh_i"};
    mssm_bbH_signals = {"bbh"};
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
    }
    mssm_ggH_signals_additional = ch::JoinStr({mssm_ggH_signals_scalar, mssm_ggH_signals_pseudoscalar});
    mssm_bbH_signals_additional = ch::JoinStr({mssm_bbH_signals_scalar, mssm_bbH_signals_pseudoscalar});
    mssm_ggH_signals = ch::JoinStr({mssm_ggH_signals_smlike, mssm_ggH_signals_scalar, mssm_ggH_signals_pseudoscalar});
    mssm_bbH_signals = ch::JoinStr({mssm_bbH_signals_smlike, mssm_bbH_signals_scalar, mssm_bbH_signals_pseudoscalar});
  }
  mssm_signals = ch::JoinStr({mssm_ggH_signals, mssm_bbH_signals});


  std::cout << "Used BSM signals: ";
  for(auto proc : mssm_signals){
    std::cout << proc << " ";
  }
  std::cout << std::endl;

  bkgs = {"EMB", "ZL", "TTL", "VVL", "jetFakes"};
  bkgs_tt = {"EMB", "ZL", "TTL", "VVL", "jetFakes", "wFakes"};
  bkgs_HWW = {"ggHWW125", "qqHWW125", "WHWW125", "ZHWW125"};
  bkgs_em = {"EMB", "W", "QCD", "ZL", "TTL", "VVL"};
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
    dout("WARNING: the EMB process is removed from backgrounds");
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "EMB"), bkgs.end());
    bkgs_em.erase(std::remove(bkgs_em.begin(), bkgs_em.end(), "EMB"), bkgs_em.end());
  }
  map<int, VString> SUSYggH_masses;
  map<int, VString> SUSYbbH_masses;

  //bool do_morph=false;
  bool do_morph=true; // TODO use an option for that
  if(do_morph) {

    // new DESY datacards should have all masses now?
    SUSYbbH_masses[2018] = {"60","80","100","120","125","130","140","160","180","200","250","300","350","400","450","500","600","700","800","900","1000","1200","1400","1600","1800","2000","2300","2600","2900","3200","3500"};
    SUSYbbH_masses[2017] = SUSYbbH_masses[2018];
    SUSYbbH_masses[2016] = {"60","80","100","120","125","130","140","160","180","200","250","350","400","450","500","600","800","900","1200","1400","1600","1800","2000","2300","2600","2900","3200","3500"};  // Missing 300,700,1000
    SUSYggH_masses[2018] = {"60","80","100","120","125","130","140","160","180","200","250","300","350","400","450","500","600","700","800","900","1000","1200","1400","1600","1800","2000","2300","2600","2900","3200","3500"};
    SUSYggH_masses[2016] = SUSYggH_masses[2018];
    SUSYggH_masses[2017] = SUSYggH_masses[2018];

  } else {
    // dont use mass morphing - need to specify a mass here
    string mass_str = "700";
    SUSYggH_masses[2016] = {mass_str};
    SUSYggH_masses[2017] = {mass_str};
    SUSYggH_masses[2018] = {mass_str};
    SUSYbbH_masses[2016] = {mass_str};
    SUSYbbH_masses[2017] = {mass_str};
    SUSYbbH_masses[2018] = {mass_str};
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
    for(auto chn : chns){
        bkg_procs[chn] = JoinStr({bkg_procs[chn],sm_signals});
    }
  }
  else if((analysis == "bsm-model-indep" && sub_analysis == "hSM-in-bg") || analysis == "bsm-model-dep-additional"){
    bkg_procs["tt"] = JoinStr({bkg_procs["tt"],main_sm_signals});
    bkg_procs["mt"] = JoinStr({bkg_procs["mt"],main_sm_signals});
    bkg_procs["et"] = JoinStr({bkg_procs["et"],main_sm_signals});
    bkg_procs["em"] = JoinStr({bkg_procs["em"],main_sm_signals,bkgs_HWW});
    if(category == "et_xxh" || category == "et_tt" || category == "et_zll" || category == "et_misc" || category == "et_emb" || category == "et_ff"){
      bkg_procs["et"] = JoinStr({bkg_procs["et"],sm_signals,main_sm_signals,bkgs_HWW});
    }
    else if(category == "mt_xxh" || category == "mt_tt" || category == "mt_zll" || category == "mt_misc" || category == "mt_emb" || category == "mt_ff"){
      bkg_procs["mt"] = JoinStr({bkg_procs["mt"],sm_signals,main_sm_signals,bkgs_HWW});
    }
    else if(category == "tt_xxh" || category == "tt_misc" || category == "tt_emb" || category == "tt_ff"){
      bkg_procs["tt"] = JoinStr({bkg_procs["tt"],sm_signals,main_sm_signals,bkgs_HWW});
    }
    else if (category == "em_xxh" || category == "em_tt" || category == "em_ss" || category == "em_misc" || category == "em_db" || category == "em_emb") {
      bkg_procs["em"] = JoinStr({bkg_procs["em"], sm_signals, main_sm_signals,bkgs_HWW});
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
  // mA is used as model parameter in case of sub_analysis "sm-like-light", for "sm-like-heavy" and "cpv" it is mHp
  if(sub_analysis == "sm-like-light")
  {
    mA.setConstant(true);
  }

  // Define MSSM CPV model-dependent mass parameters mH3, mH2, mH1 (sub_analysis "cpv")
  RooRealVar mH3("mH3", "mH3", 125., 90., 4000.);
  RooRealVar mH2("mH2", "mH2", 125., 90., 4000.);
  RooRealVar mH1("mH1", "mH1", 125., 90., 4000.);

  // Define MSSM model-independent mass parameter MH
  RooRealVar MH("MH", "MH", 125., 90., 4000.);
  MH.setConstant(true);

  // Define categories
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
        {  2, "em_DZetaLtm35"},
        { 32, "em_Nbtag0_DZetaGt30"},
        { 33, "em_Nbtag0_DZetam10To30"},
        { 34, "em_Nbtag0_DZetam35Tom10"},
        { 35, "em_NbtagGt1_DZetaGt30"},
        { 36, "em_NbtagGt1_DZetam10To30"},
        { 37, "em_NbtagGt1_DZetam35Tom10"},
    };
  }
  else if(categorization == "sm-ml-only"){
    cats["et"] = {
      { 1, "et_xxh"}, // SM Signal Category

      {13, "et_tt"},
      {15, "et_zll"},
      {16, "et_misc"},
      {20, "et_emb"},
      {21, "et_ff"}
    };

    cats["mt"] = {
      { 1, "mt_xxh"}, // SM Signal Category

      {13, "mt_tt"},
      {15, "mt_zll"},
      {16, "mt_misc"},
      {20, "mt_emb"},
      {21, "mt_ff"}
    };

    cats["tt"] = {
      { 1, "tt_xxh"}, // SM Signal Category

      {16, "tt_misc"},
      {20, "tt_emb"},
      {21, "tt_ff"}
    };

    cats["em"] = {
      { 1, "em_xxh"}, // SM Signal Category

      {13, "em_tt"},
      {14, "em_ss"},
      {16, "em_misc"},
      {19, "em_db"},
      {20, "em_emb"}
    };
  }
  else if(categorization == "with-sm-ml"){
    cats["et"] = {
        { 1, "et_xxh"}, // SM Signal Category

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
    cats["mt"] = {
        { 1, "mt_xxh"}, // SM Signal Category

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
    cats["tt"] = {
        { 1, "tt_xxh"}, // SM Signal Category

        {16, "tt_misc"},
        {20, "tt_emb"},
        {21, "tt_ff"},

        {32, "tt_Nbtag0_MHGt250"},

        {35, "tt_NbtagGt1"},
    };
    cats["em"] = {
        { 1, "em_xxh"}, // SM Signal Category

        { 2, "em_DZetaLtm35"},

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
  }
  else throw std::runtime_error("Given categorization is not known.");

  // Create combine harverster object
  ch::CombineHarvester cb;
  cb.SetFlag("workspaces-use-clone", true);


  // Introduce ordering of categories for the final discriminator in MSSM
  std::vector<int> sm_categories = {13,14,15,16,19,20,21}; // Control regions from the ML SM HTT analysis
  std::vector<int> em_control_category = {2}; // Control region for em channel
  std::vector<int> mssm_btag_categories = {35,36,37}; // b-tagged MSSM-like categories with mt_tot as discriminator
  std::vector<int> mssm_nobtag_categories = {32,33,34}; // non-btagged MSSM-like categories with mt_tot as discriminator
  std::vector<int> sm_signal_category = {1}; // category for the SM signal


  for (auto chn : chns) {
  // build category maps used for the different analyses
    Categories sm_and_btag_cats = cats[chn]; // contain 1, 2, 13-21, 35-37
    Categories mssm_btag_cats = cats[chn]; // contain 2, 35-37
    Categories mssm_cats = cats[chn]; // contain 2, 32-37
    Categories sm_signal_cat = cats[chn]; // contain 1

    for (auto catit = sm_signal_cat.begin(); catit != sm_signal_cat.end(); ++catit)
    {
      if(std::find(sm_categories.begin(), sm_categories.end(), (*catit).first) != sm_categories.end()){
        sm_signal_cat.erase(catit);
        --catit;
      }
      if(std::find(mssm_nobtag_categories.begin(), mssm_nobtag_categories.end(), (*catit).first) != mssm_nobtag_categories.end()){
        sm_signal_cat.erase(catit);
        --catit;
      }
      if(std::find(mssm_btag_categories.begin(), mssm_btag_categories.end(), (*catit).first) != mssm_btag_categories.end()){
        sm_signal_cat.erase(catit);
        --catit;
      }
      if(std::find(em_control_category.begin(), em_control_category.end(), (*catit).first) != em_control_category.end()){
        sm_signal_cat.erase(catit);
        --catit;
      }
    }

    for (auto catit = mssm_cats.begin(); catit != mssm_cats.end(); ++catit)
    {
      if(std::find(sm_categories.begin(), sm_categories.end(), (*catit).first) != sm_categories.end()){
        mssm_cats.erase(catit);
        --catit;
      }
      if(std::find(sm_signal_category.begin(), sm_signal_category.end(), (*catit).first) != sm_signal_category.end()){
        mssm_cats.erase(catit);
        --catit;
      }
    }

    for (auto catit = mssm_btag_cats.begin(); catit != mssm_btag_cats.end(); ++catit)
    {
      if(std::find(mssm_nobtag_categories.begin(), mssm_nobtag_categories.end(), (*catit).first) != mssm_nobtag_categories.end()){
        mssm_btag_cats.erase(catit);
        --catit;
      }
      if(std::find(sm_categories.begin(), sm_categories.end(), (*catit).first) != sm_categories.end()){
        mssm_btag_cats.erase(catit);
        --catit;
      }
      if(std::find(sm_signal_category.begin(), sm_signal_category.end(), (*catit).first) != sm_signal_category.end()){
        mssm_btag_cats.erase(catit);
        --catit;
      }
    }

    for (auto catit = sm_and_btag_cats.begin(); catit != sm_and_btag_cats.end(); ++catit)
    {
      if(std::find(mssm_nobtag_categories.begin(), mssm_nobtag_categories.end(), (*catit).first) != mssm_nobtag_categories.end()){
        sm_and_btag_cats.erase(catit);
        --catit;
      }
    }

    std::cout << "[INFO] Using the following categories:" << std::endl;
    std::cout << "   sm_and_btag_cats:" << std::endl;
    for (const auto i: sm_and_btag_cats)
      std::cout << "      " << i.first << ' ' << i.second << std::endl;
    std::cout  << std::endl;
    std::cout << "    mssm_cats:" << std::endl;
    for (const auto i: mssm_cats)
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

    cb.AddObservations({"*"}, {"htt"}, {era_tag}, {chn}, cats[chn]);
    cb.AddProcesses({"*"}, {"htt"}, {era_tag}, {chn}, bkg_procs[chn], cats[chn], false);

    if(analysis == "sm"){
      cb.AddProcesses({""}, {"htt"}, {era_tag}, {chn}, main_sm_signals, cats[chn], true);
    }
    else if(analysis == "bsm-model-indep"){
      // filter masses above treshold:
      // std::vector<std::string> masses_SM_bbH;
      // std::vector<std::string> masses_SM_ggH;
      // int bbH_threshold = SM_thresholds_bbH[era][chn];
      // int ggH_threshold = SM_thresholds_ggH[era][chn];
      // if(sm){
      //     MH.setRange(90., std::max(bbH_threshold,ggH_threshold));
      // }
      // std::copy_if (SUSYbbH_masses[era].begin(), SUSYbbH_masses[era].end(), std::back_inserter(masses_SM_bbH), [bbH_threshold](std::string i){return std::stoi(i)<=bbH_threshold;} );
      // std::copy_if (SUSYggH_masses[era].begin(), SUSYggH_masses[era].end(), std::back_inserter(masses_SM_ggH), [ggH_threshold](std::string i){return std::stoi(i)<=ggH_threshold;} );
      // std::cout << "[INFO] Using masses for SM: ";
      //   // printVector(masses_SM_bbH);
      //   // printVector(masses_SM_ggH);
      // cb.AddProcesses(masses_SM_bbH, {"htt"}, {era_tag}, {chn}, mssm_bbH_signals, sm_signal_cat, true);
      // cb.AddProcesses(masses_SM_ggH, {"htt"}, {era_tag}, {chn}, mssm_ggH_signals, sm_signal_cat, true);

      // Adding configured SUSY signals in all categories but SM ML HTT background categories (13-21) for bsm model-independent analyses
      cb.AddProcesses(SUSYbbH_masses[era], {"htt"}, {era_tag}, {chn}, mssm_bbH_signals, sm_signal_cat, true);
      cb.AddProcesses(SUSYggH_masses[era], {"htt"}, {era_tag}, {chn}, mssm_ggH_signals, sm_signal_cat, true);
      cb.AddProcesses(SUSYbbH_masses[era], {"htt"}, {era_tag}, {chn}, mssm_bbH_signals, mssm_cats, true);
      cb.AddProcesses(SUSYggH_masses[era], {"htt"}, {era_tag}, {chn}, mssm_ggH_signals, mssm_cats, true);
    }
    else if(analysis == "bsm-model-dep-additional" || analysis == "bsm-model-dep-full"){
      // Adding at first the additional Higgs boson signals
      Categories additional_higgses_cats;
      if(categorization == "classic") // classic categories cover full mass range for additional Higgs signals
      {
         additional_higgses_cats = cats[chn];
      }
      else if(categorization == "with-sm-ml" || categorization == "sm-ml-only")
      {
         if(sub_analysis == "sm-like-light" || sub_analysis == "cpv") // in that case, the additional Higgs bosons are relatively heavy
         {
         additional_higgses_cats = mssm_cats;
         }
         else if(sub_analysis == "sm-like-heavy") // in that case, all additional Higgs bosons are relatively light (< 200 GeV) --> don't consider them in high mass no-btag categories
         {
         additional_higgses_cats = sm_and_btag_cats;
         }
      }
      cb.AddProcesses(SUSYggH_masses[era], {"htt"}, {era_tag}, {chn}, mssm_ggH_signals_additional, additional_higgses_cats, true);
      cb.AddProcesses(SUSYbbH_masses[era], {"htt"}, {era_tag}, {chn}, mssm_bbH_signals_additional, additional_higgses_cats, true);

      if(analysis == "bsm-model-dep-full")
      {
        // Adding SM Higgs processes as signal for model-dependent analyses with full neutral Higgs modelling (since testing then against SM Higgs + BG hypothesis)
        cb.AddProcesses({""}, {"htt"}, {era_tag}, {chn}, ch::JoinStr({main_sm_signals, sm_signals}), cats[chn], true);

        // Defining categories for qqphi & ggphi: exclude mssm_nobtag_categories in case SM ML HTT categories are used because of m_sv >= 250 GeV cut
        // Defining categories for bbphi: exclude mssm_nobtag_categories and SM categories in case SM ML HTT categories are used because of m_sv >= 250 GeV cut and bad population in SM cats
        Categories qq_gg_phi_cats, bbphi_cats;
        if(categorization == "classic")
        {
          qq_gg_phi_cats = cats[chn];
          bbphi_cats = cats[chn];
        }
        else if(categorization == "with-sm-ml" || categorization == "sm-ml-only")
        {
          qq_gg_phi_cats = sm_and_btag_cats;
          bbphi_cats = mssm_btag_cats;
        }

        // Adding the qqphi process for all bsm model-dependent analyses with full neutral Higgs modelling
        cb.AddProcesses({""}, {"htt"}, {era_tag}, {chn}, qqh_bsm_signals, qq_gg_phi_cats, true);

        VString empty_masses = {""};
        VString ggH_SMlike_masses = (sm_like_hists == "sm125") ? empty_masses : SUSYggH_masses[era];
        cb.AddProcesses(ggH_SMlike_masses, {"htt"}, {era_tag}, {chn}, mssm_ggH_signals_smlike, qq_gg_phi_cats, true);
        VString bbH_SMlike_masses = (sm_like_hists == "sm125") ? empty_masses : SUSYbbH_masses[era];
        cb.AddProcesses(bbH_SMlike_masses, {"htt"}, {era_tag}, {chn}, mssm_bbH_signals_smlike, bbphi_cats, true);
      }
    }
  }

  // Add systematics TODO: update also for all BSM flavours (ggphi, bbphi, qqphi)
  dout("[INFO] Add systematics AddMSSMvsSMRun2Systematics, embedding:", ! no_emb, " sm categories:", sm);
  ch::AddMSSMvsSMRun2Systematics(cb, true, ! no_emb, true, true, true, era, mva, sm);
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
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
      input_file_base, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");

    if(analysis == "sm"){
      cb.cp().channel({chn}).process(main_sm_signals).ExtractShapes(
        input_file_base, "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC");
    }
    else if(analysis == "bsm-model-indep"){
      cb.cp().channel({chn}).process(mssm_ggH_signals).ExtractShapes(
        input_file_base, "$BIN/$PROCESS_$MASS", "$BIN/$PROCESS_$MASS_$SYSTEMATIC");
      cb.cp().channel({chn}).process(mssm_bbH_signals).ExtractShapes(
        input_file_base, "$BIN/bbH_$MASS", "$BIN/bbH_$MASS_$SYSTEMATIC");
    }
    else if(analysis == "bsm-model-dep-additional" || analysis == "bsm-model-dep-full"){
      if(sub_analysis == "sm-like-light" || sub_analysis == "sm-like-heavy"){
        cb.cp().channel({chn}).process(mssm_ggH_signals_additional).ExtractShapes(
          input_file_base, "$BIN/$PROCESS_$MASS", "$BIN/$PROCESS_$MASS_$SYSTEMATIC");
      }
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
      }
      cb.cp().channel({chn}).process(mssm_bbH_signals_additional).ExtractShapes(
        input_file_base, "$BIN/bbH_$MASS", "$BIN/bbH_$MASS_$SYSTEMATIC");

      if(analysis == "bsm-model-dep-full"){
        if(sm_like_hists == "bsm"){
          if(sub_analysis == "sm-like-light" || sub_analysis == "sm-like-heavy"){
            cb.cp().channel({chn}).process(mssm_ggH_signals_smlike).ExtractShapes(
              input_file_base, "$BIN/$PROCESS_$MASS", "$BIN/$PROCESS_$MASS_$SYSTEMATIC");
          }
          else if(sub_analysis == "cpv"){
            cb.cp().channel({chn}).process({"ggH1_t"}).ExtractShapes(
              input_file_base, "$BIN/ggh_t_$MASS", "$BIN/ggh_t_$MASS_$SYSTEMATIC");
            cb.cp().channel({chn}).process({"ggH1_b"}).ExtractShapes(
              input_file_base, "$BIN/ggh_b_$MASS", "$BIN/ggh_b_$MASS_$SYSTEMATIC");
            cb.cp().channel({chn}).process({"ggH1_i"}).ExtractShapes(
              input_file_base, "$BIN/ggh_i_$MASS", "$BIN/ggh_i_$MASS_$SYSTEMATIC");
          }
          cb.cp().channel({chn}).process(mssm_bbH_signals_smlike).ExtractShapes(
            input_file_base, "$BIN/bbH_$MASS", "$BIN/bbH_$MASS_$SYSTEMATIC");
        }
        else if(sm_like_hists == "sm125"){
          cb.cp().channel({chn}).process(mssm_ggH_signals_smlike).ExtractShapes(
            input_file_base, "$BIN/ggH125$MASS", "$BIN/ggH125$MASS_$SYSTEMATIC");
          cb.cp().channel({chn}).process(mssm_bbH_signals_smlike).ExtractShapes(
            input_file_base, "$BIN/bbH_125$MASS", "$BIN/bbH_125$MASS_$SYSTEMATIC");
        }
        cb.cp().channel({chn}).process(qqh_bsm_signals).ExtractShapes(
          input_file_base, "$BIN/qqH125$MASS", "$BIN/qqH125$MASS_$SYSTEMATIC");
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
  std::cout << "[INFO] Setting mssm signals with negative yield to 0 (excluding ggX interference).\n";
  cb.cp().process({"ggH_i","ggh_i","ggA_i", "ggH1_i", "ggH2_i", "ggH3_i"}, false).ForEachProc([mssm_signals](ch::Process *p) {
    if (std::find(mssm_signals.begin(), mssm_signals.end(), p->process()) != mssm_signals.end())
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

  cb.cp().process({"ggH_i","ggh_i","ggA_i", "ggH1_i", "ggH2_i", "ggH3_i"}, false).ForEachSyst([mssm_signals](ch::Systematic *s) {
    if (std::find(mssm_signals.begin(), mssm_signals.end(), s->process()) != mssm_signals.end())
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
  if (((syst->type().find("shape") != std::string::npos) &&
       (syst->ClonedShapeU()->Integral()==0. || syst->ClonedShapeD()->Integral() == 0.) &&

       (syst->process() == "bbH2" || syst->process() == "bbH3" || syst->process() == "bbH" || syst->process() == "bbA" ||
        syst->process() == "ggH_i" || syst->process() == "ggh_i" || syst->process() == "ggA_i" ||
        syst->process() == "ggH1_i" || syst->process() == "ggH2_i" || syst->process() == "ggH3_i")) ||

      ((syst->name().find("CMS_htt_boson_scale_met") != std::string::npos || syst->name().find("CMS_htt_boson_res_met") != std::string::npos ||
        syst->name().find("CMS_scale_e") != std::string::npos || syst->name().find("CMS_scale_t_3prong_2018") != std::string::npos) &&

       syst->ClonedShapeU()->Integral()==0 && syst->ClonedShapeD()->Integral() == 0 &&

       (syst->process() == "bbH2" || syst->process() == "bbH3" || syst->process() == "bbH" || (syst->process() == "bbA")))){
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

      if(value_u <0 || value_d <0) {
        std::cout << "[WARNING] Setting shape systematic variation to the nominal as the yields would change sign otherwise \n ";
        std::cout << ch::Systematic::PrintHeader << *syst << "\n";
        TH1D *shape_u = (TH1D*)syst->ClonedShapeU().get()->Clone();
        TH1D *shape_d = (TH1D*)syst->ClonedShapeD().get()->Clone();
        TH1D* nominal = new TH1D();
        cb.cp().ForEachProc([&](ch::Process *proc){
          bool match_proc = (MatchingProcess(*proc,*syst));
          if(match_proc) nominal = (TH1D*)proc->ClonedShape().get()->Clone();
        });
        if(value_u<0){
          syst->set_value_u(1.0);
          shape_u=(TH1D*)nominal->Clone();
        }
        if(value_d<0){
          syst->set_value_d(1.0);
          shape_d=(TH1D*)nominal->Clone();
        }
        syst->set_shapes(std::unique_ptr<TH1>(static_cast<TH1*>(shape_u)),std::unique_ptr<TH1>(static_cast<TH1*>(shape_d)),nullptr);
      }
    }
  });

  // Manual rebinning for histograms
  if(manual_rebin)
  {
    std::map<std::string, std::map<unsigned int, std::map<unsigned int, std::vector<double> > > >binning_map;
    binning_map["em"] = {};
    binning_map["et"] = {};
    binning_map["mt"] = {};
    binning_map["tt"] = {};


    binning_map["em"][1] = {};
    binning_map["em"][2] = {};
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

    binning_map["et"][1] = {};
    binning_map["et"][13] = {};
    binning_map["et"][15] = {};
    binning_map["et"][16] = {};
    binning_map["et"][20] = {};
    binning_map["et"][21] = {};

    binning_map["et"][32] = {};
    binning_map["et"][33] = {};
    binning_map["et"][35] = {};
    binning_map["et"][36] = {};


    binning_map["mt"][1] = {};
    binning_map["mt"][13] = {};
    binning_map["mt"][15] = {};
    binning_map["mt"][16] = {};
    binning_map["mt"][20] = {};
    binning_map["mt"][21] = {};

    binning_map["mt"][32] = {};
    binning_map["mt"][33] = {};
    binning_map["mt"][35] = {};
    binning_map["mt"][36] = {};


    binning_map["tt"][1] = {};
    binning_map["tt"][10] = {};
    binning_map["tt"][16] = {};
    binning_map["tt"][20] = {};
    binning_map["tt"][21] = {};

    binning_map["tt"][32] = {};
    binning_map["tt"][35] = {};

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
  // rebinning of SM categories according to ML analysis == "  // Rebin categories to predefined binning for binning
  if (rebin_sm && sm) {
    // Rebin background categories
    // for(auto chn : chns)
    // {
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
      // }
    }
  }

  std::vector<int> mssm_bins = {2,32,33,34,35,36,37};

  // Turn systematics into lnN
  std::cout << "[INFO] Transforming shape systematics for category " << category << std::endl;
  cb.cp().bin_id(mssm_bins, false).ForEachSyst([category, mssm_signals](ch::Systematic *s){
    TString sname = TString(s->name());
    if((s->type() == "shape") && (std::find(mssm_signals.begin(), mssm_signals.end(), s->process()) == mssm_signals.end()))
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

  for(auto u : jetmet_systs) ConvertShapesToLnN (cb.cp().bin_id(mssm_bins), u);

  // some FF unc1 systematics for the tt channel only affect the normalisations so can be converted to lnN:
  for (string y : {"2016","2017","2018"}) {
    ConvertShapesToLnN(cb.cp().bin_id({32,35}), "CMS_ff_total_qcd_stat_dR_unc1_tt_"+y);
    ConvertShapesToLnN(cb.cp().bin_id({32,35}), "CMS_ff_total_qcd_stat_pt_unc1_tt_"+y);
    ConvertShapesToLnN(cb.cp().bin_id({35}),    "CMS_ff_total_qcd_stat_njet1_jet_pt_low_unc1_tt_"+y);
    ConvertShapesToLnN(cb.cp().bin_id({35}),    "CMS_ff_total_qcd_stat_njet1_jet_pt_med_unc1_tt_"+y);
    ConvertShapesToLnN(cb.cp().bin_id({35}),    "CMS_ff_total_qcd_stat_njet1_jet_pt_high_unc1_tt_"+y);
  }

  // rename some fake factor systematics so that they are decorrelated between categories to match how closure corrections are measured
  for (string y : {"2016","2017","2018"}) {
    for (string u : {"unc1", "unc2"}) {

      cb.cp().bin_id({32}).channel({"tt"}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_dR_"+u+"_tt_"+y,"CMS_ff_total_qcd_stat_dR_"+u+"_tt_Nbtag0_"+y);
      cb.cp().bin_id({35}).channel({"tt"}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_dR_"+u+"_tt_"+y,"CMS_ff_total_qcd_stat_dR_"+u+"_tt_NbtagGt1_"+y);

      cb.cp().bin_id({32}).channel({"tt"}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_pt_"+u+"_tt_"+y,"CMS_ff_total_qcd_stat_pt_"+u+"_tt_Nbtag0_"+y);
      cb.cp().bin_id({35}).channel({"tt"}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_pt_"+u+"_tt_"+y,"CMS_ff_total_qcd_stat_pt_"+u+"_tt_NbtagGt1_"+y);

      for (string c : {"mt","et"}) {
        cb.cp().bin_id({32}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_"+y,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_Nbtag0_MTLt40_"+y);
        cb.cp().bin_id({33}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_"+y,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_Nbtag0_MT40To70_"+y);
        cb.cp().bin_id({35}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_"+y,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_NbtagGt1_MTLt40_"+y);
        cb.cp().bin_id({36}).channel({c}).RenameSystematic(cb,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_"+y,"CMS_ff_total_wjets_stat_extrap_"+u+"_"+c+"_NbtagGt1_MT40To70_"+y);
        
        cb.cp().bin_id({32,35}).channel({c}).RenameSystematic(cb,"CMS_ff_total_ttbar_stat_l_pt_"+u+"_"+c+"_"+y,"CMS_ff_total_ttbar_stat_l_pt_"+u+"_"+c+"_MTLt40_"+y);
        cb.cp().bin_id({33,36}).channel({c}).RenameSystematic(cb,"CMS_ff_total_ttbar_stat_l_pt_"+u+"_"+c+"_"+y,"CMS_ff_total_ttbar_stat_l_pt_"+u+"_"+c+"_MT40To70_"+y);

        cb.cp().bin_id({32,33}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_os_"+u+"_"+c+"_"+y,"CMS_ff_total_qcd_stat_os_"+u+"_"+c+"_Nbtag0_"+y);
        cb.cp().bin_id({35,36}).channel({c}).RenameSystematic(cb,"CMS_ff_total_qcd_stat_os_"+u+"_"+c+"_"+y,"CMS_ff_total_qcd_stat_os_"+u+"_"+c+"_NbtagGt1_"+y);

      }
    }
  }

  // At this point we can fix the negative bins for the remaining processes
  // We don't want to do this for the ggH i component since this can have negative bins
  std::cout << "[INFO] Fixing negative bins.\n";
  cb.cp().process({"ggH_i","ggh_i","ggA_i", "ggH1_i", "ggH2_i", "ggH3_i"}, false).ForEachProc([](ch::Process *p) {
    if (ch::HasNegativeBins(p->shape())) {
      std::cout << "[WARNING] Fixing negative bins for process: \n ";
      std::cout << ch::Process::PrintHeader << *p << "\n";
      auto newhist = p->ClonedShape();
      ch::ZeroNegativeBins(newhist.get());
      p->set_shape(std::move(newhist), false);
    }
  });

  cb.cp().process({"ggH_i","ggh_i","ggA_i", "ggH1_i", "ggH2_i", "ggH3_i"}, false).ForEachSyst([](ch::Systematic *s) {
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
      else if(analysis == "bsm-model-indep" || analysis == "bsm-model-dep-additional"){
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
      else if(analysis == "bsm-model-dep-full"){
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
    {"bbh", &mh},
    {"bbH", &mH},
    {"bbA", &mA}
  };

  std::map<std::string, std::string> process_norm_map = {
    {"ggh_t", "prenorm"}, {"ggh_b", "prenorm"}, {"ggh_i", "prenorm"},
    {"ggH_t", "prenorm"}, {"ggH_b", "prenorm"}, {"ggH_i", "prenorm"},
    {"ggA_t", "prenorm"}, {"ggA_b", "prenorm"}, {"ggA_i", "prenorm"},
    {"bbh", "norm"},
    {"bbH", "norm"},
    {"bbA", "norm"}
  };

  // Avoid morphing for the SM-like 'bbphi' process in case it is using the 125 GeV template.
  if(sm_like_hists == "smh125")
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


   
    std::cout << "[INFO] Adding aditional terms for mssm ggh NLO reweighting.\n";
    // Assuming sm fractions of t, b and i contributions of 'ggh' in model-independent analysis
    TFile fractions_sm(sm_gg_fractions.c_str());
    std::cout << "[INFO] --> Loading WS: " << sm_gg_fractions.c_str() << std::endl;
    RooWorkspace *w_sm = (RooWorkspace*)fractions_sm.Get("w");
    if(do_morph) {
      //w_sm->var("Yb_MSSM_h")->setVal(0.); // un-comment to remove bottom and top-bottom contirbutions (top only)
      //w_sm->var("Yt_MSSM_h")->setVal(0.); // un-comment to remove top and top-bottom contirbutions (bottom only)
      w_sm->var("mh")->SetName("MH");
      RooAbsReal *t_frac = w_sm->function("ggh_t_MSSM_frac");
      RooAbsReal *b_frac = w_sm->function("ggh_b_MSSM_frac");
      RooAbsReal *i_frac = w_sm->function("ggh_i_MSSM_frac");
      t_frac->SetName("ggh_t_frac");
      b_frac->SetName("ggh_b_frac");
      i_frac->SetName("ggh_i_frac");
      ws.import(MH);
      ws.import(*t_frac, RooFit::RecycleConflictNodes());
      ws.import(*b_frac, RooFit::RecycleConflictNodes());
      ws.import(*i_frac, RooFit::RecycleConflictNodes());
    }
    else{
      w_sm->var("mh")->setVal(std::stof(SUSYggH_masses[2018][0]));
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
      {"bbH1", &mH1},
      {"bbH2", &mH2},
      {"bbH3", &mH3}
    };

    process_norm_map = {
      {"ggH1_t", "prenorm"}, {"ggH1_b", "prenorm"}, {"ggH1_i", "prenorm"},
      {"ggH2_t", "prenorm"}, {"ggH2_b", "prenorm"}, {"ggH2_i", "prenorm"},
      {"ggH3_t", "prenorm"}, {"ggH3_b", "prenorm"}, {"ggH3_i", "prenorm"},
      {"bbH1", "norm"},
      {"bbH2", "norm"},
      {"bbH3", "norm"}
    };

    // Avoid morphing for the SM-like 'bbphi' process in case it is using the 125 GeV template.
    if(sm_like_hists == "smh125")
    {
      mass_var.erase("bbH1");
      process_norm_map.erase("bbH1");
    }
  }

  dout("[INFO] Prepare demo.");
  if(do_morph && (analysis == "bsm-model-indep" || analysis == "bsm-model-dep-additional" || analysis == "bsm-model-dep-full"))
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
        for (auto proc : mssm_ggH_signals)
        {
          std::string prenorm_name = bin + "_" + proc + "_morph_prenorm";
          std::string norm_name = bin + "_" + proc + "_morph_norm";
          ws.factory(TString::Format("expr::%s('@0*@1',%s, %s_frac)", norm_name.c_str(), prenorm_name.c_str(), proc.c_str()));
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
  else if(!do_morph && analysis == "bsm-model-indep"){

   //double Tfrac = ws.function("ggh_t_frac")->getVal();
   //double Bfrac = ws.function("ggh_b_frac")->getVal();
   //double Ifrac = ws.function("ggh_i_frac")->getVal();
   double Tfrac=1., Bfrac=0., Ifrac=0.; // use t-only when no morphing option is used
   //if (Ifrac<0.) {
   //  Ifrac=fabs(Ifrac);
   //  // set a constant rate parameter = -1 for the interference 
   //  cb.cp()
   //   .process({"ggh_i"})
   //   .AddSyst(cb, "rate_minus","rateParam",SystMap<>::init(-1.0));
   //  cb.GetParameter("rate_minus")->set_range(-1.0,-1.0);
   //}
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


  std::cout << "[INFO] Writing datacards to " << output_folder << std::endl;
    // We need to do this to make sure the ttbarShape uncertainty is added properly when we use a shapeU
  cb.GetParameter("CMS_htt_ttbarShape")->set_err_d(-1.);
  cb.GetParameter("CMS_htt_ttbarShape")->set_err_u(1.);


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
}
