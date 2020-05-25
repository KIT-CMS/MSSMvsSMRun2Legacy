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

int main(int argc, char **argv) {
  typedef vector<string> VString;
  typedef vector<pair<int, string>> Categories;
  using ch::syst::bin_id;
  using ch::JoinStr;

  // Define program options
  string output_folder = "output";
  string base_path = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/MSSMvsSMRun2Legacy/shapes/";
  string sm_gg_fractions = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v3.root";
  string chan = "all";
  string category = "all";
  string heavy_mass = "all";
  string light_mass = "all";
  string variable = "m_ttvisbb";
  bool auto_rebin = false;
  bool real_data = false;
  bool binomial_bbb = false;
  bool verbose = false;
  string analysis = "nmssm"; // "sm", "sm_standard"  "mssm", "mssm_btag", "mssm_vs_sm_standard", "mssm_vs_sm_qqh", "gof"
  int era = 2016; // 2016, 2017 or 2018
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()
      ("base_path", po::value<string>(&base_path)->default_value(base_path))
      ("sm_gg_fractions", po::value<string>(&sm_gg_fractions)->default_value(sm_gg_fractions))
      ("variable", po::value<string>(&variable)->default_value(variable))
      ("channel", po::value<string>(&chan)->default_value(chan))
      ("category", po::value<string>(&category)->default_value(category))
      ("heavy_mass", po::value<string>(&heavy_mass)->default_value(heavy_mass))
      ("light_mass", po::value<string>(&light_mass)->default_value(light_mass))
      ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(auto_rebin))
      ("real_data", po::value<bool>(&real_data)->default_value(real_data))
      ("verbose", po::value<bool>(&verbose)->default_value(verbose))
      ("output_folder", po::value<string>(&output_folder)->default_value(output_folder))
      ("analysis", po::value<string>(&analysis)->default_value(analysis))
      ("binomial_bbb", po::value<bool>(&binomial_bbb)->default_value(binomial_bbb))
      ("era", po::value<int>(&era)->default_value(era));
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  // Define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
  output_folder = output_folder + "_" + analysis + "_" + variable + "_" + heavy_mass + "_" + light_mass;
  std::map<string, string> input_dir;
  input_dir["mt"] = base_path;
  input_dir["et"] = base_path;
  input_dir["tt"] = base_path;
  input_dir["em"] = base_path;

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
  if (chan == "all")
    chns = {"mt", "et", "tt", "em"};

  // Define restriction to the channel defined by '--category' option
  if(category != "all"){
    std::vector<std::string> category_split;
    boost::split(category_split, category, boost::is_any_of("_"));
    chns = {category_split.at(0)};
  }

  // Define background and signal processes
  map<string, VString> bkg_procs;
  VString bkgs, bkgs_em, sm_signals, main_sm_signals, mssm_ggH_signals, mssm_bbH_signals, mssm_signals;

  sm_signals = {"WH125", "ZH125", "ttH125"};
  main_sm_signals = {"ggH125", "qqH125"};

  mssm_signals = ch::JoinStr({mssm_ggH_signals, mssm_bbH_signals});

  bkgs = {"EMB", "ZL", "TTL", "VVL", "jetFakes", "ggHWW125", "qqHWW125"};
  bkgs_em = {"EMB", "W", "QCD", "ZL", "TTL", "VVL", "ggHWW125", "qqHWW125"};

  std::cout << "[INFO] Considering the following processes as main backgrounds:\n";

  if (chan.find("em") != std::string::npos || chan.find("all") != std::string::npos) {
    std::cout << "For em channel : \n";
    for (unsigned int i=0; i < bkgs_em.size(); i++) std::cout << bkgs_em[i] << std::endl;
  }
  if (chan.find("mt") != std::string::npos || chan.find("et") != std::string::npos || chan.find("tt") != std::string::npos || chan.find("all") != std::string::npos) {
    std::cout << "For et,mt,tt channels : \n";
    for (unsigned int i=0; i < bkgs.size(); i++) std::cout << bkgs[i] << std::endl;
  }
  bkg_procs["et"] = bkgs;
  bkg_procs["mt"] = bkgs;
  bkg_procs["tt"] = bkgs;
  bkg_procs["em"] = bkgs_em;

if(analysis == "nmssm"){
    for(auto chn : chns){
        bkg_procs[chn] = JoinStr({bkg_procs[chn],sm_signals,main_sm_signals});
    }
  }

  /*
  // vector<string> heavy_masses = {"240", "280", "320", "360", "400", "450", "500", "550", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "1800", "2000", "2500", "3000"};
  vector<string> light_mass_coarse = {"60"};//, "70", "80", "90", "100", "120", "150", "170", "190", "250", "300", "350", "400", "450", "500", "550", "600", "650", "700", "800", "900", "1000", "1100", "1200", "1300", "1400", "1600", "1800", "2000", "2200", "2400", "2600", "2800"};
  vector<string> light_mass_fine = {"60"};//, "70", "75", "80", "85", "90", "95", "100", "110", "120", "130", "150", "170", "190", "250", "300", "350", "400", "450", "500", "550", "600", "650", "700", "750", "800", "850"};

  vector<string> masses;
  vector<string> light_masses_list;

  int hm;
  int lm;
  // for(int i=0;i<heavy_masses.size();i++) {
  //   sig_procs.push_back({});
  //   light_masses_list.push_back({});
  //   hm = std::stoi(heavy_masses[i]);
  //   if(hm>1001) {
  //     for(int j=0;j<light_mass_coarse.size();j++) {
  //       lm = std::stoi(light_mass_coarse[i]);
  //       if((lm+125)>hm) {
  //         continue;
  //       }
  //         sig_procs[i].push_back("NMSSM_"+heavy_masses[i]+"_125_"+light_mass_coarse[i]);
  //         light_masses_list.push_back(lm);
  //     }
  //   }
  //   else {
  //     for(int j=0;j<light_mass_fine.size();j++) {
  //       lm = std::stoi(light_mass_fine[i]);
  //       if((lm+125)>hm) {
  //         continue;
  //       }
  //         sig_procs[i].push_back("NMSSM_"+heavy_masses[i]+"_125_"+light_mass_fine[i]);
  //         light_masses_list.push_back(lm);
  //     }
  //   }
    hm = std::stoi(heavy_mass);
    if(hm>1001) {
      for(uint i=0;i<light_mass_coarse.size();i++) {
        lm = std::stoi(light_mass_coarse[i]);
        if((lm+125)>hm) {
          continue;
        }
          light_masses_list.push_back(light_mass_coarse[i]);
      }
    }
    else {
      for(uint i=0;i<light_mass_fine.size();i++) {
        lm = std::stoi(light_mass_fine[i]);
        if((lm+125)>hm) {
          continue;
        }
          light_masses_list.push_back(light_mass_fine[i]);
      }
    }
  */
  // Specify signal processes and masses
  vector<string> sig_procs;
  // STXS stage 0: ggH and VBF processes
  sig_procs = {
      "NMSSM_"+heavy_mass+"_125_"+light_mass
  };
  // Define NMSSM model-dependent mass parameters mH, mhprime, mh
  RooRealVar mH("mH", "mH", 125., 240., 3000.);
  RooRealVar mhprime("mhprime", "mhprime", 170., 60., 2800.);
  RooRealVar mh("mh", "mh", 125., 124.9, 125.1);
  mh.setConstant(true);
  mH.setConstant(true);
  mhprime.setConstant(true);


  // Define categories
  map<string, Categories> cats;
  // STXS stage 0 categories (optimized on ggH and VBF)
  if(analysis == "nmssm"){
    cats["et"] = {
      { 1 , "et_emb"},
      { 2 , "et_tt"},
      { 3 , "et_misc"},
      { 4 , "et_ff"},
      { 5 , "et_NMSSM_MH1001to1999_boosted"},
      { 6 , "et_NMSSM_MH1001to1999_unboosted"},
      { 7 , "et_NMSSM_MH2000toinfty_boosted"},
      { 8 , "et_NMSSM_MH2000toinfty_unboosted"},
      { 9 , "et_NMSSM_MH240to320_unboosted"},
      { 10, "et_NMSSM_MH321to500_boosted"},
      { 11, "et_NMSSM_MH321to500_unboosted"},
      { 12, "et_NMSSM_MH501to700_boosted"},
      { 13, "et_NMSSM_MH501to700_unboosted"},
      { 14, "et_NMSSM_MH701to1000_boosted"},
      { 15, "et_NMSSM_MH701to1000_unboosted"}
    };



    cats["mt"] = {
      { 1 , "mt_emb"},
      { 2 , "mt_tt"},
      { 3 , "mt_misc"},
      { 4 , "mt_ff"},
      { 5 , "mt_NMSSM_MH1001to1999_boosted"},
      { 6 , "mt_NMSSM_MH1001to1999_unboosted"},
      { 7 , "mt_NMSSM_MH2000toinfty_boosted"},
      { 8 , "mt_NMSSM_MH2000toinfty_unboosted"},
      { 9 , "mt_NMSSM_MH240to320_unboosted"},
      { 10, "mt_NMSSM_MH321to500_boosted"},
      { 11, "mt_NMSSM_MH321to500_unboosted"},
      { 12, "mt_NMSSM_MH501to700_boosted"},
      { 13, "mt_NMSSM_MH501to700_unboosted"},
      { 14, "mt_NMSSM_MH701to1000_boosted"},
      { 15, "mt_NMSSM_MH701to1000_unboosted"}
    };
    cats["tt"] = {
      { 1 , "tt_emb"},
      { 2 , "tt_tt"},
      { 3 , "tt_misc"},
      { 4 , "tt_ff"},
      { 5 , "tt_NMSSM_MH1001to1999_boosted"},
      { 6 , "tt_NMSSM_MH1001to1999_unboosted"},
      { 7 , "tt_NMSSM_MH2000toinfty_boosted"},
      { 8 , "tt_NMSSM_MH2000toinfty_unboosted"},
      { 9 , "tt_NMSSM_MH240to320_unboosted"},
      { 10, "tt_NMSSM_MH321to500_boosted"},
      { 11, "tt_NMSSM_MH321to500_unboosted"},
      { 12, "tt_NMSSM_MH501to700_boosted"},
      { 13, "tt_NMSSM_MH501to700_unboosted"},
      { 14, "tt_NMSSM_MH701to1000_boosted"},
      { 15, "tt_NMSSM_MH701to1000_unboosted"}
    };
    cats["em"] = {
      { 1 , "em_emb"},
      { 2 , "em_tt"},
      { 3 , "em_misc"},
      { 4 , "em_ss"},
      { 5 , "em_NMSSM_MH1001to1999_boosted"},
      { 6 , "em_NMSSM_MH1001to1999_unboosted"},
      { 7 , "em_NMSSM_MH2000toinfty_boosted"},
      { 8 , "em_NMSSM_MH2000toinfty_unboosted"},
      { 9 , "em_NMSSM_MH240to320_unboosted"},
      { 10, "em_NMSSM_MH321to500_boosted"},
      { 11, "em_NMSSM_MH321to500_unboosted"},
      { 12, "em_NMSSM_MH501to700_boosted"},
      { 13, "em_NMSSM_MH501to700_unboosted"},
      { 14, "em_NMSSM_MH701to1000_boosted"},
      { 15, "em_NMSSM_MH701to1000_unboosted"}
    };
  }
  else throw std::runtime_error("Given categorization is not known.");

  // Create combine harverster object
  ch::CombineHarvester cb;
  cb.SetFlag("workspaces-use-clone", true);

  // Add observations and processes
  std::string era_tag;
  if (era == 2016) era_tag = "2016";
  else if (era == 2017) era_tag = "2017";
  else if (era == 2018) era_tag = "2018";

  else std::runtime_error("Given era is not implemented.");

  std::vector<int> mttbb_categories = {1,2}; 

  for (auto chn : chns) {
    cb.AddObservations({"*"}, {"htt"}, {era_tag}, {chn}, cats[chn]);
    cb.AddProcesses({"*"}, {"htt"}, {era_tag}, {chn}, bkg_procs[chn], cats[chn],
                    false);
    cb.AddProcesses({"*"}, {"htt"}, {era_tag}, {chn}, sig_procs, cats[chn],
                    true);
  }

  // Add systematics
  ch::AddMSSMvsSMRun2Systematics(cb, true, true, true, true, true, era, heavy_mass, light_mass);

  // Define restriction to the desired category
  // if(category != "all"){
  //   cb = cb.bin({category});
  // }

  // Extract shapes from input ROOT files
  for (string chn : chns) {
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
          input_dir[chn] + "htt_" + chn + ".inputs-nmssm-" + era_tag + "-" + variable + ".root", "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");
    if(analysis == "nmssm"){
      cb.cp().channel({chn}).process(sig_procs).ExtractShapes(
          input_dir[chn] + "htt_" + chn + ".inputs-nmssm-" + era_tag + "-" + variable + ".root", "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");
      }
    }

  // Delete processes with 0 yield
  cb.FilterProcs([&](ch::Process *p) {
    if (std::find(sig_procs.begin(), sig_procs.end(), p->process()) != sig_procs.end())
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

  // Delete systematics with 0 yield since these result in a bogus norm error in combine
  cb.FilterSysts([&](ch::Systematic *s) {
    if (std::find(sig_procs.begin(), sig_procs.end(), s->process()) != sig_procs.end())
    {
      return false;
    }
    if (s->type() == "shape") {
      if (s->shape_u()->Integral() == 0.0) {
        std::cout << "[WARNING] Removing systematic with null yield in up shift:" << std::endl;
        std::cout << ch::Systematic::PrintHeader << *s << "\n";
        return true;
      }
      if (s->shape_d()->Integral() == 0.0) {
        std::cout << "[WARNING] Removing systematic with null yield in down shift:" << std::endl;
        std::cout << ch::Systematic::PrintHeader << *s << "\n";
        return true;
      }
    }
    return false;
  });

  // Replacing observation with the sum of the backgrounds (Asimov data)
  // useful to be able to check this, so don't do the replacement
  // for these
  if (!real_data) {
    for (auto b : cb.cp().bin_set()) {
      std::cout << "[INFO] Replacing data with asimov in bin " << b << "\n";
      auto background_shape = cb.cp().bin({b}).backgrounds().GetShape();
      auto total_procs_shape = cb.cp().bin({b}).data().GetShape();
      total_procs_shape.Scale(0.0);
      if(analysis == "nmssm"){
        auto sm_signal_shape = cb.cp().bin({b}).process(ch::JoinStr({sm_signals, main_sm_signals})).GetShape();
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
          total_procs_shape = total_procs_shape + background_shape + sm_signal_shape;
        }
        cb.cp().bin({b}).ForEachObs([&](ch::Observation *obs) {
          obs->set_shape(total_procs_shape,true);
        });
      }
    }
  }

  // At this point we can fix the negative bins
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

  // Perform auto-rebinning
  if (auto_rebin) {
    std::cout << "[INFO] Performing auto-rebinning.\n";
    auto rebin = ch::AutoRebin().SetBinThreshold(10.0).SetBinUncertFraction(0.9).SetRebinMode(1).SetPerformRebin(true).SetVerbosity(1);
    rebin.Rebin(cb, cb);
  }

  // Merge bins and set bin-by-bin uncertainties if no autoMCStats is used.
  if (binomial_bbb) {
    std::cout << "[INFO] Adding binomial bbb.\n";
    auto bbb = ch::BinomialBinByBinFactory()
                   .SetPattern("CMS_$ANALYSIS_$CHANNEL_$BIN_$ERA_$PROCESS_binomial_bin_$#")
                   .SetBinomialP(0.022)
                   .SetBinomialN(1000.0)
                   .SetFixNorm(false);
    bbb.AddBinomialBinByBin(cb.cp().channel({"em"}).process({"EMB"}), cb);
  }

  // This function modifies every entry to have a standardised bin name of
  // the form: {analysis}_{channel}_{bin_id}_{era}
  ch::SetStandardBinNames(cb, "$ANALYSIS_$CHANNEL_$BINID_$ERA");
  ch::CombineHarvester cb_obs = cb.deep().backgrounds();

  // Adding bin-by-bin uncertainties
  std::cout << "[INFO] Adding bin-by-bin uncertainties.\n";
  cb.SetAutoMCStats(cb, 0.0);
  
  // Setup morphed mssm signals for model-independent case
  RooWorkspace ws("htt", "htt");

  std::map<std::string, RooAbsReal *> mass_var = {};

  std::map<std::string, std::string> process_norm_map = {};

  /*
  
  if(analysis == "nmssm")
  {
    mass_var = {
      {"NMSSM_"+heavy_mass+"_125", &mhprime}
    };

    process_norm_map = {
      {"NMSSM_"+heavy_mass+"_125", "prenorm"}
    };
  }
  
  if(analysis == "nmssm" && false)
  {
    TFile morphing_demo(("htt_mssm_morphing_" + category+ "_demo.root").c_str(), "RECREATE");

    // Perform morphing
    // auto sig_procs = ch::JoinStr({sig_procs});
    std::cout << "[INFO] Performing template morphing for nmssm.\n";
    auto morphFactory = ch::CMSHistFuncFactory();
    morphFactory.SetHorizontalMorphingVariable(mass_var);
    morphFactory.Run(cb, ws, process_norm_map);
  

    // if(analysis == "nmssm"){
    //   // Adding 'norm' terms into workspace according to desired signals
    //   for (auto bin : cb.cp().bin_set())
    //   {
    //     for (auto proc : mssm_ggH_signals)
    //     {
    //       std::string prenorm_name = bin + "_" + proc + "_morph_prenorm";
    //       std::string norm_name = bin + "_" + proc + "_morph_norm";
    //       ws.factory(TString::Format("expr::%s('@0*@1',%s, %s_frac)", norm_name.c_str(), prenorm_name.c_str(), proc.c_str()));
    //     }
    //   }
    // }

    // Saving workspace with morphed signals
    morphing_demo.cd();
    if (verbose)
      ws.Print();
    ws.Write();
    morphing_demo.Close();
    cb.AddWorkspace(ws);
    cb.ExtractData("htt", "$BIN_data_obs");
    cb.ExtractPdfs(cb, "htt", "$BIN_$PROCESS_morph");
    std::cout << "[INFO] Finished template morphing.\n";
  }
  */
  std::cout << "[INFO] Writing datacards to " << output_folder << std::endl;

  // Decide, how to write out the datacards depending on --category option
  // Write out datacards. Naming convention important for rest of workflow. We
  // make one directory per chn-cat, one per chn and cmb. In this code we only
  // store the individual datacards for each directory to be combined later.
  ch::CardWriter writer(output_folder + "/" + era_tag + "/$TAG/$BIN.txt",
                        output_folder + "/" + era_tag + "/$TAG/common/htt_input_" + era_tag + ".root");
  std::cout << "after writer" << std::endl;
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

  if (verbose)
    cb.PrintAll();

  std::cout << "[INFO] Done producing datacards.\n";
}
