#include <map>
#include <iostream>
#include <fstream>
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/algorithm/string.hpp"
#include "TSystem.h"
#include "TH2F.h"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/ParseCombineWorkspace.h"
#include "CombineHarvester/CombineTools/interface/TFileIO.h"
#include "CombineHarvester/CombineTools/interface/Logging.h"

namespace po = boost::program_options;

using namespace std;

int main(int argc, char* argv[]) {

  // Setup program options
  std::string datacard   = "";
  std::string workspace  = "";
  std::string fit = "";
  std::string fit_name = "";
  std::string category_name     = "";
  std::string signal_mass_name = "";
  std::vector<unsigned int> signal_masses;
  std::vector<std::string> parameters_of_interest;
  std::string data       = "data_obs";

  po::options_description help_config("Help");
  help_config.add_options()
    ("help,h", "produce help message");

  po::options_description config("Configuration");
  config.add_options()
    ("workspace,w",
      po::value<string>(&workspace)->required(),
      "The input ROOT file containing the workspace [REQUIRED]")
    ("datacard,d",
      po::value<string>(&datacard)->required(),
      "The input datacard to be used for restoring the binning. Please note, that the input files mentioned in the datacard should also be accessible [REQUIRED]")
    ("fit,f",
      po::value<string>(&fit)->required(),
      "The input ROOT file containing the background-only fit [REQUIRED]")
    ("fitname,F",
      po::value<string>(&fit_name)->required(),
      "Name the background-only fit [REQUIRED]")
    ("category,c ",
      po::value<string>(&category_name)->required(),
      "The name of the category [REQUIRED]")
    ("signal-mass-name,m",
      po::value<string>(&signal_mass_name)->required(),
      "The name of the mass parameter in the workspace [REQUIRED]")
    ("signal-masses,M ",
      po::value<std::vector<unsigned int>>(&signal_masses)->required(),
      "The list of signal mass values to be considered [REQUIRED]")
    ("poi-list,P ",
      po::value<std::vector<std::string>>(&parameters_of_interest)->required(),
      "The list of parameters of interest and the values they should be given. Each element should be of the following structure: name:value [REQUIRED]")
    ("data",
      po::value<string>(&data)->default_value(data),
      "The name of observed, measured data");

  po::variables_map vm;

  // First check if the user has set the "--help" or "-h" option, and if so
  // just prin the usage information and quit
  po::store(po::command_line_parser(argc, argv)
    .options(help_config).allow_unregistered().run(), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << config << std::endl;
    std::cout << "Example usage:" << std::endl << std::endl;
    std::cout << "PostFitShapesForHEPData -w ws_htt_tt_35_2018.root -d htt_tt_35_2018.txt -c htt_tt_35_2018 \\" << std::endl;
    std::cout << "                       -P r_ggH:1 -P r_bbH:1 -m MH -f multidimfitggH.bkgOnly.bestfit.robustHesse.root -F fit_mdf \\" << std::endl;
    std::cout << "                       -M 60 -M 80 -M 95 -M 100 -M 120 -M 125 -M 130 -M 140 -M 160 -M 180 -M 200 \\" << std::endl;
    std::cout << "                       -M 250 -M 300 -M 350 -M 400 -M 450 -M 500 -M 600 -M 700 -M 800 -M 900 -M 1000 \\" << std::endl;
    std::cout << "                       -M 1200 -M 1400 -M 1600 -M 1800 -M 2000 -M 2300 -M 2600 -M 2900 -M 3200 -M 3500" << std::endl;
    return 1;
  }

  // Parse the main config options
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  // Need this to read combine workspaces
  gSystem->Load("libHiggsAnalysisCombinedLimit");

  // Get workspace
  auto inputfile = TFile::Open(workspace.c_str(), "read");
  RooWorkspace *ws = (RooWorkspace*)inputfile->Get("w");

  // Obtain background-only fit
  auto fitfile = TFile::Open(fit.c_str(), "read");
  RooFitResult* fitres = (RooFitResult*)fitfile->Get(fit_name.c_str());
  auto postfit_parameters = fitres->floatParsFinal();
  TIterator* iter(postfit_parameters.createIterator());
  std::vector<std::string> parnames = {};

  // Getting datacard + input root file for restoring the binning
  ch::CombineHarvester cmb_card;
  cmb_card.SetFlag("workspaces-use-clone",true);
  cmb_card.ParseDatacard(datacard.c_str(), "", "", "", 0, "125");
  TH1F reference_binning = cmb_card.cp().GetObservedShape();

  // Create CH instance and parse the workspace
  ch::CombineHarvester cmb;
  cmb.SetFlag("workspaces-use-clone", true);
  ch::ParseCombineWorkspace(cmb, *ws, "ModelConfig", data.c_str(), false);

  // Apply post-fit parameter values from fit
  for (TObject *parit = iter->Next(); parit != nullptr; parit = iter->Next()) {
    RooRealVar *postfitpar = dynamic_cast<RooRealVar *>(parit);
    auto par = cmb.cp().GetParameter(postfitpar->GetName());
    parnames.push_back(postfitpar->GetName());
    if(par){
      //std::cout << "Initial parameter: " << par->name() << ", value: " << par->val() <<", -1 sigma: " << par->err_d() << ", +1 sigma: " << par->err_u() << std::endl;
      par->set_val(postfitpar->getVal());
      par->set_err_d(postfitpar->getErrorLo());
      par->set_err_u(postfitpar->getErrorHi());
      //std::cout << "\tFinal parameter: " << par->name() << ", value: " << par->val() <<", -1 sigma: " << par->err_d() << ", +1 sigma: " << par->err_u() << std::endl;
    }
    else {
      //std::cout << "WARNING: Following parameter not in workspace: " << postfitpar->GetName() << std::endl;
    }
  }

  // Writing out correlations between parameters
  std::ofstream correlations;
  correlations.open(fit_name+"_correlations.csv");
  correlations << "Parameter1,Parameter2,Correlation\n";
  for (unsigned int par1=0; par1 < parnames.size(); ++par1){
    for(unsigned int par2=0; par2 < par1; ++par2){
      double correlation_val = fitres->correlation(parnames.at(par1).c_str(), parnames.at(par2).c_str());
      if(std::abs(correlation_val) >= 1e-2){ // Include only correlations >= 1 %
        correlations << parnames.at(par1) << "," << parnames.at(par2) << "," << correlation_val  << std::endl;
      }
    }
  }

  // Drop any process that has no hist/data/pdf
  cmb.FilterProcs([&](ch::Process * proc) {
    bool no_shape = !proc->shape() && !proc->data() && !proc->pdf();
    if (no_shape) {
      std::cout << "Filtering process with no shape:\n";
      std::cout << ch::Process::PrintHeader << *proc << "\n";
    }
    return no_shape;
  });

  ch::CombineHarvester cmb_bin = cmb.cp().bin({category_name.c_str()});

  // Storing all info in an output file
  TFile* out = TFile::Open((category_name+"_hepdata.root").c_str(), "recreate");
  out->cd();

  // Extracting observed data
  out->mkdir("data_obs");
  out->cd("data_obs");
  auto data_obs = cmb_bin.cp().GetObservedShape();
  data_obs = ch::RestoreBinning(data_obs, reference_binning);
  data_obs.SetName("data_obs");
  data_obs.SetTitle("data_obs");
  data_obs.Write();
  out->cd();

  // Extracting background processes
  ch::CombineHarvester cmb_bin_bg = cmb_bin.cp().backgrounds();
  auto bgs = cmb_bin_bg.cp().process_set();
  for (auto bg : bgs)  {
    out->mkdir(bg.c_str());
    out->cd(bg.c_str());
    ch::CombineHarvester cmb_bin_bgproc = cmb_bin_bg.cp().process({bg});
    auto bg_shape = cmb_bin_bgproc.cp().GetShape();
    bg_shape = ch::RestoreBinning(bg_shape, reference_binning);
    bg_shape.SetName(bg.c_str());
    bg_shape.SetTitle(bg.c_str());
    bg_shape.Write();
    unsigned int n_bins = bg_shape.GetNbinsX();
    //std::cout << "Background: " << bg << ", nominal yield: " << bg_shape.Integral() << std::endl;
    auto systematics = cmb_bin_bgproc.cp().GetParameters();
    for (auto syst : systematics) {
      if (0 == syst.err_d() && 0 == syst.err_u())
      {
        //std::cout << "\tExcluding: " << syst.name() << std::endl;
        continue; // Avoid parameters not used in the fit
      }
      else
      {
        // Determine impact of uncertainty variation
        //std::cout << "\tSystematic (central): " << syst.name() << ", value: " << syst.val() <<", -1 sigma: " << syst.err_d() << ", +1 sigma: " << syst.err_u() << std::endl;
        double initial_syst_val = syst.val();
        syst.set_val(initial_syst_val + syst.err_u());
        //std::cout << "\tSystematic (up): " << syst.name() << ", value: " << syst.val() <<", -1 sigma: " << syst.err_d() << ", +1 sigma: " << syst.err_u() << std::endl;
        cmb_bin_bgproc.cp().UpdateParameters({syst});
        auto bg_shape_syst_up = cmb_bin_bgproc.cp().GetShape();
        bg_shape_syst_up = ch::RestoreBinning(bg_shape_syst_up, reference_binning);
        //std::cout << "\tUpdated integral of background (upward): " << bg_shape_syst_up.Integral() << std::endl;
        std::string up_name;
        syst.set_val(initial_syst_val + syst.err_d());
        //std::cout << "\tSystematic (down): " << syst.name() << ", value: " << syst.val() <<", -1 sigma: " << syst.err_d() << ", +1 sigma: " << syst.err_u() << std::endl;
        cmb_bin_bgproc.cp().UpdateParameters({syst});
        auto bg_shape_syst_down = cmb_bin_bgproc.cp().GetShape();
        bg_shape_syst_down = ch::RestoreBinning(bg_shape_syst_down, reference_binning);
        //std::cout << "\tUpdated integral of background (downward): " << bg_shape_syst_down.Integral() << std::endl;
        std::string down_name;
        syst.set_val(initial_syst_val);
        //std::cout << "\tSystematic (restored): " << syst.name() << ", value: " << syst.val() <<", -1 sigma: " << syst.err_d() << ", +1 sigma: " << syst.err_u() << std::endl;
        cmb_bin_bgproc.cp().UpdateParameters({syst});

        // Check impact of systematic uncertainty, and excluding, if no impact found
        bool valid_systematic = false;
        TH1F* bg_shape_copy_up = (TH1F*) bg_shape.Clone();
        bg_shape_copy_up->Divide(&bg_shape_syst_up);
        TH1F* bg_shape_copy_down = (TH1F*) bg_shape.Clone();
        bg_shape_copy_down->Divide(&bg_shape_syst_down);
        std::vector<float> up_vals; up_vals.reserve(n_bins);
        std::vector<float> down_vals; down_vals.reserve(n_bins);
        for(unsigned int i=1; i <= n_bins; i++){
          //std::cout << "\t\tBin: " << i << " ContentUp: " << bg_shape_copy_up->GetBinContent(i) << " ContentDown: " << bg_shape_copy_down->GetBinContent(i) << std::endl;
            up_vals.push_back(bg_shape_copy_up->GetBinContent(i));
            down_vals.push_back(bg_shape_copy_down->GetBinContent(i));
            if(valid_systematic == false && (up_vals.back() != 1 || down_vals.back() != 1)){
              valid_systematic = true;
            }
        }

        // Determine type of systematic
        /*
        for (auto val : up_vals){
          std::cout << val << " ";
        }
        std::cout << "\n";
        for (auto val : down_vals){
          std::cout << val << " ";
        }
        std::cout << "\n";
        */
        bool up_vals_equal = std::all_of(up_vals.begin(), up_vals.end(), [&] (float test){return std::abs(test - up_vals[0]) < 1e-6;});
        bool down_vals_equal = std::all_of(down_vals.begin(), down_vals.end(), [&] (float test){return std::abs(test - down_vals[0]) < 1e-6;});
        if(up_vals_equal && down_vals_equal){
          up_name = bg + "_norm_" + syst.name() + "_Up";
          down_name = bg + "_norm_" + syst.name() + "_Down";
        }
        else if (std::string(syst.name()).find("prop") == std::string::npos){
          up_name = bg + "_shape_" + syst.name() + "_Up";
          down_name = bg + "_shape_" + syst.name() + "_Down";
        }
        else {
          up_name = bg + "_" + syst.name() + "_Up";
          down_name = bg + "_" + syst.name() + "_Down";
        }

        // If a proper uncertainty, write systematic variations to output
        if (valid_systematic){
          bg_shape_syst_up.SetName(up_name.c_str());
          bg_shape_syst_up.SetTitle(up_name.c_str());
          bg_shape_syst_up.Write();
          bg_shape_syst_down.SetName(down_name.c_str());
          bg_shape_syst_down.SetTitle(down_name.c_str());
          bg_shape_syst_down.Write();
        }
        else {
          //std::cout << "\tExcluding: " << syst.name() << std::endl;
        }
        //break;
      }
    }
    //break;
    out->cd();
  }

  // Extracting signal processes
  ch::CombineHarvester cmb_bin_sig = cmb_bin.cp().signals();
  auto mass = cmb_bin_sig.GetParameter(signal_mass_name);
  std::map <std::string, ch::Parameter*> pois;
  for (auto par: parameters_of_interest){
    std::vector<std::string> par_infos;
    boost::split(par_infos, par, [](char c){return c == ':';});
    pois[par_infos.at(0)] = cmb_bin_sig.GetParameter(par_infos.at(0));
    pois[par_infos.at(0)]->set_val(std::stof(par_infos.at(1)));
  }
  for (auto m : signal_masses){
    mass->set_val(m);
    auto sigs = cmb_bin_sig.cp().process_set();

    for (auto sig : sigs){
      ch::CombineHarvester cmb_bin_sigproc = cmb_bin_sig.cp().process({sig});
      auto sig_shape = cmb_bin_sigproc.cp().GetShape();
      sig_shape = ch::RestoreBinning(sig_shape, reference_binning);
      std::string sig_name = sig + "_" + std::to_string(m);
      out->mkdir(sig_name.c_str());
      out->cd(sig_name.c_str());
      sig_shape.SetName(sig_name.c_str());
      sig_shape.SetTitle(sig_name.c_str());
      sig_shape.Write();
      unsigned int n_bins = sig_shape.GetNbinsX();
      //std::cout << "Signal: " << sig_name << ", nominal yield: " << sig_shape.Integral() << std::endl;

      auto systematics = cmb_bin_sigproc.cp().GetParameters();
      for (auto syst : systematics){
        if (0 == syst.err_d() && 0 == syst.err_u())
        {
          //std::cout << "\tExcluding from variations: " << syst.name() << ", value: " << syst.val() <<", -1 sigma: " << syst.err_d() << ", +1 sigma: " << syst.err_u() << std::endl;
          continue; // Avoid parameters not used in the fit
        }
        else
        {
          // Determine impact of uncertainty variation
          //std::cout << "\tSystematic (central): " << syst.name() << ", value: " << syst.val() <<", -1 sigma: " << syst.err_d() << ", +1 sigma: " << syst.err_u() << std::endl;
          double initial_syst_val = syst.val();
          syst.set_val(initial_syst_val + syst.err_u());
          //std::cout << "\tSystematic (up): " << syst.name() << ", value: " << syst.val() <<", -1 sigma: " << syst.err_d() << ", +1 sigma: " << syst.err_u() << std::endl;
          cmb_bin_sigproc.cp().UpdateParameters({syst});
          auto sig_shape_syst_up = cmb_bin_sigproc.cp().GetShape();
          sig_shape_syst_up = ch::RestoreBinning(sig_shape_syst_up, reference_binning);
          //std::cout << "\tUpdated integral of signal (upward): " << sig_shape_syst_up.Integral() << std::endl;
          std::string up_name;
          syst.set_val(initial_syst_val + syst.err_d());
          //std::cout << "\tSystematic (down): " << syst.name() << ", value: " << syst.val() <<", -1 sigma: " << syst.err_d() << ", +1 sigma: " << syst.err_u() << std::endl;
          cmb_bin_sigproc.cp().UpdateParameters({syst});
          auto sig_shape_syst_down = cmb_bin_sigproc.cp().GetShape();
          sig_shape_syst_down = ch::RestoreBinning(sig_shape_syst_down, reference_binning);
          //std::cout << "\tUpdated integral of signal (downward): " << sig_shape_syst_down.Integral() << std::endl;
          std::string down_name;
          syst.set_val(initial_syst_val);
          //std::cout << "\tSystematic (restored): " << syst.name() << ", value: " << syst.val() <<", -1 sigma: " << syst.err_d() << ", +1 sigma: " << syst.err_u() << std::endl;
          cmb_bin_sigproc.cp().UpdateParameters({syst});

          // Check impact of systematic uncertainty, and excluding, if no impact found
          bool valid_systematic = false;
          TH1F* sig_shape_copy_up = (TH1F*) sig_shape.Clone();
          sig_shape_copy_up->Divide(&sig_shape_syst_up);
          TH1F* sig_shape_copy_down = (TH1F*) sig_shape.Clone();
          sig_shape_copy_down->Divide(&sig_shape_syst_down);
          std::vector<float> up_vals; up_vals.reserve(n_bins);
          std::vector<float> down_vals; down_vals.reserve(n_bins);
          for(unsigned int i=1; i <= n_bins; i++){
            //std::cout << "\t\tBin: " << i << " ContentUp: " << sig_shape_copy_up->GetBinContent(i) << " ContentDown: " << sig_shape_copy_down->GetBinContent(i) << std::endl;
            up_vals.push_back(sig_shape_copy_up->GetBinContent(i));
            down_vals.push_back(sig_shape_copy_down->GetBinContent(i));
            if(valid_systematic == false && (up_vals.back() != 1 || down_vals.back() != 1)){
              valid_systematic = true;
            }
          }

          // Determine type of systematic
          /*
          for (auto val : up_vals){
            std::cout << val << " ";
          }
          std::cout << "\n";
          for (auto val : down_vals){
            std::cout << val << " ";
          }
          std::cout << "\n";
          */
          bool up_vals_equal = std::all_of(up_vals.begin(), up_vals.end(), [&] (float test){return std::abs(test - up_vals[0]) < 1e-6;});
          bool down_vals_equal = std::all_of(down_vals.begin(), down_vals.end(), [&] (float test){return std::abs(test - down_vals[0]) < 1e-6;});
          if(up_vals_equal && down_vals_equal){
            up_name = sig_name + "_norm_" + syst.name() + "_Up";
            down_name = sig_name + "_norm_" + syst.name() + "_Down";
          }
          else if (std::string(syst.name()).find("prop") == std::string::npos){
            up_name = sig_name + "_shape_" + syst.name() + "_Up";
            down_name = sig_name + "_shape_" + syst.name() + "_Down";
          }
          else {
            up_name = sig_name + "_" + syst.name() + "_Up";
            down_name = sig_name + "_" + syst.name() + "_Down";
          }

          // If a proper uncertainty, write systematic variations to output
          if (valid_systematic){
            sig_shape_syst_up.SetName(up_name.c_str());
            sig_shape_syst_up.SetTitle(up_name.c_str());
            sig_shape_syst_up.Write();
            sig_shape_syst_down.SetName(down_name.c_str());
            sig_shape_syst_down.SetTitle(down_name.c_str());
            sig_shape_syst_down.Write();
          }
          else {
            //std::cout << "\tExcluding: " << syst.name() << std::endl;
          }
          //break;
        }
      }
      //break;
      out->cd();
    }
    //break;
  }
  return 0;
}
