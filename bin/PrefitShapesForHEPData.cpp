#include <map>
#include "boost/program_options.hpp"
#include "boost/format.hpp"
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
  std::string category     = "";
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
    ("category,c ",
      po::value<string>(&cat)->required(),
      "The name of the category [REQUIRED]")
    ("signal-mass-name,m",
      po::value<string>(&signal_mass_name)->required(),
      "The name of the mass parameter in the workspace [REQUIRED]")
    ("signal-masses,M ",
      po::value<std::vector<unsigned int>>(&signal_masses)->required(),
      "The list of signal mass values to be considered [REQUIRED]")
    ("poi-list,P ",
      po::value<std::vector<std::string>>(&parameters_of_interest)->required(),
      "The list of parameters of interest and the values they should be given. Each element should be of the following structure: name_value [REQUIRED]")
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
    return 1;
  }

  // Parse the main config options
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  // Need this to read combine workspaces
  gSystem->Load("libHiggsAnalysisCombinedLimit");

  // Get workspace of the analysis category htt_tt_35_2018
  std::string category_name = "htt_tt_35_2018";
  auto inputfile = TFile::Open(("ws_" + category_name +  ".root").c_str(), "read");
  RooWorkspace *ws = (RooWorkspace*)inputfile->Get("w");

  // Getting datacard + input root file for restoring the binning
  ch::CombineHarvester cmb_card;
  cmb_card.SetFlag("workspaces-use-clone",true);
  cmb_card.ParseDatacard((category_name + ".txt").c_str(), "", "", "", 0, "125");
  TH1F reference_binning = cmb_card.cp().GetObservedShape();

  // Create CH instance and parse the workspace
  ch::CombineHarvester cmb;
  cmb.SetFlag("workspaces-use-clone", true);
  ch::ParseCombineWorkspace(cmb, *ws, "ModelConfig", "data_obs", false);

  // Drop any process that has no hist/data/pdf
  cmb.FilterProcs([&](ch::Process * proc) {
    bool no_shape = !proc->shape() && !proc->data() && !proc->pdf();
    if (no_shape) {
      cout << "Filtering process with no shape:\n";
      cout << ch::Process::PrintHeader << *proc << "\n";
    }
    return no_shape;
  });

  // Loop through categories
  auto bins = cmb.cp().bin_set();
  for (auto bin : bins) {
    ch::CombineHarvester cmb_bin = cmb.cp().bin({bin});

    // Storing all info in an output file
    TFile* out = TFile::Open((bin+"_hepdata.root").c_str(), "recreate");
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
          continue; // Avoid unconstrained nuisance parameters
        }
        else
        {
          // Determine impact of uncertainty variation
          //std::cout << "\tSystematic: " << syst.name() << ", value: " << syst.val() <<", -1 sigma: " << syst.err_d() << ", +1 sigma: " << syst.err_u() << std::endl;
          syst.set_val(1);
          cmb_bin_bgproc.cp().UpdateParameters({syst});
          auto bg_shape_syst_up = cmb_bin_bgproc.cp().GetShape();
          bg_shape_syst_up = ch::RestoreBinning(bg_shape_syst_up, reference_binning);
          //std::cout << "\tUpdated integral of background (upward): " << bg_shape_syst_up.Integral() << std::endl;
          std::string up_name;
          syst.set_val(-1);
          cmb_bin_bgproc.cp().UpdateParameters({syst});
          auto bg_shape_syst_down = cmb_bin_bgproc.cp().GetShape();
          bg_shape_syst_down = ch::RestoreBinning(bg_shape_syst_down, reference_binning);
          //std::cout << "\tUpdated integral of background (downward): " << bg_shape_syst_down.Integral() << std::endl;
          std::string down_name;
          syst.set_val(0);
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
    auto MH = cmb_bin_sig.GetParameter("MH");
    auto r_ggH = cmb_bin_sig.GetParameter("r_ggH");
    auto r_bbH = cmb_bin_sig.GetParameter("r_bbH");
    r_ggH->set_val(1);
    r_bbH->set_val(1);
    std::vector<unsigned int> masses = {60, 80, 95, 100, 120, 125, 130, 140, 160, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2300, 2600, 2900, 3200, 3500};
    for (auto m : masses){
      MH->set_val(m);
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
            continue; // Avoid unconstrained parameters
          }
          else
          {
            // Determine impact of uncertainty variation
            //std::cout << "\tSystematic: " << syst.name() << ", value: " << syst.val() <<", -1 sigma: " << syst.err_d() << ", +1 sigma: " << syst.err_u() << std::endl;
            syst.set_val(1);
            cmb_bin_sigproc.cp().UpdateParameters({syst});
            auto sig_shape_syst_up = cmb_bin_sigproc.cp().GetShape();
            sig_shape_syst_up = ch::RestoreBinning(sig_shape_syst_up, reference_binning);
            //std::cout << "\tUpdated integral of signal (upward): " << sig_shape_syst_up.Integral() << std::endl;
            std::string up_name;
            syst.set_val(-1);
            cmb_bin_sigproc.cp().UpdateParameters({syst});
            auto sig_shape_syst_down = cmb_bin_sigproc.cp().GetShape();
            sig_shape_syst_down = ch::RestoreBinning(sig_shape_syst_down, reference_binning);
            //std::cout << "\tUpdated integral of signal (downward): " << sig_shape_syst_down.Integral() << std::endl;
            std::string down_name;
            syst.set_val(0);
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
  }
  return 0;
}
