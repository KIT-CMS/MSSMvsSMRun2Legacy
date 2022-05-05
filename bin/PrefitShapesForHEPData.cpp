#include <map>
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "TSystem.h"
#include "TH2F.h"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/ParseCombineWorkspace.h"
#include "CombineHarvester/CombineTools/interface/TFileIO.h"
#include "CombineHarvester/CombineTools/interface/Logging.h"

using namespace std;

int main(int argc, char* argv[]) {
  // Need this to read combine workspaces
  gSystem->Load("libHiggsAnalysisCombinedLimit");

  auto inputfile = TFile::Open("ws_htt_tt_35_2018.root", "read");
  RooWorkspace *ws = (RooWorkspace*)inputfile->Get("w");

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
      bg_shape.SetName(bg.c_str());
      bg_shape.SetTitle(bg.c_str());
      bg_shape.Write();
      unsigned int n_bins = bg_shape.GetNbinsX();
      std::cout << "Background: " << bg << ", nominal yield: " << bg_shape.Integral() << std::endl;
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
          //std::cout << "\tUpdated integral of background (upward): " << bg_shape_syst_up.Integral() << std::endl;
          std::string up_name = bg + "_" + syst.name() + "_Up";
          bg_shape_syst_up.SetName(up_name.c_str());
          bg_shape_syst_up.SetTitle(up_name.c_str());
          syst.set_val(-1);
          cmb_bin_bgproc.cp().UpdateParameters({syst});
          auto bg_shape_syst_down = cmb_bin_bgproc.cp().GetShape();
          //std::cout << "\tUpdated integral of background (downward): " << bg_shape_syst_down.Integral() << std::endl;
          std::string down_name = bg + "_" + syst.name() + "_Down";
          bg_shape_syst_down.SetName(down_name.c_str());
          bg_shape_syst_down.SetTitle(down_name.c_str());
          syst.set_val(0);
          cmb_bin_bgproc.cp().UpdateParameters({syst});

          // Check impact of systematic uncertainty, and excluding, if no impact found
          bool valid_systematic = false;
          TH1F* bg_shape_copy_up = (TH1F*) bg_shape.Clone();
          bg_shape_copy_up->Divide(&bg_shape_syst_up);
          TH1F* bg_shape_copy_down = (TH1F*) bg_shape.Clone();
          bg_shape_copy_down->Divide(&bg_shape_syst_down);
          for(unsigned int i=1; i <= n_bins; i++){
            //std::cout << "\t\tBin: " << i << " ContentUp: " << bg_shape_copy_up->GetBinContent(i) << " ContentDown: " << bg_shape_copy_down->GetBinContent(i) << std::endl;
            if(bg_shape_copy_up->GetBinContent(i) != 1 || bg_shape_copy_down->GetBinContent(i) != 1){
              valid_systematic = true;
              break;
            }
          }

          // If a proper uncertainty, write systematic variations to output
          if (valid_systematic){
            bg_shape_syst_up.Write();
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
        std::string sig_name = sig + "_" + std::to_string(m);
        out->mkdir(sig_name.c_str());
        out->cd(sig_name.c_str());
        sig_shape.SetName(sig_name.c_str());
        sig_shape.SetTitle(sig_name.c_str());
        sig_shape.Write();
        unsigned int n_bins = sig_shape.GetNbinsX();
        std::cout << "Signal: " << sig_name << ", nominal yield: " << sig_shape.Integral() << std::endl;

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
            //std::cout << "\tUpdated integral of signal (upward): " << sig_shape_syst_up.Integral() << std::endl;
            std::string up_name = sig_name + "_" + syst.name() + "_Up";
            sig_shape_syst_up.SetName(up_name.c_str());
            sig_shape_syst_up.SetTitle(up_name.c_str());
            syst.set_val(-1);
            cmb_bin_sigproc.cp().UpdateParameters({syst});
            auto sig_shape_syst_down = cmb_bin_sigproc.cp().GetShape();
            //std::cout << "\tUpdated integral of signal (downward): " << sig_shape_syst_down.Integral() << std::endl;
            std::string down_name = sig_name + "_" + syst.name() + "_Down";
            sig_shape_syst_down.SetName(down_name.c_str());
            sig_shape_syst_down.SetTitle(down_name.c_str());
            syst.set_val(0);
            cmb_bin_sigproc.cp().UpdateParameters({syst});

            // Check impact of systematic uncertainty, and excluding, if no impact found
            bool valid_systematic = false;
            TH1F* sig_shape_copy_up = (TH1F*) sig_shape.Clone();
            sig_shape_copy_up->Divide(&sig_shape_syst_up);
            TH1F* sig_shape_copy_down = (TH1F*) sig_shape.Clone();
            sig_shape_copy_down->Divide(&sig_shape_syst_down);
            for(unsigned int i=1; i <= n_bins; i++){
              //std::cout << "\t\tBin: " << i << " ContentUp: " << sig_shape_copy_up->GetBinContent(i) << " ContentDown: " << sig_shape_copy_down->GetBinContent(i) << std::endl;
              if(sig_shape_copy_up->GetBinContent(i) != 1 || sig_shape_copy_down->GetBinContent(i) != 1){
                valid_systematic = true;
                break;
              }
            }

            // If a proper uncertainty, write systematic variations to output
            if (valid_systematic){
              sig_shape_syst_up.Write();
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
