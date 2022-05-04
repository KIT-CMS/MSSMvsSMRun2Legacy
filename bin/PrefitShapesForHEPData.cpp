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
    auto data_obs = cmb_bin.cp().GetObservedShape();
    data_obs.SetName("data_obs");
    data_obs.SetTitle("data_obs");
    data_obs.Write();

    // Extracting background processes
    ch::CombineHarvester cmb_bin_bg = cmb_bin.cp().backgrounds();
    auto bgs = cmb_bin_bg.cp().process_set();
    for (auto bg : bgs)  {
      ch::CombineHarvester cmb_bin_bgproc = cmb_bin_bg.cp().process({bg});
      auto bg_shape = cmb_bin_bgproc.cp().GetShape();
      bg_shape.SetName(bg.c_str());
      bg_shape.SetTitle(bg.c_str());
      bg_shape.Write();
      std::cout << "Background: " << bg << ", nominal yield: " << bg_shape.Integral() << std::endl;
      auto systematics = cmb_bin_bgproc.cp().GetParameters();
      for (auto syst : systematics) {
        if (0 == syst.err_d() && 0 == syst.err_u())
        {
          std::cout << "\tExcluding: " << syst.name() << std::endl;
          continue; // Avoid unconstrained rate parameters
        }
        else
        {
          std::cout << "\tSystematic: " << syst.name() << ", value: " << syst.val() <<", -1 sigma: " << syst.err_d() << ", +1 sigma: " << syst.err_u() << std::endl;
          syst.set_val(1);
          cmb_bin_bgproc.cp().UpdateParameters({syst});
          auto bg_shape_syst_up = cmb_bin_bgproc.cp().GetShape();
          std::cout << "\tUpdated integral of background (upward): " << bg_shape_syst_up.Integral() << std::endl;
          std::string up_name = bg + "_" + syst.name() + "_Up";
          bg_shape_syst_up.SetName(up_name.c_str());
          bg_shape_syst_up.SetTitle(up_name.c_str());
          bg_shape_syst_up.Write();
          syst.set_val(-1);
          cmb_bin_bgproc.cp().UpdateParameters({syst});
          auto bg_shape_syst_down = cmb_bin_bgproc.cp().GetShape();
          std::cout << "\tUpdated integral of background (downward): " << bg_shape_syst_down.Integral() << std::endl;
          std::string down_name = bg + "_" + syst.name() + "_Down";
          bg_shape_syst_down.SetName(down_name.c_str());
          bg_shape_syst_down.SetTitle(down_name.c_str());
          bg_shape_syst_down.Write();
          syst.set_val(0);
          cmb_bin_bgproc.cp().UpdateParameters({syst});
          //break;
        }
      }
      //break;
    }
    out->Close();
  }


  return 0;
}
