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

void ReverseBins(TH1F & h) {
  std::vector<float> contents(h.GetNbinsX());
  std::vector<float> errors(h.GetNbinsX());
  for (int i = 0; i < h.GetNbinsX(); ++i) {
    contents[i] = h.GetBinContent(i + 1);
    errors[i] = h.GetBinError(i + 1);
  }
  for (int i = 0; i < h.GetNbinsX(); ++i) {
    h.SetBinContent(h.GetNbinsX() - i, contents[i]);
    h.SetBinError(h.GetNbinsX() - i, errors[i]);
  }
}

int main(int argc, char* argv[]) {
  // Need this to read combine workspaces
  gSystem->Load("libHiggsAnalysisCombinedLimit");

  string datacard   = "";
  string category   = "";
  string workspace  = "";
  string output     = "";
  std::string freeze_arg = "";
  string data       = "data_obs";

  po::options_description help_config("Help");
  help_config.add_options()
    ("help,h", "produce help message");

  po::options_description config("Configuration");
  config.add_options()
    ("workspace,w",
      po::value<string>(&workspace)->required(),
      "The input workspace-containing file [REQUIRED]")
    ("category,c",
      po::value<string>(&category)->required(),
      "The the cateogory for which shapes should be created [REQUIRED]")
    ("datacard,d",
      po::value<string>(&datacard),
      "The input datacard, only used for rebinning")
    ("freeze",
      po::value<string>(&freeze_arg)->default_value(freeze_arg),
      "Format PARAM1,PARAM2=X,PARAM3=Y where the values X and Y are optional")
    ("output,o ",
      po::value<string>(&output)->required(),
      "Name of the output root file to create [REQUIRED]");

  po::variables_map vm;

  // First check if the user has set the "--help" or "-h" option, and if so
  // just prin the usage information and quit
  po::store(po::command_line_parser(argc, argv)
    .options(help_config).allow_unregistered().run(), vm);
  po::notify(vm);
  // Parse the main config options
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  TFile infile(workspace.c_str());

  RooWorkspace *ws = dynamic_cast<RooWorkspace*>(gDirectory->Get("w"));

  if (!ws) {
    throw std::runtime_error(
        FNERROR("Could not locate workspace in input file"));
  }

  // Create CH instance and parse the workspace
  ch::CombineHarvester cmb;
  cmb.SetFlag("workspaces-use-clone", true);
  ch::ParseCombineWorkspace(cmb, *ws, "ModelConfig", data, false);

  // Only evaluate in case parameters to freeze are provided
  if(! freeze_arg.empty())
  {
    vector<string> freeze_vec;
    boost::split(freeze_vec, freeze_arg, boost::is_any_of(","));
    for (auto const& item : freeze_vec) {
      vector<string> parts;
      boost::split(parts, item, boost::is_any_of("="));
      if (parts.size() == 1) {
        ch::Parameter *par = cmb.GetParameter(parts[0]);
        if (par) par->set_frozen(true);
        else throw std::runtime_error(
          FNERROR("Requested variable to freeze does not exist in workspace"));
      } else {
        if (parts.size() == 2) {
          ch::Parameter *par = cmb.GetParameter(parts[0]);
          if (par) {
            par->set_val(boost::lexical_cast<double>(parts[1]));
            par->set_frozen(true);
          }
          else throw std::runtime_error(
            FNERROR("Requested variable to freeze does not exist in workspace"));
        }
      }
    }
  }

  ch::CombineHarvester cmb_card;
  cmb_card.SetFlag("workspaces-use-clone",true);
  if (datacard != "") {
    cmb_card.ParseDatacard(datacard, "", "", "", 0,"125");
  }

  // Drop any process that has no hist/data/pdf
  cmb.FilterProcs([&](ch::Process * proc) {
    bool no_shape = !proc->shape() && !proc->data() && !proc->pdf();
    if (no_shape) {
      cout << "Filtering process with no shape:\n";
      cout << ch::Process::PrintHeader << *proc << "\n";
    }
    return no_shape;
  });

  cmb = cmb.bin({category});
  cmb_card = cmb_card.bin({category});

  TFile outfile(output.c_str(), "RECREATE");
  TH1::AddDirectory(false);

  // Also create a simple map for storing total histograms, summed 
  // over all bins, in the form:
  //   pre_shapes_tot[<process>]
  map<string, TH1F> pre_shapes_tot;

  // We can always do the prefit version,
  // Loop through the bins writing the shapes to the output file
  pre_shapes_tot["data_obs"] = cmb.GetObservedShape();
  // Then fill total signal and total bkg hists
  std::cout << ">> Doing prefit: TotalBkg" << std::endl;
  pre_shapes_tot["TotalBkg"] =
      cmb.cp().backgrounds().GetShapeWithUncertainty();

  ch::CombineHarvester cmb_sm = cmb.cp();
  ch::Parameter *par_sm = cmb_sm.GetParameter("x");
  par_sm->set_val(0.0);
  std::vector<std::string> sm_signals = {"ggH125", "qqH125", "WH125", "ZH125", "ttH125"};
  for (auto proc : sm_signals) {
    cmb_sm.cp().process({proc}).PrintProcs();
    pre_shapes_tot[proc] = cmb_sm.cp().process({proc}).GetShapeWithUncertainty();
  }

  ch::CombineHarvester cmb_mssm = cmb.cp();
  ch::Parameter *par_mssm = cmb_mssm.GetParameter("x");
  par_mssm->set_val(1.0);
  std::vector<std::string> mssm_signals = {"ggh_i", "ggh_t", "ggh_b", "ggA_i", "ggA_t", "ggA_b", "ggH_i", "ggH_t", "ggH_b", "bbh", "bbH", "bbA", "qqh"};
  for (auto proc : mssm_signals) {
    cmb_mssm.cp().process({proc}).PrintProcs();
    pre_shapes_tot[proc] = cmb_mssm.cp().process({proc}).GetShapeWithUncertainty();
  }

  if (datacard != "") {
    TH1F ref = cmb_card.cp().GetObservedShape();
    for (auto & it : pre_shapes_tot) {
      it.second = ch::RestoreBinning(it.second, ref);
    }
  }

  // Can write these straight into the output file
  outfile.cd();
  for (auto& iter : pre_shapes_tot) {
    ch::WriteToTFile(&(iter.second), &outfile, category + "_prefit/" + iter.first);
  }

  // And we're done!
  outfile.Close();
  return 0;
}

