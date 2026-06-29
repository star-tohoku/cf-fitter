// ROOT macro wrapper — prefer: build/cf-fit --demo --mode both
// Loads shared library built by CMake (set CF_FITTER_BUILD before sourcing).
void fitCf(const char* input = nullptr,
           const char* histName = "hCF",
           const char* channel = "phi_proton",
           const char* scenario = "HAL",
           const char* mode = "both",
           const char* configRoot = nullptr) {
  TString cmd = Form("%s/cf-fit", gSystem->Getenv("CF_FITTER_BUILD") ? gSystem->Getenv("CF_FITTER_BUILD") : "build");
  if (!input || !input[0]) cmd += " --demo";
  else cmd += Form(" --input %s --hist %s", input, histName);
  cmd += Form(" --channel %s --scenario %s --mode %s", channel, scenario, mode);
  if (configRoot && configRoot[0]) cmd += Form(" --config %s", configRoot);
  gSystem->Exec(cmd.Data());
}
