import ROOT
output_filename = "1d_hid.root"
output_file = ROOT.TFile(output_filename, "RECREATE")
# Create a 1D histogram with 163 bins
hid = ROOT.TH1I("hid", "hid", 163, 0, 163)

# Loop over the bins and set the content to 1
for i in range(0, hid.GetNbinsX() + 1):
    hid.SetBinContent(i, 1)

hid.Write()
output_file.Close()


