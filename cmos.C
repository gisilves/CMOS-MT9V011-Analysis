#include "Riostream.h"
#include <sstream>
#include <fstream>
#include <string>

void cmos(std::string filename)
{

  TCanvas *c1 = new TCanvas();
  TString dir = gROOT->GetTutorialDir();
  dir.Append("/tree/");
  dir.ReplaceAll("/./", "/");

  std::string line;
  std::ifstream infile(filename.c_str());

  Int_t npixels = 0;
  auto f = TFile::Open("basic.root", "RECREATE");
  TH2F *h2 = new TH2F("h2", "xy distribution", 640, -0.5, 639.5, 480, -0.5, 479.5);

  int value;
  int row = 0;
  int column = 0;

  if (std::getline(infile, line))
  {
    std::istringstream iss(line);
    while (iss >> value)
    {
      if(value > 400)
      {
      cout << "Col:" << column << " Row:" << row << " Val:" << value << endl;
      }
      h2->SetBinContent(column, row, value);
      npixels++;
      column++;
      if (column == 641)
      {
        column = 0;
        row++;
      }
    }
  }

  printf(" found %d points\n", npixels);
  h2->Draw("colz");
  infile.close();

  f->Write();
}
