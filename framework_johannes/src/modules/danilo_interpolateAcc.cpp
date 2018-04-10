/*
 * Interpolate signal scan acceptance to fix bad statistics for GGM scan
 */
#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"
 
#include <TChain.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TLine.h>
#include <TTreeReader.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TColor.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <regex>

using namespace std;

static Config const &cfg=Config::get();

int ReadFile(std::string fname, std::map<std::string, float> *m) {
	int count = 0;
	if (access(fname.c_str(), R_OK) < 0)
		return -errno;

	FILE *fp = fopen(fname.c_str(), "r");
	if (!fp)
		return -errno;

	m->clear();

	char *buf = 0;
	size_t buflen = 0;

	while(getline(&buf, &buflen, fp) > 0) {
		char *nl = strchr(buf, '\n');
		if (nl == NULL)
			continue;
		*nl = 0;

		char *sep = strchr(buf, '=');
		if (sep == NULL)
			continue;
		*sep = 0;
		sep++;

		std::string s1 = buf;
		float s2 = std::stof(sep);

		(*m)[s1] = s2;

		count++;
	}

	if (buf)
		free(buf);

	fclose(fp);
	return count;
}


std::pair<int,int> getMasses(std::string fileName) {
	//~ std::cout<<fileName<<std::endl;
	std::smatch m;
	std::regex e;
	e="GGM_M1_M2_(.*)_(.*)";

	std::regex_search (fileName,m,e);
	assert(m.size()==3);
	return std::make_pair(std::stoi(m[1]),std::stoi(m[2]));
}

int interpolateAcc() {
	std::map<std::string,float> in;
	TGraph2D graphAcc;
	ReadFile("/home/home4/institut_1b/dmeuser/master_code/framework_johannes/output/stuff/accMap.txt",&in);
	
	for (auto const &map:in) {
		std::pair<int,int> const masses=getMasses(map.first);
		if (masses.first>800 && (masses.second==450)) continue;
		graphAcc.SetPoint(graphAcc.GetN(),masses.first,masses.second,map.second);		
	}
	
	for (auto const &map:in) {
		std::pair<int,int> const masses=getMasses(map.first);
		if (masses.first>800 && (masses.second==450)) {
			std::cout<<masses.first<<"  "<<masses.second<<"  "<<map.second<<"  "<<graphAcc.Interpolate(masses.first,masses.second)<<std::endl;
		}
	}
	
	io::RootFileSaver saver_hist("test.root","",false);
	saver_hist.save(graphAcc,"interAcc");
	
	return 0;
}

extern "C"

void run(){
	interpolateAcc();
}

