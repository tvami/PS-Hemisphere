#include <algorithm>
#include <iostream>
#include <cstring>
#include <vector>
#include <complex>
#include <cmath>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TLorentzVector.h"

#ifdef __CLING__
R__LOAD_LIBRARY(/Users/blackmac/Software/Delphes/libDelphes.so)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#endif

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Simple helper functions with inline keyword for ROOT JIT
inline double c_sum(const std::vector<double>& vec) {
    double sum = 0.0;
    for (size_t i = 0; i < vec.size(); ++i) {
        sum += vec[i];
    }
    return sum;
}

inline std::vector<double> calculate_power(const std::vector<double>& vec, double exponent) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = std::pow(vec[i], exponent);
    }
    return result;
}

inline std::vector<std::vector<double>> order_by_pt(const std::vector<std::vector<double>>& pN) {
    std::vector<std::vector<double>> result = pN;
    
    std::sort(result.begin(), result.end(), 
        [](const std::vector<double>& a, const std::vector<double>& b) {
            double pt_a = std::sqrt(a[1]*a[1] + a[2]*a[2]);
            double pt_b = std::sqrt(b[1]*b[1] + b[2]*b[2]);
            return pt_a > pt_b;
        });
    
    return result;
}

// Main function
inline std::vector<double> psDistance(
    std::vector<std::vector<double>> pN_evt1,
    std::vector<std::vector<double>> pN_evt2,
    bool E_plus_pz_onSimplex = false,
    double Sphere2Simplex = 4.0,
    bool l_has_chi = false,
    double Chi2Simplex = 16.0 * M_PI * M_PI
) {
    pN_evt1 = order_by_pt(pN_evt1);
    pN_evt2 = order_by_pt(pN_evt2);
    
    int N = pN_evt1.size();
    // std::cout << "N: " << N << std::endl;
    
    double E_evt1 = 0.0, pz_evt1 = 0.0;
    double E_evt2 = 0.0, pz_evt2 = 0.0;
    
    for (int i = 0; i < N; ++i) {
        E_evt1 += pN_evt1[i][0];
        pz_evt1 += pN_evt1[i][3];
        E_evt2 += pN_evt2[i][0];
        pz_evt2 += pN_evt2[i][3];
    }
    
    double chi_evt1 = (E_evt1 - pz_evt1) / (E_evt1 + pz_evt1);
    double chi_evt2 = (E_evt2 - pz_evt2) / (E_evt2 + pz_evt2);
    
    std::vector<double> uN_evt1(N), uN_evt2(N);
    std::vector<std::complex<double>> vN_evt1(N), vN_evt2(N);
    
    double sqrt_factor1, sqrt_factor2, sqrt_factor3, sqrt_factor4;
    
    if (E_plus_pz_onSimplex) {
        sqrt_factor1 = 1.0 / std::sqrt(E_evt1 + pz_evt1);
        sqrt_factor2 = 1.0 / std::sqrt(E_evt2 + pz_evt2);
        sqrt_factor3 = 1.0 / std::sqrt(E_evt1 - pz_evt1);
        sqrt_factor4 = 1.0 / std::sqrt(E_evt2 - pz_evt2);
        
        for (int i = 0; i < N; ++i) {
            uN_evt1[i] = sqrt_factor1 * std::sqrt(pN_evt1[i][0] + pN_evt1[i][3]);
            uN_evt2[i] = sqrt_factor2 * std::sqrt(pN_evt2[i][0] + pN_evt2[i][3]);
            
            vN_evt1[i] = sqrt_factor3 * std::complex<double>(pN_evt1[i][1], pN_evt1[i][2]) / 
                         std::sqrt(pN_evt1[i][0] + pN_evt1[i][3]);
            vN_evt2[i] = sqrt_factor4 * std::complex<double>(pN_evt2[i][1], pN_evt2[i][2]) / 
                         std::sqrt(pN_evt2[i][0] + pN_evt2[i][3]);
        }
    } else {
        sqrt_factor1 = 1.0 / std::sqrt(E_evt1 - pz_evt1);
        sqrt_factor2 = 1.0 / std::sqrt(E_evt2 - pz_evt2);
        sqrt_factor3 = 1.0 / std::sqrt(E_evt1 + pz_evt1);
        sqrt_factor4 = 1.0 / std::sqrt(E_evt2 + pz_evt2);
        
        for (int i = 0; i < N; ++i) {
            uN_evt1[i] = sqrt_factor1 * std::sqrt(pN_evt1[i][0] - pN_evt1[i][3]);
            uN_evt2[i] = sqrt_factor2 * std::sqrt(pN_evt2[i][0] - pN_evt2[i][3]);
            
            vN_evt1[i] = sqrt_factor3 * std::complex<double>(pN_evt1[i][1], pN_evt1[i][2]) / 
                         std::sqrt(pN_evt1[i][0] - pN_evt1[i][3]);
            vN_evt2[i] = sqrt_factor4 * std::complex<double>(pN_evt2[i][1], pN_evt2[i][2]) / 
                         std::sqrt(pN_evt2[i][0] - pN_evt2[i][3]);
        }
    }
    
    std::vector<double> rhoN_evt1 = calculate_power(uN_evt1, 2.0);
    std::vector<double> rhoN_evt2 = calculate_power(uN_evt2, 2.0);
    
    std::vector<std::complex<double>> vprimeNminus1_evt1(N-1);
    std::vector<std::complex<double>> vprimeNminus1_evt2(N-1);
    
    if (N == 2) {
        for (int i = 0; i < N-1; ++i) {
            vprimeNminus1_evt1[i] = vN_evt1[i] / std::sqrt(rhoN_evt1[1]);
            vprimeNminus1_evt2[i] = vN_evt2[i] / std::sqrt(rhoN_evt2[1]);
        }
    } else {
        for (int i = 0; i < N-1; ++i) {
            std::complex<double> vprime_i_evt1(0.0, 0.0);
            
            for (int j = 0; j < N-1; ++j) {
                if (i == j) {
                    double factor = std::pow(uN_evt1[i], 2) + uN_evt1[N-1] * 
                                   (1 - std::pow(uN_evt1[i], 2) - std::pow(uN_evt1[N-1], 2));
                    vprime_i_evt1 += factor * vN_evt1[i];
                } else {
                    double factor = uN_evt1[i] * uN_evt1[j] * (1 - uN_evt1[N-1]);
                    vprime_i_evt1 += factor * vN_evt1[j];
                }
            }
            
            vprime_i_evt1 /= (uN_evt1[N-1] * (1 - std::pow(uN_evt1[N-1], 2)));
            vprimeNminus1_evt1[i] = vprime_i_evt1;
        }
        
        for (int i = 0; i < N-1; ++i) {
            std::complex<double> vprime_i_evt2(0.0, 0.0);
            
            for (int j = 0; j < N-1; ++j) {
                if (i == j) {
                    double factor = std::pow(uN_evt2[i], 2) + uN_evt2[N-1] * 
                                   (1 - std::pow(uN_evt2[i], 2) - std::pow(uN_evt2[N-1], 2));
                    vprime_i_evt2 += factor * vN_evt2[i];
                } else {
                    double factor = uN_evt2[i] * uN_evt2[j] * (1 - uN_evt2[N-1]);
                    vprime_i_evt2 += factor * vN_evt2[j];
                }
            }
            
            vprime_i_evt2 /= (uN_evt2[N-1] * (1 - std::pow(uN_evt2[N-1], 2)));
            vprimeNminus1_evt2[i] = vprime_i_evt2;
        }
    }
    
    // Compute distances
    double l_chi = std::fabs(std::log(chi_evt1) - std::log(chi_evt2));
    
    std::vector<double> diff(N);
    for (int i = 0; i < N; ++i) {
        diff[i] = rhoN_evt2[i] - rhoN_evt1[i];
    }
    double l_simplex = std::sqrt(c_sum(calculate_power(diff, 2)));
    
    std::complex<double> dot_product1(0.0, 0.0);
    for (int i = 0; i < N-1; ++i) {
        dot_product1 += vprimeNminus1_evt1[i] * std::conj(vprimeNminus1_evt2[i]);
    }
    double dot_product1_abs = std::abs(dot_product1);
    
    if (dot_product1_abs > 1.0 && std::fabs(dot_product1_abs - 1.0) <= (1e-02 + 1e-05 * 1.0)) {
        dot_product1_abs = std::round(dot_product1_abs);
    }
    
    std::complex<double> dot_product2(0.0, 0.0);
    for (int i = 0; i < N-1; ++i) {
        dot_product2 += vprimeNminus1_evt1[i] * vprimeNminus1_evt2[i];
    }
    double dot_product2_abs = std::abs(dot_product2);
    
    if (dot_product2_abs > 1.0 && std::fabs(dot_product2_abs - 1.0) <= (1e-02 + 1e-05 * 1.0)) {
        dot_product2_abs = std::round(dot_product2_abs);
    }
    
    double l_sphere;
    if (dot_product1_abs <= 1.0 && dot_product2_abs <= 1.0) {
        l_sphere = std::acos(std::fmax(dot_product1_abs, dot_product2_abs));
    } else {
        l_sphere = -1.0;
    }
    
    double l;
    double OverallFactor;
    
    if (!l_has_chi) {
        OverallFactor = (1.0 / (4.0 * M_PI * M_PI * Sphere2Simplex)) * 
                       std::pow(Sphere2Simplex / 4.0, (N - 1.0) / (3.0 * N - 4.0));
        if (l_sphere >= 0) {
            l = std::sqrt(OverallFactor * (std::pow(l_simplex, 2) + Sphere2Simplex * std::pow(l_sphere, 2)));
        } else {
            l = -1.0;
        }
    } else {
        OverallFactor = (1.0 / (4.0 * M_PI * M_PI * Sphere2Simplex)) * 
                       std::pow(Sphere2Simplex / 4.0, N / (3.0 * N - 3.0)) * 
                       std::pow(Chi2Simplex / (16.0 * M_PI * M_PI), -1.0 / (3.0 * N - 3.0));
        if (l_sphere >= 0) {
            l = std::sqrt(OverallFactor * (std::pow(l_simplex, 2) + 
                                          Sphere2Simplex * std::pow(l_sphere, 2) + 
                                          Chi2Simplex * std::pow(l_chi, 2)));
        } else {
            l = -1.0;
        }
    }
    
    return {l, l_simplex, l_sphere, l_chi};
}

/*
 example running
 root -q -b readDelphes.C'("/Users/tav/Software/MG5_aMC_v3_5_3/vbfHToTauTau-M750/Events/run_01/tag_1_delphes_events.root")'

 */

GenParticle* getMother(TClonesArray *branchParticle, const GenParticle *particle);

void printMotherHistory(TClonesArray *branchParticle, const GenParticle *particle);

void readDelphes(const char *inputFile) {
  int debug = 0;
  gInterpreter->AddIncludePath("/Users/blackmac/Software/Delphes");
  gInterpreter->AddIncludePath("/Users/blackmac/Software/Delphes/external");
  gSystem->Load("/Users/blackmac/Software/Delphes/libDelphes.so");
  
  std::string fullpath = inputFile;
  std::istringstream iss(inputFile);
  std::string sampleName;
  
  std::getline(iss, sampleName, '/');
  int i = 0;
  while (i < 5) {
    std::getline(iss, sampleName, '/');
    ++i;
  }
  
  std::cout << "Sample used is " << sampleName << std::endl;
  

  TChain chain("Delphes");
  chain.Add(inputFile);
//  chain.Add("/Users/tav/Software/MG5_aMC_v3_5_3/vbfHToTauTau-M750/Events/run_01/tag_1_delphes_events.root");

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

    // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");

    // Histograms
  std::vector<TH1F*> hists;
  TH1F histJetInvMass("jet_mass", ";M_{inv}(j_{1}, j_{2}) [GeV];Jets / 5 GeV", 150, 0.0, 1500.0);
  hists.push_back(&histJetInvMass);
  TH1F histJetTauInvMass("jetTau_mass", ";M_{inv}(j_{tau, 1}, j_{tau, 2}) [GeV];Tau jets / 5 GeV", 150, 0.0, 1500.0);
  hists.push_back(&histJetTauInvMass);
  TH1F histJetPt("jet_pt", ";p_{T}(j) [GeV];Jets / 5 GeV", 150, 0.0, 1500.0);
  hists.push_back(&histJetPt);
  TH1F histJetTauPt("jetTau_pt", ";p_{T}(j_{tau}) [GeV];Jets / 5 GeV", 150, 0.0, 1500.0);
  hists.push_back(&histJetTauPt);
  TH1F histJetEta("jet_eta", ";#eta(j);Jets / bin", 60, -3., 3.);
  hists.push_back(&histJetEta);
  TH1F histJetTauEta("jetTau_eta", ";#eta(j_{tau}) ;Jets  / bin", 60, -3., 3.);
  hists.push_back(&histJetTauEta);
  TH1F histJetN("jet_N", ";N(j);Jets / 1", 15, -0.5, 14.5);
  hists.push_back(&histJetN);
  TH1F histJetTauN("jetTau_N", ";N(j_{tau});Jets  / bin", 15, -0.5, 14.5);
  hists.push_back(&histJetTauN);
  TH1F histGenZDauther("genZDauther", ";#tau daughter PDG IDs;Generator particle  / bin", 18, -0.5, 17.5);
  histGenZDauther.GetXaxis()->SetBinLabel(1,"ele");
  histGenZDauther.GetXaxis()->SetBinLabel(2,"#nu_{ele}");
  histGenZDauther.GetXaxis()->SetBinLabel(3,"mu");
  histGenZDauther.GetXaxis()->SetBinLabel(4,"#nu_{mu}");
  histGenZDauther.GetXaxis()->SetBinLabel(5,"#nu_{mu}");
  histGenZDauther.GetXaxis()->SetBinLabel(6,"#gamma");
  histGenZDauther.GetXaxis()->SetBinLabel(7,"#pi^{0}");
  histGenZDauther.GetXaxis()->SetBinLabel(8,"#pi^{-/+}");
  histGenZDauther.GetXaxis()->SetBinLabel(9,"K^{0}");
  histGenZDauther.GetXaxis()->SetBinLabel(10,"K^{+/-}");
  histGenZDauther.GetXaxis()->SetBinLabel(11,"K^{0}_{L}");
  histGenZDauther.GetXaxis()->SetBinLabel(12,"K^{0}_{S}");
  histGenZDauther.GetXaxis()->SetBinLabel(13,"#eta");
  histGenZDauther.GetXaxis()->SetBinLabel(14,"p^{+/-}"); 
  histGenZDauther.GetXaxis()->SetBinLabel(15,"#Sigma^{+/-}");
  histGenZDauther.GetXaxis()->SetBinLabel(16,"#Xi^{-/+}");
  histGenZDauther.GetXaxis()->SetBinLabel(17,"#Omega^{-/+}");

  hists.push_back(&histGenZDauther);
  TH1F histGenElectronPt("genElectron_pt", ";Generator Electron p_{T} [GeV];Generator Electrons / 1 GeV", 50, 0.0, 50.0);
  hists.push_back(&histGenElectronPt);

  TH1F histElectronPt("electron_pt", ";Electron p_{T} [GeV];Electrons / 1 GeV", 50, 0.0, 50.0);
  hists.push_back(&histElectronPt);
  TH1F histElectronInvariantMass("electron_invariantMass", ";Electron e^{+}e^{-} Invariant Mass [GeV];Electrons / 0.1 GeV", 100, 0.0, 200.0);
  hists.push_back(&histElectronInvariantMass);

  TH1F histEFlowTrackN("eflowTrack_N", ";EFlowTrack N;EFlowTracks / bin", 50, -0.5, 49.5);
  hists.push_back(&histEFlowTrackN);
  TH1F histEFlowTrackPt("eflowTrack_pt", ";EFlowTrack p_{T} [GeV];EFlowTracks / 1 GeV", 50, 0.0, 50.0);
  hists.push_back(&histEFlowTrackPt);
  TH1F histEFlowTrackEta("eflowTrack_eta", ";EFlowTrack #eta;EFlowTracks / bin", 60, -3., 3.);
  hists.push_back(&histEFlowTrackEta);
  TH1F histEflowTrackD0("eflowTrack_d0", ";EFlowTrack d_{0} [cm];EFlowTracks / bin", 10, -0.05, 0.05);
  hists.push_back(&histEflowTrackD0);
  TH1F histEflowTrackDz("eflowTrack_dz", ";EFlowTrack d_{z} [cm];EFlowTracks / bin", 10, -0.05, 0.05);
  hists.push_back(&histEflowTrackDz);
  TH1F histEFlowTrackPID("eflowTrack_PID", ";EFlowTrack PID;EFlowTracks / bin", 18, -0.5, 17.5);
  histEFlowTrackPID.GetXaxis()->SetBinLabel(1,"ele");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(2,"#nu_{ele}");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(3,"mu");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(4,"#nu_{mu}");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(5,"#nu_{mu}");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(6,"#gamma");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(7,"#pi^{0}");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(8,"#pi^{-/+}");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(9,"K^{0}");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(10,"K^{+/-}");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(11,"K^{0}_{L}");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(12,"K^{0}_{S}");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(13,"#eta");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(14,"p^{+/-}"); 
  histEFlowTrackPID.GetXaxis()->SetBinLabel(15,"#Sigma^{+/-}");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(16,"#Xi^{-/+}");
  histEFlowTrackPID.GetXaxis()->SetBinLabel(17,"#Omega^{-/+}");
  hists.push_back(&histEFlowTrackPID);
  TH1F histEflowInvariantMass("eflow_invariantMass", ";EFlowTrack Invariant Mass [GeV];EFlowTracks / 0.1 GeV", 100, 0.0, 200.0);
  hists.push_back(&histEflowInvariantMass);

  // Repeat all this for the track collection next:
  TH1F histTrackN("track_N", ";Track N;Tracks / bin", 50, -0.5, 49.5);
  hists.push_back(&histTrackN);
  TH1F histTrackPt("track_pt", ";Track p_{T} [GeV];Tracks / 1 GeV", 50, 0.0, 50.0);
  hists.push_back(&histTrackPt);
  TH1F histTrackEta("track_eta", ";Track #eta;Tracks / bin", 60, -3., 3.);
  hists.push_back(&histTrackEta);
  TH1F histTrackD0("track_d0", ";Track d_{0} [cm];Tracks / bin", 10, -0.05, 0.05);
  hists.push_back(&histTrackD0);
  TH1F histTrackDz("track_dz", ";Track d_{z} [cm];Tracks / bin", 10, -0.05, 0.05);
  hists.push_back(&histTrackDz);
  TH1F histTrackPID("track_PID", ";Track PID;Tracks / bin", 18, -0.5, 17.5);
  histTrackPID.GetXaxis()->SetBinLabel(1,"ele");
  histTrackPID.GetXaxis()->SetBinLabel(2,"#nu_{ele}");
  histTrackPID.GetXaxis()->SetBinLabel(3,"mu");
  histTrackPID.GetXaxis()->SetBinLabel(4,"#nu_{mu}");
  histTrackPID.GetXaxis()->SetBinLabel(5,"#nu_{mu}");
  histTrackPID.GetXaxis()->SetBinLabel(6,"#gamma");
  histTrackPID.GetXaxis()->SetBinLabel(7,"#pi^{0}");
  histTrackPID.GetXaxis()->SetBinLabel(8,"#pi^{-/+}");
  histTrackPID.GetXaxis()->SetBinLabel(9,"K^{0}");
  histTrackPID.GetXaxis()->SetBinLabel(10,"K^{+/-}");
  histTrackPID.GetXaxis()->SetBinLabel(11,"K^{0}_{L}");
  histTrackPID.GetXaxis()->SetBinLabel(12,"K^{0}_{S}");
  histTrackPID.GetXaxis()->SetBinLabel(13,"#eta");
  histTrackPID.GetXaxis()->SetBinLabel(14,"p^{+/-}"); 
  histTrackPID.GetXaxis()->SetBinLabel(15,"#Sigma^{+/-}"); 
  histTrackPID.GetXaxis()->SetBinLabel(16,"#Xi^{-/+}");
  histTrackPID.GetXaxis()->SetBinLabel(17,"#Omega^{-/+}");
  hists.push_back(&histTrackPID);
  TH1F histTrackInvariantMass("track_invariantMass", ";Track Invariant Mass [GeV];Tracks / 0.1 GeV", 100, 0.0, 200.0);
  hists.push_back(&histTrackInvariantMass);
  TH1F histTrackDr("track_dR", ";dR between tracks;Tracks / bin", 50, 0.0, 5.0);
  hists.push_back(&histTrackDr);
  TH1F histTrackDrLeading("track_dR_leading", ";dR between leading track and other tracks;Tracks / bin", 50, 0.0, 5.0);
  hists.push_back(&histTrackDrLeading);
  TH1F histTrackSameHemisphereDrLeading("track_sameHemisphere_dR_leading", ";dR between leading track and other tracks in same hemisphere;Tracks / bin", 50, 0.0, 5.0);
  hists.push_back(&histTrackSameHemisphereDrLeading);
  TH1F histTrackSameHemisphereN("track_sameHemisphere_N", ";Track N in same hemisphere as leading track;Tracks / bin", 50, -0.5, 49.5);
  hists.push_back(&histTrackSameHemisphereN);
  TH1F histTrackSameHemispherePt("track_sameHemisphere_pt", ";Track p_{T} in same hemisphere as leading track [GeV];Tracks / 1 GeV", 50, 0.0, 50.0);
  hists.push_back(&histTrackSameHemispherePt);
  TH1F histTrackSameHemisphereEta("track_sameHemisphere_eta", ";Track #eta in same hemisphere as leading track;Tracks / bin", 60, -3., 3.);
  hists.push_back(&histTrackSameHemisphereEta);
  TH1F histTrackSameHemispherePhi("track_sameHemisphere_phi", ";Track #phi in same hemisphere as leading track;Tracks / bin", 64, -3.2, 3.2);
  hists.push_back(&histTrackSameHemispherePhi);


  TH1F histTrackOppositeHemisphereDrLeading("track_oppositeHemisphere_dR_leading", ";dR between leading track and other tracks in opposite hemisphere;Tracks / bin", 50, 0.0, 5.0);
  hists.push_back(&histTrackOppositeHemisphereDrLeading);
  TH1F histTrackOppositeHemisphereN("track_oppositeHemisphere_N", ";Track N in opposite hemisphere to leading track;Tracks / bin", 50, -0.5, 49.5);
  hists.push_back(&histTrackOppositeHemisphereN);
  TH1F histTrackOppositeHemispherePt("track_oppositeHemisphere_pt", ";Track p_{T} in opposite hemisphere to leading track [GeV];Tracks / 1 GeV", 50, 0.0, 50.0);
  hists.push_back(&histTrackOppositeHemispherePt);
  TH1F histTrackOppositeHemisphereEta("track_oppositeHemisphere_eta", ";Track #eta in opposite hemisphere to leading track;Tracks / bin", 60, -3., 3.);
  hists.push_back(&histTrackOppositeHemisphereEta);
  TH1F histTrackOppositeHemispherePhi("track_oppositeHemisphere_phi", ";Track #phi in opposite hemisphere to leading track;Tracks / bin", 64, -3.2, 3.2);
  hists.push_back(&histTrackOppositeHemispherePhi);

  std::vector<TH2F*> hists_2d;
  TH2F histTrackInvariantMassVsDR("track_invariantMass_vs_dR", ";Track Invariant Mass [GeV];dR between tracks", 100, 0.0, 200.0, 50, 0.0, 5.0);
  hists_2d.push_back(&histTrackInvariantMassVsDR);
  TH2F histTrackInvariantMassVsPt("track_invariantMass_vs_Pt", ";Track Invariant Mass [GeV];p_{T} of leading track [GeV]", 100, 0.0, 200.0, 50, 0.0, 50.0);
  hists_2d.push_back(&histTrackInvariantMassVsPt);

  
  std::cout << "Running on " << numberOfEntries << " events" << std::endl;
  // Min track pt requirement
  const double cut_minTrackPt_ = 5.0; // GeV

  // numberOfEntries = 10000;
  numberOfEntries = 10;
  for (Long64_t entry = 0; entry < numberOfEntries; ++entry) {
    if (debug > 0) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Event " << entry << std::endl;
    }

    if (entry % 2000 == 0) {
      std:cout << "Processing event " << entry << std::endl;
    }
    treeReader->ReadEntry(entry);

    // // -------------------------------------------------------------------------------
    // // Loop on the jets
    // bool filled = false;
    // bool filledTau = false;
    // int numTauJets = 0;
    // int numJets = 0;
    // for (int i = 0; i < branchJet->GetEntries(); ++i) {
    //   Jet *jet = (Jet*) branchJet->At(i);
    //   if (!jet) continue;
    //   ++numJets;
    //   histJetPt.Fill(jet->PT);
    //   histJetEta.Fill(jet->Eta);

    //   if (jet->TauTag == 1) {
    //     ++numTauJets;
    //     histJetTauPt.Fill(jet->PT);
    //     histJetTauEta.Fill(jet->Eta);
    //   }
    // } // end loop on jets

    // histJetTauN.Fill(numTauJets);
    // histJetN.Fill(branchJet->GetEntries());    

    // -------------------------------------------------------------------------------
    // Loop on the gen particles
    int numGenPart = 0;
    for (int i = 0; i < branchParticle->GetEntries(); ++i) {
      GenParticle *gen = (GenParticle*) (branchParticle->At(i));
      // dont look at non-final state genparts
      if (!gen) continue;
      if (gen->Status != 1)  continue;
      // count the remaining gen parts
      ++numGenPart;
      if (debug > 2) {
        std::cout << "Status = " <<  gen->Status << ", PID: " << gen->PID <<  std::endl;
        printMotherHistory(branchParticle, gen);
      }
      GenParticle *genMom = getMother(branchParticle, gen);
      if (genMom==0) continue;
      if (abs(genMom->PID) == 23) {
        int genPID = abs(gen->PID);
        if (genPID == 11) { // ele
          histGenZDauther.Fill(0);
          histGenElectronPt.Fill(gen->PT);
        } else if (genPID == 12) { // nu ele
          histGenZDauther.Fill(1);
        } else if (genPID == 13) { // mu
          histGenZDauther.Fill(2);
        } else if (genPID == 14) { // nu mu
          histGenZDauther.Fill(3);
        } else if (genPID == 16) { // nu tau
          histGenZDauther.Fill(4);
        }  else if (genPID == 22) { // gamma
          histGenZDauther.Fill(5);
        } else if (genPID == 111) { // pi0
          histGenZDauther.Fill(6);
        } else if (genPID == 211) { // pi+/-
          histGenZDauther.Fill(7);
        } else if (genPID == 311) { // K0
          histGenZDauther.Fill(8);
        } else if (genPID == 321 || genPID == 323) { // K+
          histGenZDauther.Fill(9);
        } else if (genPID == 130) { // KL0
          histGenZDauther.Fill(10);
        } else if (genPID == 310) { // Ks
          histGenZDauther.Fill(11);
        } else if (genPID == 221) { // eta
          histGenZDauther.Fill(12);
        } else if (genPID == 2212) { // proton
          histGenZDauther.Fill(13);
        } else if (genPID == 3222 || genPID == 3112) { // sigma+/-
          histGenZDauther.Fill(14);
        } else if (genPID == 3312) { // xi-/+
          histGenZDauther.Fill(15);
        } else if (genPID == 3334) { // omega-/+
          histGenZDauther.Fill(16);
        } else {
          std::cout << "This particle comes from a tau decay, PID is " << gen->PID <<  std::endl;
        }
      } // if from Z
      // Print statements here
      // print(,gen.PID,", E: ",gen.E,", PT: ",gen.PT,", Eta: ",gen.Eta,", M: ",gen.Mass,", M1: ",gen.M1,", M2: ",gen.M2,", D1: ",gen.D1,", D2: ",gen.D2)
    } // end lopp on genPart

    // -------------------------------------------------------------------------------
    // start loop on EFlowTrack
    int numEFlowTrack = 0;
    for (int i = 0; i < branchEFlowTrack->GetEntries(); ++i) {
      ParticleFlowCandidate *eflowTrack = (ParticleFlowCandidate*) (branchEFlowTrack->At(i));
      if (!eflowTrack) continue;
      auto eflowTrackPt = eflowTrack->PT;
      auto eflowTrackEta = eflowTrack->Eta;
      auto eflowTrackD0 = eflowTrack->D0;
      auto eflowTrackDz = eflowTrack->DZ;
      auto eflowTrackCharge = eflowTrack->Charge;
      auto eflowTrackPID = abs(eflowTrack->PID);
      if (eflowTrackPt < cut_minTrackPt_) continue;
      if (std::abs(eflowTrackEta) > 2.5) continue;
      if (std::abs(eflowTrackD0) > 0.05) continue;
      if (std::abs(eflowTrackDz) > 0.05) continue;
      // fill PID histogram
      if (eflowTrackPID == 11) { // ele
        histEFlowTrackPID.Fill(0);
      } else if (eflowTrackPID == 12) { // nu ele
        histEFlowTrackPID.Fill(1);
      } else if (eflowTrackPID == 13) { // mu
        histEFlowTrackPID.Fill(2);
      } else if (eflowTrackPID == 14) { // nu mu
        histEFlowTrackPID.Fill(3);
      } else if (eflowTrackPID == 16) { // nu tau
        histEFlowTrackPID.Fill(4);
      }  else if (eflowTrackPID == 22) { // gamma
        histEFlowTrackPID.Fill(5);
      } else if (eflowTrackPID == 111) { // pi0
        histEFlowTrackPID.Fill(6);
      } else if (eflowTrackPID == 211) { // pi+/-
        histEFlowTrackPID.Fill(7);
      } else if (eflowTrackPID == 311) { // K0
        histEFlowTrackPID.Fill(8);
      } else if (eflowTrackPID == 321 || eflowTrackPID == 323) { // K+
        histEFlowTrackPID.Fill(9);
      } else if (eflowTrackPID == 130) { // KL0
        histEFlowTrackPID.Fill(10);
      } else if (eflowTrackPID == 310) { // Ks
        histEFlowTrackPID.Fill(11);
      } else if (eflowTrackPID == 221) { // eta
        histEFlowTrackPID.Fill(12);
      } else if (eflowTrackPID == 2212) { // proton
        histEFlowTrackPID.Fill(13);
      } else if (eflowTrackPID == 3222 || eflowTrackPID == 3112) { // sigma+/-
        histEFlowTrackPID.Fill(14);
      } else if (eflowTrackPID == 3312) { // xi-/+
        histEFlowTrackPID.Fill(15);
      } else if (eflowTrackPID == 3334) { // omega-/+
        histEFlowTrackPID.Fill(16);
      } else {
        std::cout << "PF PID is " << eflowTrackPID <<  std::endl;
      }
      if (eflowTrackPID != 11) continue;
      histEFlowTrackPt.Fill(eflowTrackPt);
      histEFlowTrackEta.Fill(eflowTrackEta);
      histEflowTrackD0.Fill(eflowTrackD0);
      histEflowTrackDz.Fill(eflowTrackDz);
      ++numEFlowTrack;
      // Loop again on the eflow to measure invariant mass of pairs
      for (int j = 0; j < branchEFlowTrack->GetEntries(); ++j) {
        if (j == i) continue;
        ParticleFlowCandidate *eflowTrack2 = (ParticleFlowCandidate*) (branchEFlowTrack->At(j));
        if (!eflowTrack2) continue;
        auto eflowTrack2Pt = eflowTrack2->PT;
        auto eflowTrack2Eta = eflowTrack2->Eta;
        auto eflowTrack2D0 = eflowTrack2->D0;
        auto eflowTrack2Dz = eflowTrack2->DZ;
        auto eflowTrack2Charge = eflowTrack2->Charge;
        auto eflowTrack2PID = abs(eflowTrack2->PID);
        if (eflowTrack2Pt < cut_minTrackPt_) continue;
        if (std::abs(eflowTrack2Eta) > 2.5) continue;
        if (std::abs(eflowTrack2D0) > 0.05) continue;
        if (std::abs(eflowTrack2Dz) > 0.05) continue;
        if (eflowTrack2PID != 11) continue;
        if (eflowTrackCharge == eflowTrack2Charge) continue;
        TLorentzVector p4sum = eflowTrack->P4() + eflowTrack2->P4();
        histEflowInvariantMass.Fill(p4sum.M());
        break;
        // fill histograms here if needed
      }
    } // end loop on EFlowTrack
    histEFlowTrackN.Fill(numEFlowTrack);

    // -------------------------------------------------------------------------------
    // Start loop on Electrons
    for (int i = 0; i < branchElectron->GetEntries(); ++i) {
      Electron *electron = (Electron*) (branchElectron->At(i));
      if (!electron) continue;
      histElectronPt.Fill(electron->PT);
      // Loop again on electrons to measure invariant mass of pairs
      for (int j = 0; j < branchElectron->GetEntries(); ++j) {
        if (j == i) continue;
        Electron *electron2 = (Electron*) (branchElectron->At(j));
        if (!electron2) continue;
        if (electron->Charge == electron2->Charge) continue;
        TLorentzVector p4sum = electron->P4() + electron2->P4();
        histElectronInvariantMass.Fill(p4sum.M());
        break;
        // fill histograms here if needed
      } // end loop on electrons for invariant mass
    } // end loop on Electrons


    // -------------------------------------------------------------------------------
    // Start loop on tracks
    std::vector<Track*> sorted_tracks;
    for (int i = 0; i < branchTrack->GetEntries(); ++i) {
        Track* track = (Track*)branchTrack->At(i);
        if (track) sorted_tracks.push_back(track);
    }
    // sort tracks by pt
    std::sort(sorted_tracks.begin(), sorted_tracks.end(),
              [](Track* a, Track* b) { return a->PT > b->PT; });

    int numTrack = 0;
    int numTrackSameHemisphere = 0;
    int numTrackOppositeHemisphere = 0;
    double leadingTrackEta = -99.0;
    double leadingTrackPhi = -99.0;
    for (int i = 0; i < sorted_tracks.size(); ++i) {
      Track *track = sorted_tracks[i];
      if (!track) continue;
      auto trackPt = track->PT;
      if (trackPt < cut_minTrackPt_) continue;
      auto trackEta = track->Eta;
      auto trackPhi = track->Phi;
      if (std::abs(trackEta) > 2.5) continue;
      auto trackD0 = track->D0;
      if (std::abs(trackD0) > 0.05) continue;
      auto trackDz = track->DZ;
      if (std::abs(trackDz) > 0.05) continue;
      if (leadingTrackEta == -99.0) leadingTrackEta = trackEta;
      if (leadingTrackPhi == -90.0) leadingTrackPhi = trackPhi;
      auto trackCharge = track->Charge;
      // fill PID histogram
      auto trackPID = abs(track->PID);
      if (trackPID == 11) { // ele
        histTrackPID.Fill(0);
      } else if (trackPID == 12) { // nu ele
        histTrackPID.Fill(1);
      } else if (trackPID == 13) { // mu
        histTrackPID.Fill(2);
      } else if (trackPID == 14) { // nu mu
        histTrackPID.Fill(3);
      } else if (trackPID == 16) { // nu tau
        histTrackPID.Fill(4);
      }  else if (trackPID == 22) { // gamma
        histTrackPID.Fill(5);
      } else if (trackPID == 111) { // pi0
        histTrackPID.Fill(6);
      } else if (trackPID == 211) { // pi+/-
        histTrackPID.Fill(7);
      } else if (trackPID == 311) { // K0
        histTrackPID.Fill(8);
      } else if (trackPID == 321 || trackPID == 323) { // K+
        histTrackPID.Fill(9);
      } else if (trackPID == 130) { // KL0
        histTrackPID.Fill(10);
      } else if (trackPID == 310) { // Ks
        histTrackPID.Fill(11);
      } else if (trackPID == 221) { // eta
        histTrackPID.Fill(12);
      } else if (trackPID == 2212) { // proton
        histTrackPID.Fill(13);
      } else if (trackPID == 3222 || trackPID == 3112)  { // sigma+/-
        histTrackPID.Fill(14);
      } else if (trackPID == 3312) { // xi-/+
        histTrackPID.Fill(15);
      } else if (trackPID == 3334) { // omega-/+
        histTrackPID.Fill(16);
      } else {
        std::cout << "Track PID is " << trackPID <<  std::endl;
      }
      double deltaEta_leadingTrack = leadingTrackEta - trackEta;
      double deltaPhi_leadingTrack = TVector2::Phi_mpi_pi(leadingTrackPhi - trackPhi);
      double dR_leadingTrack = std::sqrt(deltaEta_leadingTrack * deltaEta_leadingTrack + deltaPhi_leadingTrack * deltaPhi_leadingTrack);
      histTrackDrLeading.Fill(dR_leadingTrack);
      if (deltaPhi_leadingTrack < TMath::Pi()/2.0) { // same hemisphere
        histTrackSameHemisphereDrLeading.Fill(dR_leadingTrack);
        histTrackSameHemispherePt.Fill(trackPt);
        histTrackSameHemisphereEta.Fill(trackEta);
        histTrackSameHemispherePhi.Fill(trackPhi);
        numTrackSameHemisphere++;
      } else {
        histTrackOppositeHemisphereDrLeading.Fill(dR_leadingTrack);
        histTrackOppositeHemispherePt.Fill(trackPt);
        histTrackOppositeHemisphereEta.Fill(trackEta);
        histTrackOppositeHemispherePhi.Fill(trackPhi);
        numTrackOppositeHemisphere++;
      }
      // if (trackPID != 11) continue;
      histTrackPt.Fill(trackPt);
      histTrackEta.Fill(trackEta);
      histTrackD0.Fill(trackD0);
      histTrackDz.Fill(trackDz);
      ++numTrack;
      // Loop again on the tracks to measure invariant mass of pairs
      for (int j = 0; j < sorted_tracks.size(); ++j) {
        if (j == i) continue;
        Track *track2 = sorted_tracks[j];
        if (!track2) continue;
        auto track2Pt = track2->PT;
        if (track2Pt < cut_minTrackPt_) continue;
        auto track2Eta = track2->Eta;
        if (std::abs(track2Eta) > 2.5) continue;
        auto track2D0 = track2->D0;
        if (std::abs(track2D0) > 0.05) continue;
        auto track2Dz = track2->DZ;
        if (std::abs(track2Dz) > 0.05) continue;
        auto track2PID = abs(track2->PID);
        // if (track2PID != 11) continue;
        auto track2Charge = track2->Charge;
        if (trackCharge == track2Charge) continue;
        // Calculate dR distance between the two tracks
        double deltaEta = trackEta - track2Eta;
        double deltaPhi = TVector2::Phi_mpi_pi(track->Phi - track2->Phi);
        double dR = std::sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
        histTrackDr.Fill(dR);
        if (dR < 0.5) continue; // avoid double counting the same track
        // calculate invariant mass
        TLorentzVector p4sum = track->P4() + track2->P4();
        histTrackInvariantMass.Fill(p4sum.M());
        histTrackInvariantMassVsDR.Fill(p4sum.M(), dR);
        histTrackInvariantMassVsPt.Fill(p4sum.M(), trackPt);
        break;
        // fill histograms here if needed
      } // end loop on tracks for invariant mass
    } // end loop on tracks
    histTrackN.Fill(numTrack);
    histTrackSameHemisphereN.Fill(numTrackSameHemisphere);
    histTrackOppositeHemisphereN.Fill(numTrackOppositeHemisphere);






    // -------------------------------------------------------------------------------
  } // end loop on events

  std::vector<std::vector<double>> event1 = {
      {100.0, 10.0, 20.0, 95.0},  // [E, px, py, pz]
      {80.0, -5.0, 15.0, 78.0}
  };

    std::vector<std::vector<double>> event2 = {
      {100.0, 10.0, 20.0, 95.0},  // [E, px, py, pz]
      {80.0, -5.0, 15.0, 78.0}
  };

  // std::vector<std::vector<double>> event2 = {
  //     {105.0, 12.0, 18.0, 98.0},
  //     {75.0, -6.0, 14.0, 72.0}
  // };

  auto distances = psDistance(event1, event2);
  std::cout << "Phase space distances between event 1 and event 2:"
  << "\n\tl(total distance): " <<  distances[0] << ",\n\tl(simplex): " <<
  distances[1] << ",\n\tl(sphere): " <<
  distances[2] << ",\n\tl(chi): " <<
  distances[3] << std::endl;

  std::string outputFileName = "./Histos_" + sampleName + ".root";
  TFile outputFile(outputFileName.c_str(), "RECREATE");
  // Write histograms to file
  for (auto hist : hists) {
    hist->Write();
  }
  for (auto hist2d : hists_2d) {
    hist2d->Write();
  }
  outputFile.Close();
} // end of main function











// Auxiliary functions
// -------------------------------------------------------------------------------
// Functions to get mother particle
GenParticle* getMother(TClonesArray *branchParticle, const GenParticle *particle){
  GenParticle *genMother = NULL;
  if (particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }
  
    // Is this the first parent with a different ID? If yes, return, otherwise
    // go deeper into recursion
  if (particle->M1 > 0 && particle->PID != 0) {
    genMother = (GenParticle*) (branchParticle->At(particle->M1));
    if (genMother->PID ==  particle->PID) {
      return getMother(branchParticle, genMother);
    } else {
      return genMother;
    }
  }
  else {
    return NULL;
  }
}

// -------------------------------------------------------------------------------
// Print mother history
void printMotherHistory(TClonesArray *branchParticle, const GenParticle *gen){
  int lastPID = 0;
  GenParticle *genCurrent =  (GenParticle*)  gen;
  while (abs(lastPID) != 2212 || (abs(lastPID) < 6) ) {
    GenParticle *genTemp = getMother(branchParticle, genCurrent);
    if (!genTemp) break;
    lastPID = abs(genTemp->PID);
    genCurrent = genTemp;
    std::cout << "  > Mother ID = " << lastPID << " status = " << genTemp->Status << endl;
    if (lastPID == 23 )  std::cout << "    >> Found the Z!" << endl;
  }
}
