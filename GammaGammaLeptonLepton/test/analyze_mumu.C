#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test/Canvas.h"

#include "CTPPSAnalysisTools/Alignment/interface/AlignmentsFactory.h"
#include "CTPPSAnalysisTools/Reconstruction/interface/LHCConditionsFactory.h"
#include "CTPPSAnalysisTools/Reconstruction/interface/XiReconstructor.h"

#include "TH1.h"
#include "THStack.h"
#include "TGraph.h"
#include "TFile.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TPie.h"
#include "TStyle.h"

#include <iostream>
#include <unordered_map>

#define WWW_PATH "/afs/cern.ch/user/l/lforthom/www/private/dm/mumu/tmp/"
//#define WWW_PATH "/afs/cern.ch/user/l/lforthom/www/private/dm/mumu/highxicut/"

using namespace std;

enum class Region { noZ = 0, Z = 1, excess = 2, full = 3 };
std::ostream& operator<<( std::ostream& os, const Region& reg ) {
  switch ( reg ) {
    case Region::Z: return os << "Z region";
    case Region::noZ: return os << "outside Z region";
    case Region::excess: return os << "excess region";
    case Region::full: return os << "#it{#font[12]{l}}^{+}#it{#font[12]{l}}^{-} preselection";
  }
  return os;
}

void analyze_mumu()
{
  const vector<double> miss_mass_bins = { 200., 400., 800., 1000., 1200., 1400., 1800., 2500. };
  const vector<double> xangles = { 120., 130., 140., 150. };

  //const double min_z_mass = 76., max_z_mass = 110.; // "my" values
  const double min_z_mass = 88., max_z_mass = 94.; // Diego/Nicola values

  vector<TH1D>
    h_m_pair_all( 3, TH1D( "h_m_pair_all", "m(#mu^{+}#mu^{-})\\Events\\GeV", 70, 50., 400. ) ),
    h_pt_pair_all( 3, TH1D( "h_pt_pair_all", "p_{T}(#mu^{+}#mu^{-})\\Events\\GeV", 50, 25., 275. ) ),
    h_pp_mass_all( 3, TH1D( "h_pp_mass_all", "m_{pp}\\Events\\GeV", 60, 100., 2500. ) ),
    h_rap_pair_all( 3, TH1D( "h_rap_pair_all", "Y(#mu^{+}#mu^{-})\\Events", 20, -3., 3. ) ),
    h_acop_all( 3, TH1D( "h_acop_all", "1-|#Delta#phi(#mu^{+}#mu^{-})/#pi|\\Events", 50, 0., 1. ) ),
    h_xangle_all( 4, TH1D( "h_xangle_all", "Crossing angle\\Events\\#murad", 50, 100., 200. ) ),
    //-----
    h_m_pair( 2, TH1D( "h_m_pair", "m(#mu^{+}#mu^{-})\\Events\\GeV", 70, 50., 400. ) ),
    h_rap_pair( 2, TH1D( "h_rap_pair", "Y(#mu^{+}#mu^{-})\\Events", 20, -3., 3. ) ),
    h_pt_pair( 2, TH1D( "h_pt_pair", "p_{T}(#mu^{+}#mu^{-})\\Events\\GeV", 50, 25., 275. ) ),
    h_pt_leadmu( 2, TH1D( "h_pt_leadmu", "Leading #mu^{#pm} p_{T}\\Events\\GeV", 50, 0., 250. ) ),
    h_pt_subleadmu( 2, TH1D( "h_pt_subleadmu", "Subleading #mu^{#pm} p_{T}\\Events\\GeV", 50, 0., 250. ) ),
    h_eta_leadmu( 2, TH1D( "h_eta_leadmu", "Leading #mu^{#pm} #eta\\Events", 50, -2.5, 2.5 ) ),
    h_eta_subleadmu( 2, TH1D( "h_eta_subleadmu", "Subleading #mu^{#pm} #eta\\Events", 50, -2.5, 2.5 ) ),
    h_acop( 2, TH1D( "h_acop", "1-|#Delta#phi(#mu^{+}#mu^{-})/#pi|\\Events", 50, 0., 1. ) ),
    h_miss_mass( 2, TH1D( "h_miss_mass", "Missing mass\\Events\\GeV", 63, 0., 2520. ) ),
    h_miss_mass_zoom( 2, TH1D( "h_miss_mass_zoom", "Missing mass\\Events\\GeV", 25, 1150., 1650. ) ),
    h_pp_mass( 2, TH1D( "h_pp_mass", "m_{pp}\\Events\\GeV", 60, 100., 2500. ) ),
    //-----
    h_pp_mass_permmb( miss_mass_bins.size()-1, TH1D( "h_pp_mass_permmb", "m_{pp}\\Events\\GeV", 56, 200., 3000. ) );
  TGraph g_xi1xi2;

  //----- x alignment tools
  ctpps::AlignmentsFactory align_fac;
  ostringstream align_file_path;
  //align_file_path << getenv( "CMSSW_BASE" ) << "/src/CTPPSAnalysisTools/Alignment/data/2017/alignments_30jan2017.txt";
  align_file_path << getenv( "CMSSW_BASE" ) << "/src/CTPPSAnalysisTools/Alignment/data/2017/alignments_21aug2018.txt";
  align_fac.feedAlignments( align_file_path.str().c_str() );

  //----- crossing angle (+ other useful LHC information) tools
  ctpps::LHCConditionsFactory cond_fac;
  ostringstream cond_file1_path, cond_file2_path;
  cond_file1_path << getenv( "CMSSW_BASE" ) << "/src/CTPPSAnalysisTools/Reconstruction/data/2017/xangle_tillTS2.csv";
  cond_file2_path << getenv( "CMSSW_BASE" ) << "/src/CTPPSAnalysisTools/Reconstruction/data/2017/xangle_afterTS2.csv";
  cond_fac.feedConditions( cond_file1_path.str().c_str() );
  cond_fac.feedConditions( cond_file2_path.str().c_str() );

  //----- x-to-xi reconstruction tools
  ctpps::XiReconstructor reco;
  ostringstream disp_file_path; disp_file_path << getenv( "CMSSW_BASE" ) << "/src/CTPPSAnalysisTools/Reconstruction/data/2017/dispersions.txt";
  reco.feedDispersions( disp_file_path.str().c_str() );
  reco.interpolateCrossingAngle( 150., 0.46 );

  double int_lumi = 0.;

  //auto file = TFile::Open( filename );
  vector<const char*> runs_list;
  TChain ch( "analysis/ntp1" );
  ch.Add( "/eos/cms/store/group/phys_pps/dilepton/DoubleMuon/mumu_ntuple-Run2017B-ptgt50_v1.root" ); runs_list.emplace_back( "B" ); int_lumi += 2.367133065086;
  ch.Add( "/eos/cms/store/group/phys_pps/dilepton/DoubleMuon/mumu_ntuple-Run2017C-ptgt50_v1_partial.root" ); runs_list.emplace_back( "C" ); int_lumi += 8.680670867306;
  ch.Add( "/eos/cms/store/group/phys_pps/dilepton/DoubleMuon/mumu_ntuple-Run2017D-ptgt50_v1_partial.root" ); runs_list.emplace_back( "D" ); int_lumi += 4.142654224602;
  ch.Add( "/eos/cms/store/group/phys_pps/dilepton/DoubleMuon/mumu_ntuple-Run2017E-ptgt50_v1_partial.root" ); runs_list.emplace_back( "E" ); int_lumi += 9.061579665808;
  ch.Add( "/eos/cms/store/group/phys_pps/dilepton/DoubleMuon/mumu_ntuple-Run2017F-ptgt50_v1_partial.root" ); runs_list.emplace_back( "F" ); int_lumi += 13.920866597190;

  PaveText pt_runs( 0.82, 0.64, 0.85, 0.73 );
  pt_runs.AddText( "#mu^{+}#mu^{-} channel" );
  pt_runs.AddText( runs_list.size() > 1 ? Form( "Runs %s-%s", *runs_list.begin(), *runs_list.rbegin() ) : Form( "Run %s", *runs_list.begin() ) );

  const unsigned short num_samples = ch.GetListOfFiles()->GetEntries();
  vector<map<double,unsigned int> > reparts_per_xangle( num_samples, map<double,unsigned int>{ { -1., 0 } } );

  //ggll::AnalysisEvent out_evt;
  //auto out_file = TFile::Open( "skimmed_dset.root", "recreate" );
  //TTree out_tree( "analysis/ntp2", "A preskimmed analysis tree" );
  //out_evt.attach( &out_tree, ggll::DiMuon, false );

  TFile *current_file = nullptr;
  TIter fileIter( ch.GetListOfFiles() );
  unsigned short file_id = 0;
  while ( ( current_file = (TFile*)fileIter.Next() ) ) {
    auto file = TFile::Open( current_file->GetTitle() );
    cout << "opening " << current_file->GetTitle() << endl;
    auto tree = dynamic_cast<TTree*>( file->Get( "analysis/ntp1" ) );
    tree->SetBranchStatus( "*", false );
    tree->SetBranchStatus( "Run", 1 ); tree->SetBranchStatus( "LumiSection", 1 ); tree->SetBranchStatus( "EventNum", 1 );
    tree->SetBranchStatus( "nHLT", 1 ); /*tree->SetBranchStatus( "HLT_Name", 1 );*/ tree->SetBranchStatus( "HLT_Accept", 1 ); tree->SetBranchStatus( "HLT_Prescl", 1 );
    tree->SetBranchStatus( "nLocalProtCand", 1 );
    tree->SetBranchStatus( "LocalProtCand_arm", 1 ); tree->SetBranchStatus( "LocalProtCand_station", 1 ); tree->SetBranchStatus( "LocalProtCand_pot", 1 );
    tree->SetBranchStatus( "LocalProtCand_x", 1 );
    tree->SetBranchStatus( "nPair", 1 );
    tree->SetBranchStatus( "Pair_lepton1", 1 ); tree->SetBranchStatus( "Pair_lepton2", 1 ); tree->SetBranchStatus( "Pair_mass", 1 );
    tree->SetBranchStatus( "Pair_pt", 1 ); tree->SetBranchStatus( "Pair_dphi", 1 );
    tree->SetBranchStatus( "MuonCand_pt", 1 ); tree->SetBranchStatus( "MuonCand_eta", 1 ); tree->SetBranchStatus( "MuonCand_phi", 1 ); tree->SetBranchStatus( "MuonCand_e", 1 );
    tree->SetBranchStatus( "MuonCand_charge", 1 ); tree->SetBranchStatus( "MuonCand_istight", 1 );
    ggll::AnalysisEvent evt;
    evt.load( tree, ggll::DiMuon, false );

    const unsigned long long num_entries = tree->GetEntriesFast()/1;
    for ( unsigned long long i = 0; i < num_entries; ++i ) {
      tree->GetEntry( i );
      if ( i % 25000 == 0 )
        cout << "processing event " << i << " / " << num_entries << endl;

      bool passes_trigger = false;
      for ( unsigned int j = 0; j < evt.nHLT; ++j )
        if ( evt.HLT_Accept[j] == 1 && evt.HLT_Prescl[j] == 1. )
          passes_trigger = true;
      if ( !passes_trigger ) continue;

      const edm::EventID event_id( evt.Run, evt.LumiSection, evt.EventNum );
      const auto cond = cond_fac.get( event_id );

      h_xangle_all[(unsigned short)Region::full].Fill( cond.crossing_angle );
      if ( find( xangles.begin(), xangles.end(), cond.crossing_angle ) != xangles.end() )
        reparts_per_xangle[file_id][cond.crossing_angle] += 1;
      else
        reparts_per_xangle[file_id][-1.] += 1;

      unsigned short num_arm0 = 0, num_arm1 = 0;
      double xi_arm0 = -1., xi_arm1 = -1.;
      for ( unsigned int j = 0; j < evt.nLocalProtCand; ++j ) {
        if ( num_arm0 > 1 || num_arm1 > 1 ) break; //FIXME
        //if ( evt.LocalProtCand_pot[j] == 6 ) continue;
        //if ( evt.LocalProtCand_station[j] != 0 ) continue; // near pot
        if ( evt.LocalProtCand_pot[j] != 3 ) continue; // only keep pixels
        double xi = -1., xi_err = -1.;
        const unsigned short pot_id = evt.LocalProtCand_arm[j]*100+evt.LocalProtCand_station[j]*10+evt.LocalProtCand_pot[j];
        const auto align = align_fac.get( event_id, pot_id );
        const double x_aligned = evt.LocalProtCand_x[j]*100.+align.x_align; // in cm
        try {
          reco.reconstruct( cond.crossing_angle, pot_id, x_aligned, xi, xi_err );
        } catch ( const std::exception& e ) {
          //cerr << ">>>" << e.what() << endl;
          continue;
        }
        if ( evt.LocalProtCand_arm[j] == 0 ) { // sector 4-5
//          if ( xi > 0.2 ) continue;
          xi_arm0 = xi;
          ++num_arm0;
        }
        if ( evt.LocalProtCand_arm[j] == 1 ) { // sector 5-6
//          if ( xi > 0.25 ) continue;
          xi_arm1 = xi;
          ++num_arm1;
        }
      }
      //--- only keep events with at most 1 track per pot/side
      if ( num_arm0 != 1 || num_arm1 != 1 ) continue;
      const TLorentzVector pp1( 0., 0., -6500.*xi_arm0, 6500.*xi_arm0 ), pp2( 0., 0., 6500.*xi_arm1, 6500.*xi_arm1 );
      //--- only keep events with at most 1 lepton pair
      unsigned short dil_cand = 0;
      TLorentzVector pl1, pl2;
      unsigned short r = 100;
      int id_pair = -1;
      double max_mass = -1.;
      for ( unsigned int j = 0; j < evt.nPair; ++j ) {
        if ( evt.Pair_pt[j] < 60. ) continue;
        const unsigned int l1 = evt.Pair_lepton1[j], l2 = evt.Pair_lepton2[j];
        if ( evt.MuonCand_pt[l1] < 35. || evt.MuonCand_pt[l2] < 35. ) continue;
        if ( fabs( evt.MuonCand_eta[l1] ) > 2.4 || fabs( evt.MuonCand_eta[l2] ) > 2.4 ) continue;
        if ( evt.MuonCand_charge[l1]*evt.MuonCand_charge[l2] >= 0 ) continue;
        if ( !evt.MuonCand_istight[l1] || !evt.MuonCand_istight[l2] ) continue;

        if ( max_mass > 0. && evt.Pair_mass[j] < max_mass ) continue;
        max_mass = evt.Pair_mass[j];

        r = (unsigned short)( ( evt.Pair_mass[j] < min_z_mass || evt.Pair_mass[j] > max_z_mass ) ? Region::noZ : Region::Z );
        pl1.SetPtEtaPhiE( evt.MuonCand_pt[l1], evt.MuonCand_eta[l1], evt.MuonCand_phi[l1], evt.MuonCand_e[l1] );
        pl2.SetPtEtaPhiE( evt.MuonCand_pt[l2], evt.MuonCand_eta[l2], evt.MuonCand_phi[l2], evt.MuonCand_e[l2] );

        /*if ( evt.Pair_extratracks0p5mm[j] != 0 ) continue;
        if ( fabs( evt.KalmanVertexCand_z[j] ) > 15. ) continue;
        if ( 1.-fabs( evt.Pair_dphi[j] )/M_PI > 0.009 ) continue;*/

        /*if ( r == (Region)Z )
          cout << "CANDIDATE!!!" << endl;*/

        id_pair = j;
        dil_cand++;
      }
      if ( id_pair < 0 ) continue;
      if ( dil_cand != 1 ) continue;

      const double miss_mass = ( pp1+pp2-pl1-pl2 ).M();

      for ( unsigned short k = 0; k < miss_mass_bins.size()-1; ++k )
        if ( miss_mass_bins[k] < miss_mass && miss_mass < miss_mass_bins[k+1] ) {
          h_pp_mass_permmb[k].Fill( ( pp1+pp2 ).M() );
          break;
        }

      h_miss_mass[r].Fill( miss_mass );
      h_miss_mass_zoom[r].Fill( miss_mass );
      h_pp_mass[r].Fill( ( pp1+pp2 ).M() );
      const TLorentzVector pl_lead = ( pl1.Pt() > pl2.Pt() ) ? pl1 : pl2;
      const TLorentzVector pl_sublead = ( pl1.Pt() > pl2.Pt() ) ? pl2 : pl1;
      h_pt_leadmu[r].Fill( pl_lead.Pt() );
      h_pt_subleadmu[r].Fill( pl_sublead.Pt() );
      h_eta_leadmu[r].Fill( pl_lead.Eta() );
      h_eta_subleadmu[r].Fill( pl_sublead.Eta() );
      h_xangle_all[r].Fill( cond.crossing_angle );

      const TLorentzVector pll = pl1+pl2;
      h_m_pair[r].Fill( pll.M() );
      h_rap_pair[r].Fill( pll.Rapidity() );
      h_pt_pair[r].Fill( pll.Pt() );
      h_acop[r].Fill( 1.-fabs( evt.Pair_dphi[id_pair] )/M_PI );
      h_pt_pair_all[r].Fill( pll.Pt() );
      h_m_pair_all[r].Fill( pll.M() );
      h_pp_mass_all[r].Fill( ( pp1+pp2 ).M() );
      h_rap_pair_all[r].Fill( pll.Rapidity() );
      h_acop_all[r].Fill( 1.-fabs( evt.Pair_dphi[id_pair] )/M_PI );
      if ( r == (unsigned short)Region::Z && miss_mass > 1200. && miss_mass < 1350 ) {
        h_pt_pair_all[(unsigned short)Region::excess].Fill( pll.Pt() );
        h_m_pair_all[(unsigned short)Region::excess].Fill( pll.M() );
        h_pp_mass_all[(unsigned short)Region::excess].Fill( ( pp1+pp2 ).M() );
        h_rap_pair_all[(unsigned short)Region::excess].Fill( pll.Rapidity() );
        h_acop_all[(unsigned short)Region::excess].Fill( 1.-fabs( evt.Pair_dphi[id_pair] )/M_PI );
        h_xangle_all[(unsigned short)Region::excess].Fill( cond.crossing_angle );
        g_xi1xi2.SetPoint( g_xi1xi2.GetN(), xi_arm0, xi_arm1 );
      }

      //--- store a copy of the event
      //out_evt = evt;
      //out_tree.Fill();
    }
    ++file_id;
  }

  //out_tree.Write();
  //delete out_file;
  const string upper_label = Form( "CMS Preliminary 2017, L = %.2f fb^{-1}, #sqrt{s} = 13 TeV", int_lumi );

  unordered_map<string,vector<TH1D>&> plots_1d = {
    { "cand_mpair", h_m_pair }, { "cand_mpair_all", h_m_pair_all },
    { "cand_ptpair", h_pt_pair }, { "cand_ptpair_all", h_pt_pair_all },
    { "cand_pp_mass", h_pp_mass }, { "cand_pp_mass_all", h_pp_mass_all },
    { "cand_rappair", h_rap_pair }, { "cand_rappair_all", h_rap_pair_all },
    { "cand_ptlead", h_pt_leadmu }, { "cand_ptsublead", h_pt_subleadmu },
    { "cand_etalead", h_eta_leadmu }, { "cand_etasublead", h_eta_subleadmu },
    { "cand_acop", h_acop }, { "cand_acop_all", h_acop_all },
    { "cand_miss_mass", h_miss_mass }, { "cand_miss_mass_zoom", h_miss_mass_zoom },
    { "xing_angle_all", h_xangle_all },
  };
  vector<int> colours = { kBlack, kRed+1, kGreen+2, kBlue-2, kOrange };
  vector<int> markers = { 1, 20, 24, 21, 25, 22 };
  gStyle->SetOptStat( 0 );
  for ( auto& hp : plots_1d ) {
    Canvas c( hp.first.c_str(), upper_label.c_str(), true );
    c.SetLegendX1( 0.55 );
    unsigned short i = 0;
    Canvas::HistsMap hm; // typedef std::vector< std::pair<std::string,TH1*> > HistsMap;
    THStack hs;
    for ( auto& h : hp.second ) {
      if ( i > 0 )
        h.SetMarkerStyle( markers[i] );
      h.SetLineColor( colours[i] );
      h.SetLineStyle( i+1 );
      h.SetLineWidth( 3 );
      std::ostringstream os; os << (Region)i;
      /*auto hwo = (TH1D*)WithOverflow<TH1D>( &h ).Clone();
      hs.Add( hwo );
      c.AddLegendEntry( hwo, os.str().c_str(), "l" );
      hm.emplace_back( pair<string,TH1*>{ os.str(), hwo } );*/
      hs.Add( &h );
      c.AddLegendEntry( &h, os.str().c_str(), "l" );
      hm.emplace_back( pair<string,TH1*>{ os.str(), &h } );
      i++;
    }
    hs.Draw( "nostack" );
    hs.GetHistogram()->SetTitle( hp.second.begin()->GetTitle() );
    c.Prettify( hs.GetHistogram() );
    //c.RatioPlot( hm, 0., 5.2, 1., false );
    //c.RatioPlot( hm, 0., 8.4, 1., false );
    //c.RatioPlot( hm, 0., 16.8, 1., false );
    //c.RatioPlot( hm, 0., 11.8, 1., false );
    c.RatioPlot( hm, 0., 2.8, 1., false );
    if ( hp.first == "cand_mpair"
      || hp.first == "cand_mpair_all"
      || hp.first == "xing_angle_all" )
      c.GetPad( 1 )->SetLogy();
    pt_runs.Draw();
    c.Save( "pdf,png,root", WWW_PATH );
  }
  {
    Canvas c( "mpp_per_miss_mass_bin", upper_label.c_str() );
    c.SetLegendX1( 0.4 );
    c.SetLegendSizeX( 0.48 );
    THStack hs;
    for ( unsigned short k = 0; k < miss_mass_bins.size()-1; ++k ) {
      h_pp_mass_permmb[k].SetLineColor( colours[k/2] );
      h_pp_mass_permmb[k].SetLineWidth( 3 );
      h_pp_mass_permmb[k].SetLineStyle( 1+k%2 );
      hs.Add( &h_pp_mass_permmb[k] );
      //c.AddLegendEntry( &h_pp_mass_permmb[k], Form( "%g < m_{miss} < %g GeV", miss_mass_bins[k], miss_mass_bins[k+1] ) );
      c.AddLegendEntry( &h_pp_mass_permmb[k], Form( "[%g,%g]", miss_mass_bins[k], miss_mass_bins[k+1] ), "l" );
    }
    c.GetLegend()->SetNColumns( 2 );
    hs.Draw( "hist,nostack" );
    hs.GetHistogram()->SetTitle( h_pp_mass_permmb.begin()->GetTitle() );
    pt_runs.Draw();
    c.Prettify( hs.GetHistogram() );
    c.Save( "pdf,png,root", WWW_PATH );
  }
  {
    Canvas c( "cand_xi1vsxi2", upper_label.c_str() );
    g_xi1xi2.Draw( "ap" );
    g_xi1xi2.SetMarkerStyle( 24 );
    g_xi1xi2.GetHistogram()->SetTitle( ";#xi(sector 4-5);#xi(sector 5-6)" );
    g_xi1xi2.GetXaxis()->SetLimits( 0., 0.4 );
    g_xi1xi2.GetYaxis()->SetRangeUser( 0., 0.4 );
    auto pt_runs_bis = (PaveText*)pt_runs.Clone();
    pt_runs_bis->SetY1( 0.79 ); pt_runs_bis->SetY2( 0.87 );
    pt_runs_bis->Draw();
    c.Prettify( g_xi1xi2.GetHistogram() );
    c.Save( "pdf,png,root", WWW_PATH );
  }
  {
    Canvas c( "repart_xangle", upper_label.c_str() );
    c.Divide( ceil( 0.5*num_samples ), 2 );
    vector<TPie*> pies;
    for ( unsigned short i = 0; i < num_samples; ++i ) {
      c.cd( i+1 );
      unsigned short j = 0;
      pies[i] = new TPie( Form( "pie_%d", i ), Form( "Run %s", runs_list[i] ), reparts_per_xangle[i].size() );
      for ( const auto& val_vs_xangle : reparts_per_xangle[i] ) {
        pies[i]->SetEntryLabel( j, val_vs_xangle.second > 0 ? ( val_vs_xangle.first > 0. ? Form( "%g #murad", val_vs_xangle.first ) : "Other value" ) : "" );
        pies[i]->SetEntryVal( j, val_vs_xangle.second );
        pies[i]->SetEntryFillColor( j, colours[j] );
        ++j;
      }
      pies[i]->SetRadius( 0.25 );
      pies[i]->SetLabelsOffset( -0.05 );
      pies[i]->Draw( "nol <,c" );
    }
    c.Save( "pdf,png,root", WWW_PATH );
  }
}

