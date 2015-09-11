/** file: flashggPuppiProducer.cc ** C/C++ **
 *
 * Author: yhaddad <yhaddad@cern.ch>
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <memory> 
#include <algorithm>    // std::search
#include <iomanip>    // std::setw()

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/PtEtaPhiMass.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "CommonTools/PileupAlgos/interface/PuppiContainer.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/VertexCandidateMap.h"

#include "fastjet/PseudoJet.hh"

using namespace edm;
using namespace std;
using namespace flashgg;

//namespace flashgg {
class flashggPuppiProducer : public EDProducer
{
    
public:
  flashggPuppiProducer( const ParameterSet & );
  virtual ~flashggPuppiProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  // --- puppi stuff
  typedef math::XYZTLorentzVector           LorentzVector;
  typedef std::vector<LorentzVector>        LorentzVectorCollection;
  typedef reco::VertexCollection            VertexCollection;
  typedef edm::View<reco::Candidate >       CandidateView;
  typedef std::vector<reco::PFCandidate >   PFInputCollection;
  typedef std::vector<reco::PFCandidate >   PFOutputCollection;
  typedef edm::View<reco::PFCandidate>      PFView;

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
    
  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  virtual void endRun(edm::Run&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  
  int puppi_std_index(edm::Ptr<pat::PackedCandidate> cand, double dzcut, bool print);
  
  EDGetTokenT<View<reco::Vertex> >               vertexToken_;
  EDGetTokenT<View<flashgg::DiPhotonCandidate> > diPhotonsToken_;
  EDGetTokenT<View<pat::PackedCandidate> >       pfcandidateToken_;
  EDGetTokenT< VertexCandidateMap >              vertexCandidateMapToken_;
    

  // --- puppi stuff
  std::string     fPuppiName;
  std::string     fPFName;
  std::string     fPVName;
  bool            fPuppiDiagnostics;
  bool            fPuppiForLeptons;
  bool            fUseDZ;
  float           fDZCut;
  std::unique_ptr< PuppiContainer>     fPuppiContainer;
  std::vector< RecoObj>                fRecoObjCollection;
  std::auto_ptr< PFOutputCollection >  fPuppiCandidates;
  
  // --- muti vertex stuff
  unsigned     indexVtx_;
  bool         debug_;
  unsigned int eventNb;
};

int flashggPuppiProducer::puppi_std_index(edm::Ptr<pat::PackedCandidate> cand, double dzcut, bool print)
{
  if(cand->charge() == 0) return 0;
  if(fabs(cand->charge()) > 0){
    if(cand->fromPV() == 0) return 2;
    if(cand->fromPV() == 3) return 1;
    if(cand->fromPV() == 1 || cand->fromPV() == 2){
      if (print) std::cout << "stdpuppi::dz()" << cand->dz() << " : " << dzcut << std::endl; 
      if(fabs(cand->dz()) < dzcut ) return 1;
      if(fabs(cand->dz()) > dzcut ) return 2;
    }
  }
  return 0;
}

flashggPuppiProducer::flashggPuppiProducer( const ParameterSet &iConfig ) :
  vertexToken_     (consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag>("vertexName"))),
  diPhotonsToken_  (consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag>("diPhotonTag"))),
  pfcandidateToken_(consumes<View<pat::PackedCandidate> >( iConfig.getParameter<InputTag>("candName"))),
  vertexCandidateMapToken_(consumes<VertexCandidateMap>( iConfig.getParameter<InputTag>("VertexCandidateMapTag")))
{
  
  indexVtx_         = iConfig.getUntrackedParameter<unsigned> ( "vertexIndex",  0 );
  debug_            = iConfig.getUntrackedParameter<bool>     ( "debug"      ,  false );
  
  fPuppiDiagnostics = iConfig.getParameter<bool>  ("puppiDiagnostics");
  fPuppiForLeptons  = iConfig.getParameter<bool>  ("puppiForLeptons");
  fUseDZ            = iConfig.getParameter<bool>  ("UseDeltaZCut");
  fDZCut            = iConfig.getParameter<double>("DeltaZCut");
    
  fPuppiContainer = std::unique_ptr<PuppiContainer> ( new PuppiContainer(iConfig) );
  
  produces<edm::ValueMap<float> > ();
  produces<edm::ValueMap<LorentzVector> > ();
  produces< edm::ValueMap<reco::CandidatePtr> >(); 
  produces<PFOutputCollection>();
    
  if (fPuppiDiagnostics){
    produces<double> ("PuppiNAlgos");
    produces<std::vector<double>> ("PuppiRawAlphas");
    produces<std::vector<double>> ("PuppiAlphas");
    produces<std::vector<double>> ("PuppiAlphasMed");
    produces<std::vector<double>> ("PuppiAlphasRms");
  }
  
  
  if (debug_) std::cout << "  :: process puppi for vertex :: "<< indexVtx_ << std::endl;
  eventNb = 0 ;
}

flashggPuppiProducer::~flashggPuppiProducer(){;}

void
flashggPuppiProducer::produce( Event &evt , const EventSetup & )
{
  
  if (debug_) std::cout << "------  event :: " << eventNb
			<< "  with smart index :: "<< indexVtx_ << std::endl;
  // --- diphotons
  Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
  evt.getByToken( diPhotonsToken_, diPhotons );
  // const PtrVector<flashgg::DiPhotonCandidate>& diPhotonPointers= diPhotons->ptrVector();
  
  // --- Primary Vertices in case no diphoton found
  Handle<View<reco::Vertex> > primaryVertices;
  evt.getByToken( vertexToken_, primaryVertices );
  //const PtrVector<reco::Vertex>& pvPtrs = primaryVertices->ptrVector();
  int npv  = 0;
  for (unsigned int i=0; i < primaryVertices->size(); i++) {
    if (!primaryVertices->ptrAt(i)->isFake()  &&
	primaryVertices->ptrAt(i)->ndof()>=4  &&
	fabs(primaryVertices->ptrAt(i)->z())<=24)
      npv++;
  }
  
  if( debug_ ) std::cout << setw( 13 ) << " (n good pv) = " << npv << std::endl;
  // --- packed candicates
  Handle<View<pat::PackedCandidate> > pfCandidates;
  evt.getByToken( pfcandidateToken_, pfCandidates );
  // const PtrVector<pat::PackedCandidate>& pfPtrs = pfCandidates->ptrVector();
  
  // --- flashgg vertex map
  Handle<VertexCandidateMap> vtxmap;
  evt.getByToken( vertexCandidateMapToken_, vtxmap );
    
  // ---  Non-const copy so we can sort by candidates
  //      (which may be associated to more than one vertex)
  VertexCandidateMap vtxmapbycand( *vtxmap );
  std::stable_sort( vtxmapbycand.begin(), vtxmapbycand.end(), flashgg::compare_by_cand() );
    
  edm::Ptr<reco::Vertex> choosenVertex;
    
  // force the use of the MiniAOD default vertex
  int  diphoton_index  = -1;
  //int  counterFindVtx  =  0;
  bool keep_event      = false;
  
    
  if( indexVtx_ == 0 ) {
    choosenVertex =  primaryVertices->ptrAt( 0 );
    keep_event = true;
  }
  
  // Use the other vertices :
  // we are interested only by the di-photon candidates verticies
  // then retain only the vertex of the di-photon that matters
  else {
    for( unsigned int diPhoLoop = 0; diPhoLoop < diPhotons->size() ; diPhoLoop++ ) {
      // we only have a problem if the mutliple diphotons haev different vertices...
      if(debug_) std::cout << "\t *diphoton["<< diPhoLoop
			   << "] --> "       << diPhotons->ptrAt( diPhoLoop )->jetCollectionIndex()
			   << std::endl;
      if( diPhotons->ptrAt( diPhoLoop )->jetCollectionIndex() == indexVtx_ ) {
	choosenVertex  = diPhotons->ptrAt( diPhoLoop )->vtx();
	diphoton_index = diPhoLoop;
	keep_event = true;
	break;
      }
    }
  }
  
  if( debug_ ) std::cout << std::setw( 13 ) << " nbr(diphoton) = " << diPhotons->size()
			 << std::setw( 13 ) << " nbr(pfcandid) = " << pfCandidates->size()
			 << std::setw( 13 ) << " diphoton inx  = " << diphoton_index << std::endl;
  if( debug_ ) std::cout << std::setw( 13 ) << " (keep_event)  = " << keep_event << std::endl;
  
  
  // -----------------------------------
  // fill the fRecoObject for puppi
  // -----------------------------------
  fRecoObjCollection.clear();
  if (debug_) std::cout <<" ------------ flashgg puppi  --------------"  << std::endl;
  if (keep_event){
    for( unsigned int pfCandLoop = 0 ; pfCandLoop < pfCandidates->size() ; pfCandLoop++ ) {
      edm::Ptr<pat::PackedCandidate> cand = pfCandidates->ptrAt( pfCandLoop );
      
      RecoObj pReco;
      pReco.pt       = cand->pt();
      pReco.eta      = cand->eta();
      pReco.phi      = cand->phi();
      pReco.m        = cand->mass();
      pReco.rapidity = cand->rapidity();
      pReco.charge   = cand->charge();
      pReco.vtxId    = -9999;

      edm::Ptr<reco::Vertex> closestVertex;
      
      auto mapRange  = std::equal_range( vtxmapbycand.begin(),
					 vtxmapbycand.end(), cand,
					 flashgg::compare_with_cand() );
      
      // setup the default values
      pReco.id = 2; // 0 previously
      pReco.dZ = cand->dz  (choosenVertex->position());
      pReco.d0 = cand->dxy (choosenVertex->position());
      //unsigned int vtx_found = 0;
      int nvtx_found = 0;
      for( auto pair_iter = mapRange.first ; pair_iter != mapRange.second ; pair_iter++ ) {
	edm::Ptr<reco::Vertex> currentVertex       = pair_iter->first;
	if( choosenVertex == currentVertex ) {
	  // the vertex atteched is the one choosen one
	  // this is equivalent to fromPV() == 3
	  closestVertex = currentVertex;
	  pReco.id = 1;// new 
	}
	nvtx_found ++;
      }
      
      //if (closestVertex.isNull()){
      //	pReco.id = 2;
      //	// In case we pass the dz cut 
      //	// The candidate should be kept
      //	///if (fabs(pReco.dZ) < fDZCut) pReco.id = 1;
      //}else{
      //	// in the case the the candidate is not attached
      //	// no vertex, this is equivalent to fromPV()=0
      //	pReco.id = 1;
      //}
      
      // nautrals
      if (fabs(pReco.charge) == 0){
	pReco.id = 0;
	pReco.dZ = 0;
	pReco.d0 = 0;
      }
      
      fRecoObjCollection.push_back(pReco);
      int puppi_index = flashggPuppiProducer::puppi_std_index(cand, fDZCut,false);
      
      if(debug_ && indexVtx_ == 0 && fabs(cand->charge()) > 0 && puppi_index != pReco.id ){
      	std::cout << "\033[1;31m";
      	std::cout << std::setw(8) <<"cand= "       << std::setw(4)  << pfCandLoop
      		  << std::setw(8) <<" fromPV= "    << std::setw(2)  << cand->fromPV()
      		  << std::setw(8) <<" std_id= "    << std::setw(2)  << puppi_index
      		  << std::setw(8) <<" new_id= "    << std::setw(2)  << pReco.id
      		  << std::setw(8) <<" isvtx= "     << std::setw(2)  << !(closestVertex.isNull())
      		  << std::setw(8) <<" dz()= "      << std::setw(11) << cand->dz()
      		  << std::setw(8) <<" dz(vtx)= "   << std::setw(11) << pReco.dZ
      		  << std::endl;
      	std::cout << "\033[0m";
      }else if (debug_ && indexVtx_ == 0 && fabs(cand->charge()) > 0 && puppi_index == pReco.id){
      	std::cout << "\033[0;34m";
      	std::cout << std::setw(8) <<"cand= "       << std::setw(4)  << pfCandLoop
      		  << std::setw(8) <<" fromPV= "    << std::setw(2)  << cand->fromPV()
      		  << std::setw(8) <<" std_id= "    << std::setw(2)  << puppi_index
      		  << std::setw(8) <<" new_id= "    << std::setw(2)  << pReco.id
      		  << std::setw(8) <<" isvtx= "     << std::setw(2)  << !(closestVertex.isNull())
      		  << std::setw(8) <<" dz()= "      << std::setw(12) << cand->dz()
      		  << std::setw(8) <<" dz(vtx)= "   << std::setw(12) << pReco.dZ
      		  << std::endl;
      	std::cout << "\033[0m";
      }
    }
  }
  // -----------------------------------
  // call puppi algorithm
  // -----------------------------------
  
  std::auto_ptr<edm::ValueMap<reco::CandidatePtr> > pfMap_p(new edm::ValueMap<reco::CandidatePtr>());
  // for diagnostic
  std::auto_ptr<std::vector<double> > theAlphas   (new std::vector<double>());
  std::auto_ptr<std::vector<double> > theAlphasMed(new std::vector<double>());
  std::auto_ptr<std::vector<double> > theAlphasRms(new std::vector<double>());
  std::auto_ptr<std::vector<double> > alphas      (new std::vector<double>());
  std::auto_ptr<double>               nalgos      (new double());
  
  std::auto_ptr<edm::ValueMap<float> >         lPupOut (new edm::ValueMap<float>());
  std::auto_ptr<edm::ValueMap<LorentzVector> > p4PupOut(new edm::ValueMap<LorentzVector>());
  std::vector<reco::CandidatePtr>              values  (pfCandidates->size());
  
  fPuppiCandidates.reset( new PFOutputCollection );
  
  if (keep_event){
    fPuppiContainer->initialize(fRecoObjCollection);
    fPuppiContainer->setNPV( npv );
    
    //Compute the weights
    const std::vector<double> lWeights = fPuppiContainer->puppiWeights();
    //Fill it into the event
    
    edm::ValueMap<float>::Filler  lPupFiller(*lPupOut);
    lPupFiller.insert(pfCandidates,lWeights.begin(),lWeights.end());
    lPupFiller.fill();
    
    // This is a dummy to access the "translate" method which is a
    // non-static member function even though it doesn't need to be.
    // Will fix in the future.
    static const reco::PFCandidate dummySinceTranslateIsNotStatic;
    
    //Fill a new PF Candidate Collection and write out the ValueMap of the new p4s.
    // Since the size of the ValueMap must be equal to the input collection, we need
    // to search the "puppi" particles to find a match for each input. If none is found,
    // the input is set to have a four-vector of 0,0,0,0
    const std::vector<fastjet::PseudoJet> lCandidates = fPuppiContainer->puppiParticles();
    LorentzVectorCollection puppiP4s;
    
    for ( auto i0 = pfCandidates->begin(),
	    i0begin = pfCandidates->begin(),
	    i0end = pfCandidates->end(); i0 != i0end; ++i0 ) {
      
      auto id = dummySinceTranslateIsNotStatic.translatePdgIdToType(i0->pdgId());
      const reco::PFCandidate *pPF = dynamic_cast<const reco::PFCandidate*>(&(*i0));
      reco::PFCandidate        pCand( pPF ? *pPF : reco::PFCandidate(i0->charge(), i0->p4(), id) );
      LorentzVector pVec = i0->p4();
      int val  = i0 - i0begin;
      
      // Find the Puppi particle matched to the input collection using the "user_index" of the object.
      auto puppiMatched = find_if( lCandidates.begin(),
				   lCandidates.end(),
				   [&val]( fastjet::PseudoJet const & i )
				   { return i.user_index() == val; } );
      
      if ( puppiMatched != lCandidates.end() ) {
	pVec.SetPxPyPzE(puppiMatched->px(),puppiMatched->py(),puppiMatched->pz(),puppiMatched->E());
	fPuppiCandidates->push_back(pCand);
	//if(debug_) std::cout << std::setw(12) <<"cand= "      << val
	//		     << std::setw(12) <<"charge= "    << pCand.charge()
	//		     << std::setw(12) <<"pt= "        << pCand.pt()
	//		     << std::setw(12) <<"eta= "       << pCand.eta()
	//		     << std::setw(12) <<"phi= "       << pCand.phi()
	//		     << std::setw(12) <<"weight= "    << lWeights[val]
	//		     << std::endl;
      } else {
	pVec.SetPxPyPzE( 0, 0, 0, 0);
      }
      
      pCand.setP4(pVec);
      puppiP4s.push_back( pVec );
      //fPuppiCandidates->push_back(pCand);
    }
    
    //Compute the modified p4s
    edm::ValueMap<LorentzVector>::Filler  p4PupFiller(*p4PupOut);
    p4PupFiller.insert(pfCandidates,puppiP4s.begin(), puppiP4s.end() );
    p4PupFiller.fill();
    
    edm::ValueMap<reco::CandidatePtr>::Filler filler(*pfMap_p);
    filler.insert(pfCandidates, values.begin(), values.end());
    filler.fill();
    if (fPuppiDiagnostics){
      //theAlphas    = (fPuppiContainer->puppiAlphas());
      //theAlphasMed = (fPuppiContainer->puppiAlphasMed());
      //theAlphasRms = (fPuppiContainer->puppiAlphasRMS());
      //alphas       = (fPuppiContainer->puppiRawAlphas());
      //nalgos = fPuppiContainer->puppiNAlgos();
    }
  }
  std::cout << "flashgg::puppiCand == " << fPuppiCandidates->size()
	    << "         inputCand == " << fRecoObjCollection.size()
	    << std::endl;
  
  evt.put(lPupOut);
  evt.put(p4PupOut);
  edm::OrphanHandle<reco::PFCandidateCollection> oh = evt.put( fPuppiCandidates );
  if (keep_event){
    for(unsigned int ic=0, nc = oh->size(); ic < nc; ++ic) {
      reco::CandidatePtr pkref( oh, ic );
      values[ic] = pkref;
    }
  }
  
  evt.put(pfMap_p);
  //////////////////////////////////////////////
  if (fPuppiDiagnostics){
    // all the different alphas per particle
    // THE alpha per particle
    evt.put(alphas      ,"PuppiRawAlphas");
    evt.put(nalgos      ,"PuppiNAlgos");
    evt.put(theAlphas   ,"PuppiAlphas");
    evt.put(theAlphasMed,"PuppiAlphasMed");
    evt.put(theAlphasRms,"PuppiAlphasRms");
  }
  // -----------------------------------
  // save on the event
  // -----------------------------------
  
  eventNb++;
}

  

  

  
// ------------------------------------------------------------------------------------------
void flashggPuppiProducer::beginJob() {
}
// ------------------------------------------------------------------------------------------
void flashggPuppiProducer::endJob() {
}
// ------------------------------------------------------------------------------------------
void flashggPuppiProducer::beginRun(edm::Run&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void flashggPuppiProducer::endRun(edm::Run&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void flashggPuppiProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void flashggPuppiProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void flashggPuppiProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE( flashggPuppiProducer );
//}// end flashgg::namespace

