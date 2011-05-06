#include "MGMPiecewiseRegions.hh" 
#include <iostream>

ClassImp( MGMPiecewiseRegions )

Bool_t MGMPiecewiseRegions::InsertNewRegion(Double_t beginning, Double_t end)
{
  if (beginning >= end) {
    std::cout << "Region mal-formed, beginning >= end" << std::endl; 
    return false;
  }
  // Try to insert it 
  std::pair<RegionSet::iterator, bool> ret; 
  ret = fSetOfRegions.insert(MGMRegion(beginning, end));
  if (!ret.second) {
    std::cout << "Attempt to insert overlapping regions." << std::endl; 
    return false;
  }
  return true;
} 

Bool_t MGMPiecewiseRegions::IsInAcceptedRegion(Double_t val) const
{
  MGMRegion temp(val, val);
  if ( fSetOfRegions.find( temp ) != fSetOfRegions.end() ) return true; 
  return false;
}

void MGMPiecewiseRegions::InitializeRegionIterator(Double_t beginning, Double_t end) const
{
  if (end < beginning) {
    std::cout << "Entry mal-formed." << std::endl; 
    return;
  }
  fRequestedRegion.fBegin = beginning; 
  fRequestedRegion.fEnd   = end; 

  // Now find the region iterator
  // First check to see if we are within a range
  fCurrentIter = fSetOfRegions.lower_bound(MGMRegion(beginning, beginning));
  fMaxIter     = fSetOfRegions.lower_bound(MGMRegion(end, end));
}

const MGMPiecewiseRegions::MGMRegion* 
  MGMPiecewiseRegions::GetNextRegion() const
{
  static MGMRegion tempRegion;
  const MGMRegion* tempRet = NULL;
  if ( fCurrentIter == fSetOfRegions.end() ) return NULL;
  tempRet = &(*fCurrentIter);
  if ( fCurrentIter != fMaxIter ) {
    if ( tempRet->fBegin < fRequestedRegion.fBegin ) {
      tempRegion.fBegin = fRequestedRegion.fBegin;
      tempRegion.fEnd = tempRet->fEnd;
      tempRet = &tempRegion;
    }
    fCurrentIter++;
    return tempRet;
  }
  // This means the current Iter == fMaxIter
  if ( tempRet->fEnd > fRequestedRegion.fEnd ) {
    tempRegion.fBegin = tempRet->fBegin;
    tempRegion.fEnd = fRequestedRegion.fEnd;
    tempRet = &tempRegion;
  }
  fCurrentIter = fSetOfRegions.end();
  return tempRet;

}

