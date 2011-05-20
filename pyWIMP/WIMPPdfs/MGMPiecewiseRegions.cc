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

Double_t MGMPiecewiseRegions::GetAcceptanceOverRegion(Double_t beg, Double_t end) const
{
  if (beg >= end) {
    std::cout << "Region mal-formed, beginning >= end" << std::endl; 
    return 0.0;
  }
  
  InitializeRegionIterator(beg, end);
  Double_t sum = 0.0;
  const MGMPiecewiseRegions::MGMRegion* reg;
  while ( (reg = GetNextRegion()) ) sum += reg->fEnd - reg->fBegin;
  return sum/(end - beg);

}

Bool_t MGMPiecewiseRegions::IsInAcceptedRegion(Double_t val) const
{
  MGMRegion temp(val, val);
  return ( fSetOfRegions.find( temp ) != fSetOfRegions.end() ); 
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
  if ( fCurrentIter == fSetOfRegions.end() ) return NULL;
  tempRegion = *fCurrentIter;
  if ( tempRegion.fBegin < fRequestedRegion.fBegin ) {
    tempRegion.fBegin = fRequestedRegion.fBegin;
  }
  if ( tempRegion.fEnd > fRequestedRegion.fEnd ) {
    tempRegion.fEnd = fRequestedRegion.fEnd;
  }
  
  if ( fCurrentIter != fMaxIter ) fCurrentIter++;
  else fCurrentIter = fSetOfRegions.end();
  return &tempRegion;

}

