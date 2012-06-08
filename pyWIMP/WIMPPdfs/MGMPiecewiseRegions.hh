#ifndef _MGMPiecewiseRegions_hh_
#define _MGMPiecewiseRegions_hh_
#include "TObject.h"
#include <set>
 
class MGMPiecewiseRegions : public TObject 
{
  public:
    MGMPiecewiseRegions() {}
    ~MGMPiecewiseRegions() {}

    struct MGMRegion {
        Double_t fBegin;
        Double_t fEnd;
        MGMRegion(Double_t abeg = 0, Double_t anend = 0): fBegin(abeg), fEnd(anend) {}
    };
 
    Bool_t InsertNewRegion(Double_t beginning, Double_t end);
    inline void ClearRegions() { fSetOfRegions.clear(); }

    bool IsInAcceptedRegion(Double_t val) const;

    // Use the following for looping over a set of regions.
    void InitializeRegionIterator(Double_t beginning, Double_t end) const;
    const MGMRegion* GetNextRegion() const;

    // Returns how much is accepeted (integral) over a region)
    Double_t GetAcceptanceOverRegion(Double_t beginning, Double_t end) const;

  protected:
   
    struct MGMRegionCompare {
        // This function allows us to make sure we have no overlapping
        // regions
        bool operator() (const MGMRegion& lhs, const MGMRegion& rhs) const
        { return lhs.fEnd < rhs.fBegin; }
    };

    typedef std::set<MGMRegion, MGMRegionCompare> RegionSet;
    RegionSet fSetOfRegions; 
    mutable RegionSet::iterator fCurrentIter; //!
    mutable RegionSet::iterator fMaxIter;     //!
    mutable MGMRegion fRequestedRegion;       //!
  
    ClassDef(MGMPiecewiseRegions, 1)

};

#endif
