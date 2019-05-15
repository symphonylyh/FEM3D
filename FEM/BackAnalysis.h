/**
 * @file BackAnalysis.h
 * Derived class from Analysis for the back-analysis scheme.
 *
 * @author Haohang Huang
 * @date September 16, 2018
 */

#ifndef BackAnalysis_h
#define BackAnalysis_h

#include "Analysis.h"

/* Derived class for solving back analysis problems.
 */
class BackAnalysis : public Analysis
{
  public:
    /* See the documentation of base class Analysis.
     */
    BackAnalysis(Mesh & meshInfo); // ctor cannot be inherited, should explicitly call base class's ctor in derived class's ctor
    ~BackAnalysis();
    void solve();

};

#endif /* BackAnalysis_h */
