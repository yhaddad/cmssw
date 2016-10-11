#ifndef Z_GENERATOR_LINE_SHAPE
#define Z_GENERATOR_LINE_SHAPE

#include "Riostream.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "TH1F.h"
#include "TFile.h"

class ZGeneratorLineShape : public RooAbsPdf {

public:
  ZGeneratorLineShape() {} ; 
  ZGeneratorLineShape(const char *name, const char *title,
		      RooAbsReal& _m,
                      const char* genfile=0, const char* histoName= "Mass"
		      );

  ZGeneratorLineShape(const ZGeneratorLineShape& other, const char* name);
  inline virtual TObject* clone(const char* newname) const { return new ZGeneratorLineShape(*this,newname);}
  inline ~ZGeneratorLineShape(){};
  ClassDef(ZGeneratorLineShape,2)
    Double_t evaluate() const;  
 protected:
  RooRealProxy m ;
  RooDataHist* dataHist;
};

#endif
