// The code implements a "Semi-Internal Process" for the Monte Carlo
// event generator PYTHIA 8. It is not part of PYTHIA 8. However, it
// is based on PYTHIA 8. Cf. PYTHIA 8's files SigmaQCD.h, SigmaQCD.cc,
// SigmaEW.cc .

// Copyright 2021, 2022 Jens Salomon, Pablo González Durán, Torbjorn Sjostrand.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


// "main.cc" file contains the LO  calculation 
// done by Federico Deganutti.
// #include "main.cc"


namespace Pythia8 {

// QCD Pomeron base class.
class Sigma2QCDPomeron : public Sigma2Process {

public:
  // Initialize process if necessary.
  //virtual void initProc() {}

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

protected:
  // Constructor.
  Sigma2QCDPomeron(bool switchNLL=false) : sigma(), NLL(switchNLL) {}

  double sigma;
  bool NLL;

};


class Sigma2qq2qqQCDPomeron : public Sigma2QCDPomeron {

public:

  // Constructor.
  Sigma2qq2qqQCDPomeron(bool switchNLL=false) : Sigma2QCDPomeron(switchNLL) {}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q q =Pomeron=> q q";}
  // The IPROC code in HERWIG 6 is 2400.
  virtual int code() const {return 24001;}
  virtual string inFlux() const {return "qq";}
};

class Sigma2qg2qgQCDPomeron : public Sigma2QCDPomeron {

public:

  // Constructor.
  Sigma2qg2qgQCDPomeron(bool switchNLL=false) : Sigma2QCDPomeron(switchNLL) {}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q g =Pomeron=> q g";}
  virtual int code() const {return 24002;}
  virtual string inFlux() const {return "qg";}
};

class Sigma2gg2ggQCDPomeron : public Sigma2QCDPomeron {

public:

  // Constructor.
  Sigma2gg2ggQCDPomeron(bool switchNLL=false) : Sigma2QCDPomeron(switchNLL) {}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "g g =Pomeron=> g g";}
  virtual int code() const {return 24003;}
  virtual string inFlux() const {return "gg";}
};

void Sigma2QCDPomeron::sigmaKin() {
  double deltaEta = log(sH/-tH);
  double pT = sqrt(pT2);

  const double Nc = 3.;
  const double CA = Nc;
  const double CF = (Nc*Nc-1)/(2*Nc);
  const double impactFactorLO = pow(alpS*CA,2)/(Nc*Nc-1);

  double alpSbar = alpS*Nc/M_PI;

  if (NLL)
  {
    // These parameters are taken from C. Royon's HWHSNM-2.f
    const double para = 47.414, parb = -121.50, parc = 119.93, pard = 5.9833,
      pare = -17.199, parf = 0.0072066, parg = -0.29812, parh = 0.55726,
      pari = 10.385, parl = 1.5660, parm = -3.1149, parn = 1.3812;

    // Formula from C. Royon's code. It is not consistent with the
    // paper Phys.Rev.D83,034036 due to an extra factor 3/pi.
    // double zVar = (1.5/M_PI)*alpSbar*deltaEta;
    double zVar = alpSbar * deltaEta / 2;
    double zVar2 = pow(zVar, 2);
    double zVar3 = pow(zVar, 3);

    // Compute amplitude squared.
    double asq = 0;
    asq = para + parf*pT + parl*sqrt(pT);
    asq += (parb+parg*pT+parm*sqrt(pT))*zVar;
    asq += (parc + parh*pT)*zVar2;
    asq += (pari + parn*sqrt(pT))*zVar3;
    asq += exp(pard + pare*zVar);
    asq = asq/(4.*M_PI) * pow(alpS, 4.);
    //asq = asq*sH2/tH2*(16*M_PI);

    asq *= pow(CA, 4);

    //sigma = asq/(sH2*16*M_PI);
    sigma = asq/tH2;
  }
  else if (!NLL)
  {
		// See line 22 if LO is desired.
    // 1/tH2 is the qFactor.
    // sigma = 1/M_PI*pow(impactFactorLO,2)
    //   *1/tH2*pow(BFKL::LO::gluonGreenFunction(deltaEta, alpSbar), 2);
  }

  // Adapt colour factor to incoming particles other than two gluons.
  if ( inFlux() == "qg" )
    sigma *= pow(CF/CA, 2);
  else if ( inFlux() == "qq" )
    sigma *= pow(CF/CA, 4);
}

void Sigma2qq2qqQCDPomeron::setIdColAcol() {
  // Flavours are trivial.
  setId( id1, id2, id1, id2);

  // Colour flow topologies. Swap when antiquarks.
  // This is implemented on the basis of Sigma2ff2fftgmZ::setIdColAcol.
  if (id1*id2 > 0)
    setColAcol( 1, 0, 2, 0, 1, 0, 2, 0);
  else
    setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  if ( id1 < 0 )
    swapColAcol();
}

void Sigma2qg2qgQCDPomeron::setIdColAcol() {
  // Flavours are trivial.
  setId( id1, id2, id1, id2);

  setColAcol( 1, 0, 2, 3, 1, 0, 2, 3);
  if (id1 == 21) swapCol1234();
  if (id1 < 0 || id2 < 0) swapColAcol();
}

void Sigma2gg2ggQCDPomeron::setIdColAcol() {
  // Flavours are trivial.
  setId( id1, id2, id1, id2);

  setColAcol( 1, 2, 3, 4, 1, 2, 3, 4);
}


} // end namespace Pythia 8
