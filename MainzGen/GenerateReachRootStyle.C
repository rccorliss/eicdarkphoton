

const float alpha=1/137.036;
const float pi=3.14159;
const float epsilon0=1.;






//need:
//function to load branching ratio from file
class BranchingRatio{
public:
  std::vector<float> m,ratio,integral;
  float InterpolatedValue(float x, std::vector<float> source){
    for (int i=0;i<m.size()-1;i++){
      if(x>=m[i] && x<m[i+1]) {
	//interpolate:
	double fractionalpos=(x-m[i])/(m[i+1]-m[i]);
	return source[i]*(1-fractionalpos)+source[i+1]*fractionalpos;
      }
    }
    return source[m.size()-1];//if we're too high, return the last value;  Shouldn't happen.
  };

  public:
  BranchingRatio(const char* filename){
    std::ifstream in(filename);
    char dummy;
    float x,y,xold,yold,intold;
    m.push_back(0);      ratio.push_back(1); integral.push_back(0); //set a zero as a lower bound
    x=0;y=1;
    while (!in.eof()) {
      xold=x;yold=y;intold=integral[integral.size()-1];
      in>> x>>dummy>>y;
      m.push_back(x);    ratio.push_back(y); integral.push_back(intold+(yold+y)/2*(x-xold));
      printf("i=%d: x=%f,y=%f,int=%f\n",(int)m.size()-1,x,y,integral[integral.size()-1]);
    }
    xold=x;yold=y;intold=integral[integral.size()-1];
    x=100;y=100;
    m.push_back(100); ratio.push_back(y); integral.push_back(intold+(yold+y)/2*(x-xold)); //assume it stays flat above our highest point
    in.close();
    return;
  };
  float At(float x){ return InterpolatedValue(x,ratio);};
  float Average(float xlow,float xhigh){
    return (InterpolatedValue(xhigh,integral)-InterpolatedValue(xlow,integral))/(xhigh-xlow);};
};

//function to fit a gaussian
float GetGaussianWidth(TH2F*hist,int xbin, float *err);

//function to optimize the window
float GetOptimalWindow(float sigma);
//function to read a 2d hist out of the .dat files
void LoadFromFile2D(TH2F*hist,const char* filename, const char* histlabel);
//function to read a 1d hist out of the .dat files
void LoadFromFile1D(TH1F*hist,const char* filename, const char* histlabel);


void GenerateReachRootStyle(){


  BranchingRatio *br=new BranchingRatio("BranchingRatioExtracted.csv");
  float GeV=1.0;
  int nBinsX=100,nBinsY=100;
  float minTrueMass=0, maxTrueMass=10*GeV;
  float minSmearMassDiff=-0.1*GeV, maxSmearMassDiff=0.1*GeV;
  float binstep=(maxTrueMass-minTrueMass)/(1.*nBinsX);
  printf("binstep=%f\n",binstep);

  
  TH2F *hTrueMassVsReco=new TH2F("hTrueMassVsReco","True Mass vs Reco Mass (relative to true);true (GeV);reco-true (GeV)",nBinsX,minTrueMass,maxTrueMass,nBinsY,minSmearMassDiff,maxSmearMassDiff);
  LoadFromFile2D(hTrueMassVsReco,"Signal_EIC1x1.dat","# Histogram Gen and smeared Aprime mass");

  TH1F *hBranching=new TH1F("hBranching","Average Branching Ratio per bin;true mass[GeV]",nBinsX,minTrueMass,maxTrueMass);

  for (int i=0;i<nBinsX;i++){
    float mlow=minTrueMass+binstep*i;
    float mhigh=mlow+binstep;
    float inthigh=br->InterpolatedValue(mhigh,br->integral);
    float intlow=br->InterpolatedValue(mlow,br->integral);
    //printf("inthigh=%f (m=%f),intlow=%f (m=%f),int=%f, ave=%f\n",inthigh,mhigh,intlow,mlow,inthigh-intlow, (inthigh-intlow)/(mhigh-mlow));
    hBranching->Fill(minTrueMass+brstep*(i+0.5),br->Average(mlow,mhigh));
  }
  TH1F *hBgMassReco=new TH1F("hBgMassReco","Reconstructed Background Mass;reco (GeV)",nBinsX,0.0,10.);
  LoadFromFile1D(hBgMassReco,"Integration_EIC1x1.dat","# Histogram smeared Aprime mass");

  TCanvas *c=new TCanvas("c","",1000,500);
  c->Divide(3,1);
  c->cd(1);
  hTrueMassVsReco->Draw("colz");
  c->cd(2);
  hBgMassReco->Draw();
  c->cd(3);
  hBranching->Draw("hist");
  
  //we now have signal, background, and the tools for the normalization factor
  for (int i=0;i<nBinsX;i++){
    float binlow=minTrueMass+binstep*i;
    float binmid=binlow+binstep*0.5;
    float binhigh=binlow+binstep;

    int bin=hTrueMassVsReco->GetXaxis()->FindBin(binmid)
    TH1F *hSlice=hist->ProjectionY("_yslice",bin,bin); //momentum distribution around this nominal mass;

    //find the gaussian width of the smeared signal, use that to optimize the window:
    float sigmaErr;
    float sigma=GetGaussianWidth(hTrueMassVsReco, &sigmaErr);
    float binwidth=GetOptimalWindow(sigma);

    

    //the Mainz scale factor, not included in the IntegrateQED calc because it depends on our chosen binning, not just the generator:
    float N=br->Average(binlow,binhigh);//get the average branching ratio.
    float signal_scale_factor=(3*pi)/(2*N) * (epsilon0*epsilon0)/(alpha) * (binmid)/(binstep);

    //sum and rescale the signal in our window:
    int sigbinmin=hSlice->GetXaxis()->FindBin(-binwidth);
    int sigbinmax=hSlice->GetXaxis()->FindBin(binwidth);
    float totalsignal=signal_scale_factor*hSlice->Integrate(sigbinmin,sigbinmax);

    
    float totalbackground=hBgMassReco

    
      }
  
  return;
}

float GetGaussianWidth(TH1F*hist, float *err){
  //maybe eventually this would be a fit?
  //    TF1 *fGaus=new TF1("fGaus","gaus(0)",0,200);

  *err=hist->GetStdDevError();
  return hist->GetStdDev();
}

  
float GetOptimalWindow(float sigma){
  //ratio=sig/sqrt(bg)
  //we assume bg is flat, so ratio =K*sig_in_window/sqrt(window width);
  //this is a partial integral of a gaussian = an erf function
  //and has no analytical optimal solution, so we scan:
  static TF1 r("analyticalratio","TMath::Erf(x/[0])/sqrt(x)",0.05,10);
  r.SetParameter(0,sigma);
  float bestw=r.GetMaximumX(0.005,1);
  new TCanvas("junk","junk canvas",400,600);
  r.Draw();
  TLine j;
  j.DrawLine(bestw,0,bestw,100);
  return bestw;
}  

  
void LoadFromFile2D(TH2F*hist,const char* filename, const char* histlabel){

  std::ifstream in(filename);
  //seek the right line:
   while (!in.eof()) {
     std::string token;
    std::getline(in,token);
    //printf("token=%s\n",token.c_str());
    if (!strcmp(token.c_str(), histlabel)){
     std::getline(in,token); //when we find our thing, we need one more line
      break;
    }
   }
   if (in.eof()){printf("broken file!\n"); exit(1); return;}

   TAxis *yax=hist->GetYaxis();
   float ymin=yax->GetXmin();
   float ymax=yax->GetXmax();
   int ybins=yax->GetNbins();
   float ystep=(ymax-ymin)/ybins;
   TAxis *xax=hist->GetXaxis();
   float xmin=xax->GetXmin();
   float xmax=xax->GetXmax();
   int xbins=xax->GetNbins();
   float xstep=(xmax-xmin)/xbins;

   float data;
   for (int j=0;j<ybins;j++) {
     for (int i=0;i<xbins;i++) {
       in >> data;
       if (data>0){
	 hist->Fill(xmin+(i+0.5)*xstep,ymin+(j+0.5)*ystep,data);
       }
     }
   }
   in.close();
   return;
}


void LoadFromFile1D(TH1F*hist,const char* filename, const char* histlabel){

  std::ifstream in(filename);
  //seek the right line:
  while (!in.eof()) {
    std::string token;
    std::getline(in,token);
    //printf("token=%s\n",token.c_str());
    if (!strcmp(token.c_str(), histlabel)){
      std::getline(in,token); //when we find our thing, we need one more line
      printf("found the leader line\n");
      break;
    }
  }
  if (in.eof()){printf("broken file!\n"); exit(1); return;}

  TAxis *xax=hist->GetXaxis();
  float xmin=xax->GetXmin();
  float xmax=xax->GetXmax();
  int xbins=xax->GetNbins();
  float xstep=(xmax-xmin)/xbins;

  //printf("expecting %d bins\n",xbins);
  float x,data;
  for (int i=0;i<xbins;i++) {
    in >> x >> data;
    //printf("found bin at %1.2f, data=%E\n",x,data);
    hist->Fill(x,data);

  }
  in.close();
  return;
}
