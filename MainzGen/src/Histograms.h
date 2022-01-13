//
//
//

using namespace std;

class Hist {

  double *data;
  int    nx, ny;
  double minx, maxx, miny, maxy;

public:
 const char *title, *labelx, *labely, *labelz, *unitx, *unity, *unitz;

  ~Hist() {delete[] data;};
  Hist(const char *title, 
       const char *labelx, const char *labely, 
       const char *unitx,  const char *unity,
       int nx, double minx, double maxx) {
    Hist::title=title; Hist::labelx=labelx; Hist::labely=labely;
    Hist::unitx=unitx; Hist::unity=unity;   Hist::minx=minx; Hist::maxx=maxx;
    Hist::nx=nx; 
    data = new double[nx];
    for(int i=0;i<nx;i++) data[i] = 0; 
    ny = 0;
  }
  Hist(const char *title, 
       const char *labelx, const char *labely, const char *labelz,
       const char *unitx,  const char *unity, const char *unitz,
       int nx, double minx, double maxx, int ny, double miny, double maxy) {
    Hist::title=title;
    Hist::labelx=labelx; Hist::labely=labely; Hist::labelz=labelz;
    Hist::unitx=unitx;   Hist::unity=unity;   Hist::unitz=unitz;
    Hist::minx=minx;     Hist::maxx=maxx;
    Hist::miny=miny;     Hist::maxy=maxy;
    Hist::nx=nx;         Hist::ny=ny; 
    data = new double[nx*ny];
    for(int i=0; i < nx * ny; i++) data[i] = 0; 
  }
  void fill2d(double x, double y, double weight=1) {
    if (x < minx || x >= maxx || y < miny || y >= maxy) return; 
    if (isnan(weight)) {return;}
    int ix = (x - minx) / (maxx-minx) * nx;
    int iy = (y - miny) / (maxy-miny) * ny;
    if (ix<0 || ix>=nx || iy<0 || iy >= ny) return;
    data[ix*ny + iy] += weight;
  }

  void fill(double x, double weight=1) {
    if (x < minx || x >= maxx) return; 
    if (isnan(weight)) return;
    int ix = (x - minx) / (maxx-minx) * nx;
    if (ix<0 || ix>=nx) { return;}//cerr << title <<" "<<x<<endl; return;}
    data[(int)((x - minx) / (maxx-minx) * nx)] += weight;
  }

  friend ostream &operator<<(ostream &out, const Hist id) {
    out << "# Histogram " << id.title << endl
	<< "# "<<endl;
    if (id.ny) {
      for (int j=0;j<id.ny;j++) {
	for (int i=0;i<id.nx;i++) 
	  out <<" "<< id.data[i*id.ny + j];
	out << endl;
      }
    } else {
      for (int i=0;i<id.nx;i++) 
	out << id.minx+(i+0.5)*(id.maxx-id.minx)/id.nx <<" "<<id.data[i]<<endl;
    }
    out << endl << endl;
    return out;
  }

  void writeGnuplotCommands(ostream &out, int i, char *fn) {
    out << "##### " <<title<<" ##############"
	<< endl << "set xlabel '" << labelx;
    if (unitx[0]) out << " [" << unitx << "]";
    out   << "'\nset ylabel '" << labely;
    if (unity[0]) out << " [" << unity << "]";
    out   << "'\nset title '" << title << "' font 'Helvetica,14'"<<"\n";
    
    if (ny)
      out <<"set style data image" <<endl
	  << "plot  ["<<minx<<":"<<maxx<<"] ["<<miny<<":"<<maxy<<"] '" 
	  << fn << "' index " << i 
	  << " using ("<<minx<<"+$1*"<<(maxx-minx)/nx
	  <<      "):("<<miny<<"+$2*"<<(maxy-miny)/ny
	  <<"):3 matrix notitle\n";    
    else 
      out <<"set style data histeps" <<endl
	  << "plot ["<<minx<<":"<<maxx<<"] '" << fn << "' index " << i 
	  << " linecolor 0 linetype 1 notitle\n";    
  }
};


