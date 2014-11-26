//this is vtkScratch/samplePlane.cxx
#include"util.h"

int main()
{
  const float a=1,b=1,c=1,d=2,cx=0,cy=0,zstd=0.25,xystd=3.0;
  const unsigned npts = 100;
  vcl_vector<vgl_homg_point_3d<float> > pts =
    samplePlanarPoints<float, vgl_homg_point_3d>(a,b,c,d,cx,cy,zstd,xystd,npts);

  scatterPoints(pts);

  drawAxes();
  vtkBoilerPlate();
  
  return EXIT_SUCCESS;
}


