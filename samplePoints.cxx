//this is vtkScratch/samplePoints.cxx

#include"util.h"

int main()
{
  // vcl_vector<vgl_homg_point_3d<float> > pts =
  //   samplePlanarPoints<float, vgl_homg_point_3d>();

  // scatterPoints(pts);

  scatterPoints(sampleCircle<float,vgl_homg_point_3d>(25));

  drawAxes();
  
  vtkBoilerPlate();

  return EXIT_SUCCESS;
}
