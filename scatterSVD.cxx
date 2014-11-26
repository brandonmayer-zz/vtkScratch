//this is /vtkScratch/scatterSVD.cxx
#include"util.h"
int main()
{
  const float a=1,b=1,c=1,d=2,pcx=0,pcy=0,zstd=0.25,xystd=3.0;
  const unsigned npts = 100;
  vcl_vector<vgl_homg_point_3d<float> > pts =
    samplePlanarPoints<float, vgl_homg_point_3d>(a,b,c,d,pcx,pcy,zstd,xystd,npts);

  scatterPoints(pts);

  vnl_matrix<float> scatterMatrix;
  vgl_norm_trans_3d<float> norm;

  getScatterMatrix(pts, scatterMatrix, norm);

  vnl_vector_fixed<float,3> centerOfMass = norm.get_translation_vector();

  vcl_cout << "center of mass: " << centerOfMass << vcl_endl;

  float cm[3];
  center_of_mass(pts,cm);
  vcl_cout << "cm: " << cm << vcl_endl;
  
  vcl_cout << scatterMatrix << vcl_endl;

  vnl_svd<float> svd(scatterMatrix);

  const vnl_matrix<float>& V = svd.V();
  
  vcl_cout << "svd: " << vcl_endl
           << svd << vcl_endl;

  //transform points back to real world
  vnl_matrix_fixed<float,4,4> ntransp = norm.get_matrix().transpose();
  vnl_matrix_fixed<float,4,4> eigenvectors;
  eigenvectors = ntransp*V;

  //eigenvectors transformed back into world coordinates
  //each column is an eigenvector sorted by decreasing eigenvalues
  //the normal of the best fitting plane is the last column vector
  vcl_cout << "eigenvectors: " << vcl_endl << eigenvectors << vcl_endl;

  //sanity check. last eigenvector should be equal to default implementation
  vgl_fit_plane_3d<float> planeFitter(pts);
  planeFitter.fit(5);
  vcl_cout << "planeFitter: " << planeFitter.get_plane() << vcl_endl;

  scatterPoints(pts);

  double planeNormal[3] = {eigenvectors(0,3), eigenvectors(1,3), eigenvectors(2,3)};
  {
    const double planeNormalLength = l2norm(planeNormal);
    multiply(planeNormal, 1.0/planeNormalLength, planeNormal);
  }

  vcl_cout << "planeNormal: " << planeNormal << vcl_endl;
  
  double startPoint[3];
  startPoint[0] = centerOfMass[0];
  startPoint[1] = centerOfMass[1];
  startPoint[2] = centerOfMass[2];

  double endPoint[3];
  multiply(planeNormal, 2, endPoint);
  add(startPoint, endPoint, endPoint);

  drawArrow(startPoint, endPoint);

  vcl_cout << "startPoint: " << startPoint << vcl_endl;
  vcl_cout << "endPoint: " << endPoint << vcl_endl;

  {
    vcl_vector<vgl_homg_point_3d<float> > centerOfMassPt;
    centerOfMassPt.push_back(
      vgl_homg_point_3d<float>(centerOfMass[0], centerOfMass[1], centerOfMass[2]));
    scatterPoints(centerOfMassPt,0,0.5,0.5);
  }
  drawAxes();
  vtkBoilerPlate();
  
  return EXIT_SUCCESS;
}
