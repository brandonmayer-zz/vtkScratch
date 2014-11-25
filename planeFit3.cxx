//this is /scratch/planeFit/planeFit3.cxx
#include<vcl_iostream.h>
#include<vcl_vector.h>
#include<vcl_cstdlib.h>
#include<vcl_cmath.h>
#include<vcl_algorithm.h>
#include<vcl_limits.h>
#include<random>
#include<chrono>

#include<vgl/vgl_box_3d.h>
#include<vgl/vgl_point_3d.h>
#include<vgl/vgl_homg_point_3d.h>
#include<vgl/vgl_homg_plane_3d.h>
#include<vgl/vgl_vector_3d.h>

#include<vgl/algo/vgl_fit_plane_3d.h>
#include<vgl/algo/vgl_norm_trans_3d.h>


#include<vnl/vnl_vector.h>
#include<vnl/vnl_matrix.h>
#include<vnl/vnl_matrix_fixed.h>
#include<vnl/algo/vnl_svd.h>

#include"vtkActor.h"
#include"vtkAxesActor.h"
#include"vtkCellArray.h"
#include"vtkCubeSource.h"
#include"vtkDataSetMapper.h"
#include"vtkGlyph3D.h"
#include"vtkOrientationMarkerWidget.h"
#include"vtkOutlineSource.h"
#include"vtkPlaneSource.h"
#include"vtkPlaneWidget.h"
#include"vtkPointData.h"
#include"vtkPolyData.h"
#include"vtkPolyDataMapper.h"
#include"vtkProperty.h"
#include"vtkRenderWindow.h"
#include"vtkRenderWindowInteractor.h"
#include"vtkRenderer.h"
#include"vtkSmartPointer.h"
#include"vtkSphereSource.h"
#include"vtkTransform.h"
#include"vtkTransformPolyDataFilter.h"
#include"vtkType.h"
#include"vtkVersion.h"
#include"vtkVertexGlyphFilter.h"

const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

//Create a renderer, render window and interactor
vtkSmartPointer<vtkRenderer> renderer =
  vtkSmartPointer<vtkRenderer>::New();

vtkSmartPointer<vtkRenderWindow> renderWindow =
  vtkSmartPointer<vtkRenderWindow>::New();

  
vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
  vtkSmartPointer<vtkRenderWindowInteractor>::New();

template<class FieldType, template<typename> class PointType>
vcl_vector<PointType<FieldType> > planarPoints();

void vtkBoilerPlate();

template<class FieldType, template<typename> class PointType>
void scatterPoints(const vcl_vector<PointType<FieldType> >& points,
                   const float r = 1.0,
                   const float g = 0.0,
                   const float b = 0.0);

template<class FieldType, template<typename> class PointType>
vgl_homg_plane_3d<FieldType> fitPlane(const vcl_vector<PointType<FieldType> >& pts);

template<class FieldType, template<typename> class PointType>
vgl_homg_plane_3d<FieldType> myFitPlane(const vcl_vector<PointType<FieldType> >& pts);

template<class FieldType, template<typename> class PointType>
vcl_vector<PointType<FieldType> >projectToPlane(const vcl_vector<PointType<FieldType> >& v,
                                     const vgl_homg_plane_3d<FieldType>& p);


template<class FieldType, template<typename> class PointType>
void getCorners(const vcl_vector<PointType<FieldType> >& projPts,
                FieldType& xmin, FieldType& xmax,
                FieldType& ymin, FieldType& ymax,
                FieldType& zmin, FieldType& zmax);

void drawPlane(const vcl_vector<vgl_homg_point_3d<float> >& v,
               const vgl_homg_plane_3d<float>& p);


template<class FieldType, template<typename> class PointType>
vcl_ostream& operator<<(vcl_ostream& os, const vcl_vector<PointType<FieldType> >& v)
{
  for(typename vcl_vector<PointType<FieldType> >::const_iterator
        vitr = v.begin(); vitr != v.end(); ++vitr)
    os << *vitr << vcl_endl;

  return os;
}

int main()
{
  vcl_vector<vgl_homg_point_3d<float> > pts = planarPoints<float, vgl_homg_point_3d>();
  scatterPoints(pts);

  vgl_homg_plane_3d<float> plane = fitPlane(pts);
  vcl_cout << "plane: " << plane << vcl_endl;
  vcl_cout << "normal: " << plane.normal() << vcl_endl;
  vcl_cout << "D: " << plane.d() /
    vcl_sqrt(plane.a() * plane.a() +
             plane.b() * plane.b() +
             plane.c() * plane.c()) << vcl_endl;

  vcl_cout << "pts:\n " << pts << vcl_endl;

  vgl_homg_plane_3d<float> myPlane = myFitPlane(pts);

  {
    vcl_vector<vgl_homg_point_3d<float> > projPts = projectToPlane(pts, plane);

    scatterPoints(projPts, 0, 1, 0);
  }
  
  drawPlane(pts, plane);
  
  vtkBoilerPlate();
  return EXIT_SUCCESS;
}

template<class FieldType, template<typename> class PointType>
void getCorners(const vcl_vector<PointType<FieldType> >& projPts,
                FieldType& xmin, FieldType& xmax,
                FieldType& ymin, FieldType& ymax,
                FieldType& zmin, FieldType& zmax)
{
  xmin = vcl_numeric_limits<FieldType>::max();
  xmax = -vcl_numeric_limits<FieldType>::max();
  ymin = vcl_numeric_limits<FieldType>::max();
  ymax = -vcl_numeric_limits<FieldType>::max();
  zmin = vcl_numeric_limits<FieldType>::max();
  zmax = -vcl_numeric_limits<FieldType>::max();

  for(typename vcl_vector<PointType<FieldType> >::const_iterator
        pitr = projPts.begin(); pitr != projPts.end(); ++pitr)
  {
    if(xmin > pitr->x())
      xmin = pitr->x();
    
    if(xmax < pitr->x())
      xmax = pitr->x();

    if(ymin > pitr->y())
      ymin = pitr->y();

    if(ymax < pitr->y())
      ymax = pitr->y();

    if(zmin > pitr->z())
      zmin = pitr->z();

    if(zmax < pitr->z())
      zmax = pitr->z();
  }
}

void drawPlane(const vcl_vector<vgl_homg_point_3d<float> >& v,
               const vgl_homg_plane_3d<float>& p)
{
  
  float cx, cy, cz;
  for(typename vcl_vector<vgl_homg_point_3d<float> >::const_iterator
        vitr = v.begin(); vitr != v.end(); ++vitr)
  {
    cx += vitr->x();
    cy += vitr->y();
    cz += vitr->z();
  }
  cx/=v.size();
  cy/=v.size();
  cz/=v.size();

  vcl_cout << "c = " << cx << " " << cy << " " << cz << vcl_endl;

  float xmin, xmax, ymin, ymax, zmin, zmax;
  getCorners(v, xmin, xmax, ymin, ymax, zmin, zmax);

  vcl_cout << xmin << " " << xmax << " "
           << ymin << " " << ymax << " "
           << zmin << " " << zmax << vcl_endl;
  
  vtkSmartPointer<vtkPlaneSource> planeSource =
    vtkSmartPointer<vtkPlaneSource>::New();


  vgl_vector_3d<double> normal = p.normal();
  planeSource->SetOrigin(0,0,0);
  planeSource->SetPoint1(xmin,ymin,0);
  planeSource->SetPoint2(xmax,ymax,0);
  // planeSource->SetXResolution(10);
  // planeSource->SetYResolution(340);
  planeSource->SetCenter(cx,cy,cz);
  // planeSource->SetNormal(0,0,1);
  planeSource->SetNormal(normal.x_, normal.y_, normal.z_);



  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(planeSource->GetOutputPort());
    
  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();

  actor->SetMapper(mapper);
    
  renderer->AddActor(actor);
}

template<class FieldType, template<typename> class PointType>
vcl_vector<PointType<FieldType> >
projectToPlane(const vcl_vector<PointType<FieldType> >& v,
               const vgl_homg_plane_3d<FieldType>& p)
{
  vcl_vector<PointType<FieldType> > ret(v.size());

  const FieldType a = p.a(), b = p.b(), c = p.c(), d = p.d();
  typename vcl_vector<PointType<FieldType> >::iterator ritr = ret.begin();
  for(typename vcl_vector<PointType<FieldType> >::const_iterator
        vitr = v.begin(); vitr != v.end(); ++vitr, ++ritr)
    ritr->set(vitr->x(), vitr->y(), (-d - vitr->x()*a - vitr->y()*b)/c);

  return ret;
}

template<class FieldType, template<typename> class PointType>
vcl_vector<PointType<FieldType> > planarPoints()
{
  vcl_vector<vgl_homg_point_3d<float> > ret;
  std::normal_distribution<float> noise(0,1);
  const float delt = 0.5;
  const float xmin=-5, ymin=-3, xmax=5, ymax=3;
  const int xnpts = vcl_floor((xmax-xmin)/delt);
  const int ynpts = vcl_floor((ymax-ymin)/delt);

  const float a = 1, b = 1, c = 1, d = 2;
  for(unsigned i = 0; i <= xnpts; ++i)
    for(unsigned j = 0; j <= ynpts; ++j)
    {
      const float x = (xmin + i*delt) + noise(generator);
      const float y = (ymin + j*delt) + noise(generator);
      const float z = (((-d - x*a - y*b)/c) + noise(generator));
      ret.push_back(vgl_homg_point_3d<float>(x,y,z,1));
    }

  return ret;
}

template<class FieldType, template<typename> class PointType>
vgl_homg_plane_3d<FieldType> myFitPlane(const vcl_vector<PointType<FieldType> >& pts)
{ 
  const unsigned long long N = pts.size();
  vnl_matrix<FieldType> A(N,3);
  vgl_norm_trans_3d<FieldType> norm;
  norm.compute_from_points(pts);
  {
    unsigned i = 0;
    for(typename vcl_vector<PointType<FieldType> >::const_iterator
          vitr = pts.begin(); vitr != pts.end(); ++vitr, ++i)
    {
      vgl_homg_point_3d<FieldType> pt_norm = norm(*vitr);
      A(i,0) = pt_norm.x();
      A(i,1) = pt_norm.y();
      A(i,2) = pt_norm.z();
    }
  }
  
  vnl_svd<FieldType> svd(A);
  const FieldType min =  svd.sigma_min();
  vcl_cout << "svd.sigma_min() = " << min << vcl_endl;
  vcl_cout << "W_ = " << svd.W() << vcl_endl;
  vcl_cout << "nullvector = " << svd.nullvector() << vcl_endl;

  //retransform points
  vnl_matrix_fixed<FieldType, 4, 4> Hnorm_transp = norm.get_matrix().transpose();
  vnl_vector<FieldType> norm_transformed = Hnorm_transp * svd.nullvector();
  vcl_cout << "norm_transformed = " << norm_transformed << vcl_endl;
  vgl_homg_plane_3d<FieldType> ret(norm_transformed.get(0), norm_transformed.get(1),
                                norm_transformed.get(2), norm_transformed.get(3));
   
  return ret;
}

template<class FieldType, template<typename> class PointType>
vgl_homg_plane_3d<FieldType> fitPlane(const vcl_vector<PointType<FieldType> >& pts)
{
  vgl_fit_plane_3d<FieldType> planeFitter;
  typename vcl_vector<PointType<FieldType> >::const_iterator
    vitr = pts.begin(), vend = pts.end();
  for(;vitr != vend; ++vitr)
    planeFitter.add_point(vitr->x(), vitr->y(), vitr->z());

  planeFitter.fit(10.0, &vcl_cout);

  return planeFitter.get_plane();
}

template<class FieldType, template<typename> class PointType>
void scatterPoints(const vcl_vector<PointType<FieldType> >& points,
                   const float r = 1.0,
                   const float g = 0.0,
                   const float b = 0.0)
{
  vtkSmartPointer<vtkPoints> points_sptr =
    vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkCellArray> cellArray =
    vtkSmartPointer<vtkCellArray>::New();
  
  points_sptr->SetDataTypeToFloat();

  typename vcl_vector<PointType<FieldType> >::const_iterator
    vitr, vend = points.end();
  for(vitr = points.begin(); vitr != vend; ++vitr)
  {
    const vtkIdType pid =
      points_sptr->InsertNextPoint(vitr->x(), vitr->y(), vitr->z());
    cellArray->InsertNextCell(1);
    cellArray->InsertCellPoint(pid);
  }
  
  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();

  polydata->SetPoints(points_sptr);
  polydata->SetVerts(cellArray);

  vtkSmartPointer<vtkGlyph3D> glyph =
    vtkSmartPointer<vtkGlyph3D>::New();

  vtkSmartPointer<vtkSphereSource> sphere =
    vtkSmartPointer<vtkSphereSource>::New();
  sphere->SetRadius(0.25);
  
  glyph->SetInputData(polydata);
  glyph->SetSourceConnection(sphere->GetOutputPort());
  
  vtkSmartPointer<vtkDataSetMapper> mapper=
    vtkSmartPointer<vtkDataSetMapper>::New();

  mapper->SetInputConnection(glyph->GetOutputPort());
  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(r,g,b);

  renderer->AddActor(actor);
}

void vtkBoilerPlate()
{
  renderWindow->AddRenderer(renderer);
  
  renderWindowInteractor->SetRenderWindow(renderWindow);
  vtkSmartPointer<vtkAxesActor> axes = 
    vtkSmartPointer<vtkAxesActor>::New();

  vtkSmartPointer<vtkOrientationMarkerWidget> widget =
    vtkSmartPointer<vtkOrientationMarkerWidget>::New();
  widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
  widget->SetOrientationMarker( axes );
  widget->SetInteractor( renderWindowInteractor );
  widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
  widget->SetEnabled( 1 );
  // widget->InteractiveOn();
  widget->InteractiveOff();
  
  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();
}
