//this is /scratch/planeFit/planeFit2.cxx
#include<vcl_iostream.h>
#include<vcl_vector.h>
#include<vcl_cstdlib.h>
#include<vcl_cmath.h>
#include<vcl_algorithm.h>
#include<random>
#include<chrono>

#include<vgl/vgl_box_3d.h>
#include<vgl/vgl_point_3d.h>
#include<vgl/vgl_homg_point_3d.h>
#include<vgl/vgl_homg_plane_3d.h>
#include<vgl/algo/vgl_fit_plane_3d.h>
#include<vgl/vgl_vector_3d.h>

#include"vtkActor.h"
#include"vtkAxesActor.h"
#include"vtkCellArray.h"
#include"vtkCubeSource.h"
#include"vtkDataSetMapper.h"
#include"vtkGlyph3D.h"
#include"vtkPlaneSource.h"
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
#include"vtkOrientationMarkerWidget.h"
#include"vtkOutlineSource.h"

#define tol 1e-3

template<class T>
double d_signed(const vgl_homg_plane_3d<T>& plane)
{
  return plane.d()/
    vcl_sqrt(plane.a()*plane.a() +
             plane.b()*plane.b() +
             plane.c()*plane.c());
}

template<class T>
double z_from_plane(const T x,
                    const T y,
                    const vgl_homg_plane_3d<T>& plane,
                    const T t = 1e-2)
{
  if(plane.c() < t)
    return 0;
  return (-plane.d() - x*plane.a() - y*plane.b())/plane.c();
}

template<class T>
T z_from_plane(const T x,
               const T y,
               const T a,
               const T b,
               const T c,
               const T d)
{
  return (-d - x*a -y*b)/c;
}

template<class T>
vcl_ostream& operator<<(vcl_ostream& os, const vcl_vector<T>& v)
{
  typename vcl_vector<T>::const_iterator
    vitr = v.begin(), vend = v.end();
  for(;vitr != vend; ++vitr)
    os << *vitr << " ";
  return os;
}

template<class FieldType, template<typename> class PointType>
void scatter_points(const vcl_vector<PointType<FieldType> >& points,
                    vtkSmartPointer<vtkRenderer> renderer,
                    const float r = 1.0,
                    const float g = 1.0,
                    const float b = 1.0)
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
  sphere->SetRadius(0.05);
  
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


template<int D, class FieldType, template<typename> class PointType>
struct pointCompLess
{
  bool operator()(const PointType<FieldType>& rhs, const PointType<FieldType>& lhs)
  {
    return rhs.x() < lhs.x();
  }
};

template<class FieldType, template<typename> class PointType>
struct pointCompLess<1, FieldType, PointType>
{
  bool operator()(const PointType<FieldType>& rhs, const PointType<FieldType>& lhs)
  {
    return rhs.y() < lhs.y();
  }
};

template<class FieldType, template<typename> class PointType>
struct pointCompLess<2, FieldType, PointType>
{
  bool operator()(const PointType<FieldType>& rhs, const PointType<FieldType>& lhs)
  {
    return rhs.z() < lhs.z();
  }
};

template<class FieldType, template<typename> class PointType>
vcl_ostream& operator<<(vcl_ostream& os, const vcl_vector<PointType<FieldType> >& v)
{
  typename vcl_vector<PointType<FieldType> >::const_iterator
    vitr = v.begin(), vend = v.end();
  for(;vitr != vend; ++vitr)
    os << *vitr << "\n";

  return os;
}

template<class FieldType, template<typename> class PointType>
vcl_vector<PointType<FieldType> >
  getCorners(const vcl_vector<PointType<FieldType> >& points)
{
    pointCompLess<0, FieldType, PointType> CompLessX;
    pointCompLess<1, FieldType, PointType> CompLessY;
    pointCompLess<2, FieldType, PointType> CompLessZ;
 
    vcl_vector<PointType<FieldType> > ret;
    
    
    ret.push_back(
      *vcl_max_element(points.begin(), points.end(),
                       pointCompLess<0,FieldType,PointType>()));
    
    ret.push_back(
      *vcl_max_element(points.begin(), points.end(),
                       pointCompLess<1,FieldType,PointType>()));
    
    ret.push_back(
      *vcl_max_element(points.begin(), points.end(),
                       pointCompLess<2,FieldType,PointType>()));

    ret.push_back(
      *vcl_min_element(points.begin(), points.end(),
                       pointCompLess<0, FieldType, PointType>()));
    
    ret.push_back(
      *vcl_min_element(points.begin(), points.end(),
                       pointCompLess<1, FieldType, PointType>()));
    ret.push_back(
      *vcl_min_element(points.begin(), points.end(),
                       pointCompLess<2, FieldType, PointType>()));

    return ret;
}

template<class FieldType, template<typename> class PointType>
vcl_vector<PointType<FieldType> > projectedPoints(const vcl_vector<PointType<FieldType> >& pts, const vgl_homg_plane_3d<FieldType>& plane)
{
  vcl_vector<PointType<FieldType> > ret;
  typename vcl_vector<PointType<FieldType> >::const_iterator
    vitr, vend = pts.end();

  vgl_vector_3d<double> normal = plane.normal();

  const double D = d_signed(plane);
  vcl_cout << "D = " << D << vcl_endl
           <<"plane.d() = " << plane.d() << vcl_endl;
  const double aCoeff = normal.x()*D;
  const double bCoeff = normal.y()*D;
  const double cCoeff = normal.z()*D;

  for(vitr = pts.begin(); vitr != vend; ++vitr)
    ret.push_back(PointType<FieldType>(vitr->x() - aCoeff, vitr->y() - bCoeff, vitr->z() - cCoeff));

  return ret;
}

int main()
{
    //Create a renderer, render window and interactor
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  
  vgl_fit_plane_3d<float> plane_fitter;

  const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<float> noise(0,1);

  const float delta=0.5;
  const float xmin=-5, ymin=-3, xmax=5, ymax=3;
  const int xnpts = vcl_floor((xmax-xmin)/delta) + 1;
  const int ynpts = vcl_floor((ymax-ymin)/delta) + 1;
  vcl_vector<float> xpts(xnpts), ypts(ynpts);
  vcl_vector<vgl_point_3d<float> > pts;
  
  vcl_cout << "xpts = " << xpts << vcl_endl;
  vcl_cout << "ypts = " << ypts << vcl_endl;
  vgl_box_3d<float> bbox;
  for(unsigned i = 0; i < xpts.size(); ++i)
    for(unsigned j = 0; j < ypts.size(); ++j)
    {
      const float x = xmin+delta*i;
      const float y = ymin+delta*j;
      const float z = z_from_plane(
        x, y, (float)(1.0), (float)(1.0), (float)(1.0), (float)(2.0))
        + noise(generator);
      pts.push_back(vgl_point_3d<float>(x,y,z));
      bbox.add(vgl_point_3d<float>(x,y,z));
      vcl_cout << pts[i] << vcl_endl;
    }
  vcl_cout << "npoints = " << pts.size() << vcl_endl;

  vcl_vector<vgl_point_3d<float> > corners =
    getCorners<float, vgl_point_3d>(pts);

  vcl_cout << "corners = " << corners << vcl_endl;
  
  //bounding box
  vtkSmartPointer<vtkCubeSource> cubeSource =
    vtkSmartPointer<vtkCubeSource>::New();
  cubeSource->SetBounds(bbox.min_x(), bbox.max_x(),
                        bbox.min_y(), bbox.max_y(),
                        bbox.min_z(), bbox.max_z());
                        

  vcl_cout << "bbox = " << bbox << vcl_endl;

  vtkSmartPointer<vtkPolyDataMapper> cubeMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  cubeMapper->SetInputConnection(cubeSource->GetOutputPort());

  vtkSmartPointer<vtkActor> cubeActor =
    vtkSmartPointer<vtkActor>::New();
  cubeActor->SetMapper(cubeMapper);
  cubeActor->GetProperty()->SetRepresentationToWireframe();


  vgl_fit_plane_3d<float> planeFitter;
  vcl_vector<vgl_point_3d<float> >::const_iterator
    pitr, pend = pts.end();
  for(pitr = pts.begin(); pitr != pend; ++pitr)
    planeFitter.add_point(pitr->x(), pitr->y(), pitr->z());

  planeFitter.fit(tol);

  vcl_vector<vgl_homg_point_3d<float> > projPts =
    projectedPoints(planeFitter.get_points(), planeFitter.get_plane());

  scatter_points(projPts, renderer, 0, 1, 0);
  
  vgl_vector_3d<double> normal = planeFitter.get_plane().normal();
  vcl_cout << "plane = " << planeFitter.get_plane() << vcl_endl;
  vcl_cout << "normal = " << normal << vcl_endl;
  
  vtkSmartPointer<vtkPlaneSource> planeSource =
    vtkSmartPointer<vtkPlaneSource>::New();

  planeSource->SetNormal(normal.x(), normal.y(), normal.z());
  // planeSource->SetOrigin(bbox.centroid_x(), bbox.centroid_y(), bbox.centroid_z());
  planeSource->SetPoint1(bbox.min_x(), bbox.min_y(), bbox.min_z());
  planeSource->SetPoint2(bbox.max_x(), bbox.max_y(), bbox.max_z());
  
    
  vtkSmartPointer<vtkPolyDataMapper> planeMapper=
    vtkSmartPointer<vtkPolyDataMapper>::New();
  planeMapper->SetInputConnection(planeSource->GetOutputPort());

  vtkSmartPointer<vtkActor> planeActor =
    vtkSmartPointer<vtkActor>::New();
  planeActor->SetMapper(planeMapper);
  
  scatter_points(pts,renderer);

  renderer->AddActor(cubeActor);
  renderer->AddActor(planeActor);

  vtkSmartPointer<vtkAxesActor> axes = 
    vtkSmartPointer<vtkAxesActor>::New();

  vtkSmartPointer<vtkOrientationMarkerWidget> widget =
    vtkSmartPointer<vtkOrientationMarkerWidget>::New();
  widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
  widget->SetOrientationMarker( axes );
  widget->SetInteractor( renderWindowInteractor );
  widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
  widget->SetEnabled( 1 );
  widget->InteractiveOn();
  
  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();
  
  return EXIT_SUCCESS;
}
