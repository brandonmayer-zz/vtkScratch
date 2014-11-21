//this is /scratch/planeFit/planeFit1.cxx
#include<vcl_iostream.h>
#include<vcl_vector.h>
#include<vcl_cstdlib.h>
#include<vcl_cmath.h>
#include<random>
#include<chrono>

#include<vgl/vgl_point_3d.h>
#include<vgl/vgl_homg_point_3d.h>
#include<vgl/vgl_homg_plane_3d.h>
#include<vgl/algo/vgl_fit_plane_3d.h>

#include"vtkActor.h"
#include"vtkAxesActor.h"
#include"vtkCellArray.h"
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
#include"vtkTransform.h"
#include"vtkTransformPolyDataFilter.h"
#include"vtkVersion.h"
#include"vtkVertexGlyphFilter.h"
#include"vtkOrientationMarkerWidget.h"
#include"vtkSphereSource.h"
#include"vtkDataSetMapper.h"

#define tol 1e-3

template<class T>
double d_signed(const vgl_homg_plane_3d<T>& plane)
{
  return plane.d()/
    vcl_sqrt(plane.a()*plane.a() +
             plane.b()*plane.b() +
             plane.c()*plane.c());
}

template<class FieldType, template<typename> class PointType = vgl_homg_point_3d>
void scatterPoints(const vcl_vector<PointType<FieldType> >& points,
                    vtkSmartPointer<vtkRenderer> renderer,
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

  const unsigned npts = 20;
  vgl_fit_plane_3d<float> plane_fit;

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<float> dist(-10,10);
  for(unsigned i = 0; i < npts; ++i)
    plane_fit.add_point(dist(generator), dist(generator), 0);

  plane_fit.fit(tol);

  scatterPoints(plane_fit.get_points(), renderer, 1, 0, 0);

  vgl_vector_3d<double> normal = plane_fit.get_plane().normal();
  vcl_cout << "normal = " << normal << vcl_endl;
  vcl_cout << "plane = " << plane_fit.get_plane() << vcl_endl;

  vtkSmartPointer<vtkPlaneSource> plane =
    vtkSmartPointer<vtkPlaneSource>::New();
  plane->SetNormal(normal.x(), normal.y(), normal.z());
  plane->Update();
  
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(plane->GetOutputPort());

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(1,0,0);
  
  plane_fit.clear();
  for(unsigned i = 0; i < npts; ++i)
    plane_fit.add_point(0, dist(generator), dist(generator));

  scatterPoints(plane_fit.get_points(), renderer, 0, 1, 0);
  plane_fit.fit(tol);
  normal = plane_fit.get_plane().normal();
  vcl_cout << "normal = " << normal << vcl_endl;
  vcl_cout << "plane = " << plane_fit.get_plane() << vcl_endl;

  plane = vtkSmartPointer<vtkPlaneSource>::New();
  plane->SetNormal(normal.x(), normal.y(), normal.z());
  plane->Update();

  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(plane->GetOutputPort());

  vtkSmartPointer<vtkActor> actor2 =
    vtkSmartPointer<vtkActor>::New();
  actor2->SetMapper(mapper);
  actor2->GetProperty()->SetColor(0,1,0);

  //plane with offset in pos direction
  plane_fit.clear();
  for(unsigned i = 0; i < npts; ++i)
    plane_fit.add_point(dist(generator), dist(generator), 5);

  scatterPoints(plane_fit.get_points(), renderer, 0, 0, 1);

  plane_fit.fit(tol);
  normal = plane_fit.get_plane().normal();
  vcl_cout << "normal = " << normal << vcl_endl;
  vcl_cout << "plane  = " << plane_fit.get_plane() << vcl_endl;
  
  float d = vcl_abs(plane_fit.get_plane().d()) /
    vcl_sqrt(plane_fit.get_plane().a()*plane_fit.get_plane().a() +
     plane_fit.get_plane().b()*plane_fit.get_plane().b() +
     plane_fit.get_plane().c()*plane_fit.get_plane().c());
  
  vcl_cout << "d = " << d << vcl_endl;

  plane = vtkSmartPointer<vtkPlaneSource>::New();
  plane->SetNormal(normal.x(), normal.y(), normal.z());
  plane->SetCenter(normal.x()*d, normal.y()*d, normal.z()*d);
  vcl_cout << "center = " << normal.x()*d << ", "
           << normal.y()*d << ", "
           << normal.z()*d << vcl_endl;

  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(plane->GetOutputPort());
  vtkSmartPointer<vtkActor> actor3 =
    vtkSmartPointer<vtkActor>::New();
  actor3->SetMapper(mapper);
  actor3->GetProperty()->SetColor(0,0,1);

  //plane with offset in negative direction
  plane_fit.clear();
  for(unsigned i = 0; i < npts; ++i)
    plane_fit.add_point(dist(generator), dist(generator), -3);

  scatterPoints(plane_fit.get_points(), renderer, 0.5, 0.5, 0);
  plane_fit.fit(tol);
  normal = plane_fit.get_plane().normal();
  vcl_cout << "normal = " << normal << vcl_endl;
  vcl_cout << "plane = " << plane_fit.get_plane() << vcl_endl;
  d = d_signed(plane_fit.get_plane());
  vcl_cout << "d = " << d << vcl_endl;
  d = vcl_abs(d);
  vcl_cout << "center = " << normal.x()*d << ", "
           << normal.y()*d << ", "
           << normal.z()*d << vcl_endl;

  plane = vtkSmartPointer<vtkPlaneSource>::New();
  plane->SetNormal(normal.x(), normal.y(), normal.z());
  plane->SetCenter(normal.x()*d, normal.y()*d, normal.z()*d);

  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(plane->GetOutputPort());

  vtkSmartPointer<vtkActor> actor4 =
    vtkSmartPointer<vtkActor>::New();
  actor4->SetMapper(mapper);
  actor4->GetProperty()->SetColor(.5,.5,0);

  renderer->AddActor(actor);
  renderer->AddActor(actor2);
  renderer->AddActor(actor3);
  renderer->AddActor(actor4);
  
  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();
  return EXIT_SUCCESS;
}
