//this is /vtkScratch/util.h
#ifndef UTIL_H_
#define UTIL_H_

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
vcl_ostream& operator<<(vcl_ostream& os, const vcl_vector<PointType<FieldType> >& v)
{
  for(typename vcl_vector<PointType<FieldType> >::const_iterator
        vitr = v.begin(); vitr != v.end(); ++vitr)
    os << *vitr << vcl_endl;

  return os;
}

template<class FieldType, template<typename> class PointType>
vcl_vector<PointType<FieldType> >
sampleCircle(const unsigned npts=25, const float radius=2.0, const float z=0)
{
  vcl_vector<PointType<FieldType> > ret(npts);
  std::normal_distribution<FieldType> normal(0,1);

  for(typename vcl_vector<PointType<FieldType> >::iterator
        vitr = ret.begin(); vitr != ret.end(); ++vitr)
    vitr->set(radius*normal(generator), radius*normal(generator), z, 1);


  return ret;
}

template<class FieldType, template<typename> class PointType>
vcl_vector<PointType<FieldType> >
samplePlanarPoints(const float a = 1,
                   const float b = 1,
                   const float c = 1,
                   const float d = 0,
                   const float cx = 0,
                   const float cy = 0,
                   const float zstd = 1.0,
                   const float xystd = 1.0,
                   const unsigned npts=50)
{
  vcl_vector<PointType<FieldType> > ret(npts);
  
  std::normal_distribution<FieldType> xynoise(0,xystd), znoise(0,zstd);
  for(typename vcl_vector<PointType<FieldType> >::iterator
        vitr = ret.begin(); vitr != ret.end(); ++vitr)
  {
    vitr->set(cx+xynoise(generator),
              cy+xynoise(generator),
              (-d - vitr->x()*a - vitr->y()*b)/c + znoise(generator));
  }

  return ret;
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

void drawAxes(const double cylinderRadius=0.05,
              const double totalLength=5.0)
{
  vtkSmartPointer<vtkAxesActor> axes =
    vtkSmartPointer<vtkAxesActor>::New();
  axes->SetShaftTypeToCylinder();
  axes->SetCylinderRadius(cylinderRadius);
  axes->SetTotalLength(totalLength,totalLength,totalLength);

  renderer->AddActor(axes);
}


//should always be called last
void vtkBoilerPlate()
{
  renderWindow->AddRenderer(renderer);
  
  renderWindowInteractor->SetRenderWindow(renderWindow);

  
  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();
}

#endif
