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
#include<vgl/algo/vgl_norm_trans_3d.h>
#include<vgl/algo/vgl_fit_plane_3d.h>

#include<vnl/algo/vnl_svd.h>
#include<vnl/vnl_matrix.h>
#include<vnl/vnl_vector.h>
#include<vnl/vnl_vector_fixed.h>

#include"vtkActor.h"
#include"vtkAxesActor.h"
#include"vtkArrowSource.h"
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
#include"vtkMatrix4x4.h"

const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

//Create a renderer, render window and interactor
vtkSmartPointer<vtkRenderer> renderer =
  vtkSmartPointer<vtkRenderer>::New();

vtkSmartPointer<vtkRenderWindow> renderWindow =
  vtkSmartPointer<vtkRenderWindow>::New();

  
vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
  vtkSmartPointer<vtkRenderWindowInteractor>::New();

inline double l2norm(const double v[3])
{
  return vcl_sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

inline void subtract(const double a[3], const double b[3], double out[3])
{
  for(unsigned i = 0; i < 3; ++i)
    out[i] = a[i] - b[i];
}

inline void add(const double a[3], const double b[3], double out[3])
{
  for(unsigned i = 0; i < 3; ++i)
    out[i] = a[i] + b[i];
}

inline void cross(const double x[3], const double y[3], double out[3])
{
  out[0] = x[1] * y[2] - x[2] * y[1];
  out[1] = x[2] * y[0] - x[0] * y[2];
  out[2] = x[0] * y[1] - x[1] * y[0];
}

inline void multiply(const double v[3], const double s, double out[3])
{
  for(unsigned i = 0; i < 3; ++i)
    out[i] = v[i]*s;
}

template<class FieldType, template<typename> class PointType>
void center_of_mass(const vcl_vector<PointType<FieldType> >& pts,
                    FieldType out[3])
{
  for(typename vcl_vector<PointType<FieldType> >::const_iterator
        vitr = pts.begin(); vitr != pts.end(); ++vitr)
  {
    out[0] += vitr->x();
    out[1] += vitr->y();
    out[2] += vitr->z();
  }

  out[0] /= pts.size();
  out[1] /= pts.size();
  out[2] /= pts.size();
}

vcl_ostream& operator<<(vcl_ostream& os, const double v[3])
{
  for(unsigned i = 0; i < 3; ++i)
    os << v[i] << " ";
  return os;
}

template<class FieldType, template<typename> class PointType>
vcl_ostream& operator<<(vcl_ostream& os, const vcl_vector<PointType<FieldType> >& v)
{
  for(typename vcl_vector<PointType<FieldType> >::const_iterator
        vitr = v.begin(); vitr != v.end(); ++vitr)
    os << *vitr << vcl_endl;

  return os;
}

template<class FieldType>
void getScatterMatrix(const vcl_vector<vgl_homg_point_3d<FieldType> >& pts,
                      vnl_matrix<FieldType>& coeff_matrix,
                      vgl_norm_trans_3d<FieldType>& norm)
{
  if(!norm.compute_from_points(pts))
  {
    vcl_cout << "Problem with norm transform." << vcl_endl;
    exit(1);
  }
  
  // compute the matrix A of Ax=b
  FieldType A=0, B=0, C=0, D=0, E=0, F=0, G=0, H=0, I=0;
  for(typename vcl_vector<vgl_homg_point_3d<FieldType> >::const_iterator
        vitr = pts.begin(); vitr != pts.end(); ++vitr)
  {
    const vgl_homg_point_3d<FieldType> pt = norm(*vitr);//normalize
    const FieldType x = pt.x()/pt.w();
    const FieldType y = pt.y()/pt.w();
    const FieldType z = pt.z()/pt.w();
    A += x;
    B += y;
    C += z;
    D += x*x;
    E += y*y;
    F += z*z;
    G += x*y;
    H += y*z;
    I += x*z;
  }

  coeff_matrix.set_size(4,4);
  coeff_matrix(0, 0) = D;
  coeff_matrix(0, 1) = G;
  coeff_matrix(0, 2) = I;
  coeff_matrix(0, 3) = A;

  coeff_matrix(1, 0) = G;
  coeff_matrix(1, 1) = E;
  coeff_matrix(1, 2) = H;
  coeff_matrix(1, 3) = B;

  coeff_matrix(2, 0) = I;
  coeff_matrix(2, 1) = H;
  coeff_matrix(2, 2) = F;
  coeff_matrix(2, 3) = C;

  coeff_matrix(3, 0) = A;
  coeff_matrix(3, 1) = B;
  coeff_matrix(3, 2) = C;
  coeff_matrix(3, 3) = (FieldType)(pts.size());
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
                   const float d = 2,
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
    vitr->set(cx+xynoise(generator),
              cy+xynoise(generator),
              (-d - vitr->x()*a - vitr->y()*b)/c + znoise(generator));


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

void drawPlane(const double origin[3],
               const double normal[3],
               const double ax1[3],
               const double ax2[3],
               const double color[3])
{
  vtkSmartPointer<vtkPlaneSource> planeSource =
    vtkSmartPointer<vtkPlaneSource>::New();

  planeSource->SetOrigin(0, 0, 0);
  planeSource->SetPoint1(ax1[0], ax1[1], ax1[2]);
  planeSource->SetPoint2(ax2[0], ax2[1], ax2[2]);
  planeSource->SetCenter(origin[0], origin[1], origin[2]);
  planeSource->SetNormal(normal[0], normal[1], normal[2]);

  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(planeSource->GetOutputPort());

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();

  renderer->AddActor(actor);
}

void drawPlane(const double origin[3],
               const double normal[3],
               const double ax1[3],
               const double ax2[3])
{
  double color[3] = {1.0,0.0,0.0};
  drawPlane(origin, normal, ax1, ax2, color);
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

// Mostly taken from:
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/GeometricObjects/OrientedArrow
// including vtkMath gave my compiler a hard time, so I re-wrote basic
// functions like cross, add, subtract etc.
void drawArrow(const double begin[3],
               const double end[3],
               const double r = 0.45,
               const double g = 0.45,
               const double b = 0.45)
{

  vtkSmartPointer<vtkArrowSource> arrowSource =
    vtkSmartPointer<vtkArrowSource>::New();

  // Compute a basis
  double normalizedX[3];
  double normalizedY[3];
  double normalizedZ[3];

  // The X axis is a vector from start to end
  subtract(end,begin,normalizedX);
  const double length = l2norm(normalizedX);
  
  {
    //z axis is any vector cross the normalizedX
    double v[3] = {1, 0, 0};
    cross(normalizedX, v, normalizedZ);
  }

  //y is z cross x
  cross(normalizedZ, normalizedX, normalizedY);

  vtkSmartPointer<vtkMatrix4x4> matrix =
    vtkSmartPointer<vtkMatrix4x4>::New();

  matrix->Identity();
  for (unsigned int i = 0; i < 3; i++)
  {
    matrix->SetElement(i, 0, normalizedX[i]);
    matrix->SetElement(i, 1, normalizedY[i]);
    matrix->SetElement(i, 2, normalizedZ[i]);
  }
  

  vtkSmartPointer<vtkTransform> transform =
    vtkSmartPointer<vtkTransform>::New();
  transform->Translate(begin);
  transform->Concatenate(matrix);
  transform->Scale(length,length,length);

#if 0
  // Transform the polydata
  vtkSmartPointer<vtkTransformPolyDataFilter> transformPD = 
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformPD->SetTransform(transform);
  transformPD->SetInputConnection(arrowSource->GetOutputPort());
#endif
 
  //Create a mapper and actor for the arrow
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();

  mapper->SetInputConnection(arrowSource->GetOutputPort());
  actor->SetUserMatrix(transform->GetMatrix());

  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(r,g,b);

  renderer->AddActor(actor);
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
