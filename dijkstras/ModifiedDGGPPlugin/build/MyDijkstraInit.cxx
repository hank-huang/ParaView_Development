/* #undef PARAVIEW_USE_UNIFIED_BINDINGS */
#define NO_PYTHON_BINDINGS_AVAILABLE
#ifdef NO_PYTHON_BINDINGS_AVAILABLE
#undef PARAVIEW_USE_UNIFIED_BINDINGS
#endif
#ifdef PARAVIEW_USE_UNIFIED_BINDINGS
#include "vtkPython.h"
#include "vtkPythonInterpreter.h"
#endif

#include "vtkClientServerInterpreter.h"

#ifndef PARAVIEW_BUILD_SHARED_LIBS
#define PARAVIEW_BUILD_SHARED_LIBS
#endif
#if defined(PARAVIEW_BUILD_SHARED_LIBS) && defined(_WIN32)
# define VTK_WRAP_CS_EXPORT __declspec(dllexport)
#else
# define VTK_WRAP_CS_EXPORT
#endif

#ifdef PARAVIEW_USE_UNIFIED_BINDINGS
extern "C" void real_initMyDijkstraPython(const char *modulename);

void initMyDijkstraPython()
{
  static const char modulename[] = "MyDijkstraPython";
  real_initMyDijkstraPython(modulename);
}
#endif

extern void vtkMyDGGP_Init(vtkClientServerInterpreter* csi);
extern void vtkDijkstraGraphGeodesicPath1_Init(vtkClientServerInterpreter* csi);


extern "C" void VTK_WRAP_CS_EXPORT MyDijkstra_Initialize(
  vtkClientServerInterpreter *csi)
{
  // silence unreferenced var warning.
  (void) csi;

#ifdef PARAVIEW_USE_UNIFIED_BINDINGS
  if (!vtkPythonInterpreter::IsInitialized())
    {
    vtkPythonInterpreter::Initialize();
    }

  static bool initialized = false;

  if (!initialized)
    {
    initialized = true;
    PyImport_AppendInittab("MyDijkstraPython", initMyDijkstraPython);
    }

  csi->Load("MyDijkstra");
#endif

  vtkMyDGGP_Init(csi);
  vtkDijkstraGraphGeodesicPath1_Init(csi);

}
