/*=========================================================================

   Program: ParaView
   Module:    pqParaViewPlugin.cxx.in

   Copyright (c) 2005,2006 Sandia Corporation, Kitware Inc.
   All rights reserved.

   ParaView is a free software; you can redistribute it and/or modify it
   under the terms of the ParaView license version 1.2.

   See License_v1.2.txt for the full ParaView license.
   A copy of this license can be obtained by contacting
   Kitware Inc.
   28 Corporate Drive
   Clifton Park, NY 12065
   USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

========================================================================*/
#include "SurfaceCut_Plugin.h"

#include "vtkObjectFactory.h"
#include "vtkPVPluginLoader.h"

#include "/cis/home/che/Documents/SurfaceCutPlugin/build/vtkSMXML_SurfaceCut.h"

#include "/cis/home/che/Documents/SurfaceCutPlugin/build/SurfaceCut_doc.h"

#if 0 // PLUGIN_HAS_EULA
#include ""
#endif

namespace
{
  // This ensures that when the shared library for this plugin is
  // unloaded during finalization sequence, it notifies the vtkPVPluginLoader
  // so it does not double-dlclose() an already unloaded plugin.
  // This does not affect static builds and hence we don't need to worry about
  // making sure this instance gets created in static builds.
  static class SurfaceCut_PluginCleaner
    {
  public:
    SurfaceCut_PluginCleaner()
      {
      }
    ~SurfaceCut_PluginCleaner()
      {
      // The plugin library is being unloaded.
      // Let the plugin loader know so it doesn't try to unload it again.
      vtkPVPluginLoader::PluginLibraryUnloaded("SurfaceCut");
      }
    } SurfaceCut_PluginCleaner_Instance;
}

//-----------------------------------------------------------------------------
// Used to push a string returns from a function to a vector.
template <class T, class F>
void PushBack(T& vector, F& fun)
{
  char* text = fun();
  vector.push_back(text);
  delete []text;
}

//-----------------------------------------------------------------------------
#ifdef plugin_type_servermanager
# if defined(INITIALIZE_WRAPPING) || defined(INITIALIZE_EXTRA_CS_MODULES)
#include "vtkClientServerInterpreter.h"
#include "vtkClientServerInterpreterInitializer.h"

#   if defined(INITIALIZE_WRAPPING)
extern "C" void SurfaceCut_Initialize(vtkClientServerInterpreter *);
#   endif



extern "C" void VTK_EXPORT SurfaceCut_CombinedInitialize(
  vtkClientServerInterpreter *interp)
{
#   if defined(INITIALIZE_WRAPPING)
  SurfaceCut_Initialize(interp);
#   endif

  // Now initialize any extra kits as requested.
  
}
# endif
#endif

//-----------------------------------------------------------------------------
void SurfaceCut_Plugin::GetBinaryResources(
  std::vector<std::string>& resources)
{
  PushBack(resources, SurfaceCut_doc);

  (void)resources;
}

//-----------------------------------------------------------------------------
#ifdef plugin_type_gui


#endif

//-----------------------------------------------------------------------------
#ifdef plugin_type_servermanager
vtkClientServerInterpreterInitializer::InterpreterInitializationCallback
SurfaceCut_Plugin::GetInitializeInterpreterCallback()
{
#if defined(INITIALIZE_WRAPPING) || defined(INITIALIZE_EXTRA_CS_MODULES)
  return SurfaceCut_CombinedInitialize;
#else
  return NULL;
#endif
}

//-----------------------------------------------------------------------------
void SurfaceCut_Plugin::GetXMLs(std::vector<std::string> &xmls)
{
  PushBack(xmls, SurfaceCutsurfaceCutInterfaces);

  (void)xmls;
}
#endif

//-----------------------------------------------------------------------------
#ifdef plugin_type_gui
QObjectList SurfaceCut_Plugin::interfaces()
{
  QObjectList ifaces;
#ifdef PUSH_BACK_PV_INTERFACES
  PUSH_BACK_PV_INTERFACES(ifaces);
#endif
  return ifaces;
}
#endif

//-----------------------------------------------------------------------------
#ifdef plugin_type_python


void SurfaceCut_Plugin::GetPythonSourceList(std::vector<std::string>& modules,
  std::vector<std::string>& sources,
  std::vector<int> &package_flags)
{
  const char *moduleNames[] = {
    
  };
  char *moduleSources[] = {
    
  };
  const int packageFlags[] = {
    
  };

  int num_modules = sizeof(moduleNames)/sizeof(const char *);
  for (int cc=0; cc < num_modules; cc++)
    {
    modules.push_back(moduleNames[cc]);
    sources.push_back(moduleSources[cc]);
    package_flags.push_back(packageFlags[cc]);

    // free allocated memory.
    delete moduleSources[cc];
    moduleSources[cc] = NULL;
    }
}
#endif

//-----------------------------------------------------------------------------
SurfaceCut_Plugin::SurfaceCut_Plugin()
{
#if !defined(BUILD_SHARED_LIBS) && defined(plugin_type_gui)
  // For static builds, initialize the Qt resources as well as the Qt plugin.
  
  // Qt 5
  Q_IMPORT_PLUGIN(SurfaceCut_Plugin)

#endif
}

//-----------------------------------------------------------------------------
const char* SurfaceCut_Plugin::GetEULA()
{
#if 0 // PLUGIN_HAS_EULA
  return SurfaceCut_EULA;
#else
  return nullptr;
#endif
}

//-----------------------------------------------------------------------------
// Mark this as a ParaView-ServerManager plugin.
PV_PLUGIN_EXPORT(SurfaceCut, SurfaceCut_Plugin)
