/*=========================================================================

   Program: ParaView
   Module:    pqParaViewPlugin.h.in

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
/// This is an auto-generated file. Please do not edit.

/// This file is used to wrap a ParaView plugin.
/* #undef plugin_type_python */
#define plugin_type_servermanager
/* #undef plugin_type_gui */

/// this is set when the plugin has source headers that are wrapped.
#define INITIALIZE_WRAPPING
/* #undef INITIALIZE_EXTRA_CS_MODULES */

#include "vtkConfigure.h"
#include "vtkPVPlugin.h"

#ifdef plugin_type_gui
# include "vtkPVGUIPluginInterface.h"
# include <QObject>
# include <QtPlugin>
#endif

#ifdef plugin_type_servermanager
# include "vtkPVServerManagerPluginInterface.h"
#endif

#ifdef plugin_type_python
# include "vtkPVPythonPluginInterface.h"
#endif

class DijkstraGraphGeodesicPathFilter_Plugin :
#ifdef plugin_type_gui
  public QObject,
  public vtkPVGUIPluginInterface,
#endif

  public vtkPVPlugin

#ifdef plugin_type_servermanager
  , public vtkPVServerManagerPluginInterface
#endif
#ifdef plugin_type_python
  , public vtkPVPythonPluginInterface
#endif

{
#ifdef plugin_type_gui
  Q_OBJECT
  Q_INTERFACES(vtkPVGUIPluginInterface)
  Q_PLUGIN_METADATA(IID "com.kitware/paraview/DijkstraGraphGeodesicPathFilter_Plugin")
//#endif
#endif
public:
  DijkstraGraphGeodesicPathFilter_Plugin();

   // Description:
  // Returns the name for this plugin.
  const char* GetPluginName() override
    {return "DijkstraGraphGeodesicPathFilter"; }

  // Description:
  // Returns the version for this plugin.
  const char* GetPluginVersionString() override
    { return "1.0"; }

  // Description:
  // Returns true if this plugin is required on the server.
  bool GetRequiredOnServer() override
    { return 1; }

  // Description:
  // Returns true if this plugin is required on the client.
  bool GetRequiredOnClient() override
    { return 1; }

  // Description:
  // Returns a ';' separated list of plugin names required by this plugin.
  const char* GetRequiredPlugins() override
    {
    return "";
    }

  const char* GetEULA() override;

  // Description:
  // Provides access to binary resources compiled into the plugin.
  // This is primarily used to compile in icons and compressed help project
  // (qch) files into plugins.
  void GetBinaryResources(std::vector<std::string>& resources) override;

#ifdef plugin_type_servermanager
  // Description:
  // Obtain the server-manager configuration xmls, if any.
  void GetXMLs(std::vector<std::string> &xmls) override;

  // Description:
  // Returns the callback function to call to initialize the interpretor for the
  // new vtk/server-manager classes added by this plugin. Returning NULL is
  // perfectly valid.
  vtkClientServerInterpreterInitializer::InterpreterInitializationCallback
    GetInitializeInterpreterCallback() override;
#endif

#ifdef plugin_type_gui
  /// Returns the list of ParaView-Interfaces provided by this plugin.
  QObjectList interfaces() override;
#endif

#ifdef plugin_type_python
  void GetPythonSourceList(std::vector<std::string>& modules,
    std::vector<std::string>& sources,
    std::vector<int> &package_flags) override;
#endif
};
