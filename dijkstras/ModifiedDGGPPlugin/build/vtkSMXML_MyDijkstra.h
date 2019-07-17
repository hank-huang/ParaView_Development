// Loadable modules
//
// Generated by /export/bofur/che/paraview/paraview_build/bin/vtkkwProcessXML-pv5.6
//
#ifndef __vtkSMXML_MyDijkstra_h
#define __vtkSMXML_MyDijkstra_h

#include <string.h>
#include <cassert>
#include <algorithm>


// From file /cis/home/hhuang/ParaView_Files/dijkstras/ModifiedDGGPPlugin/vtkMyDGGP.xml
static const char* const MyDijkstravtkMyDGGPInterfaces0 =
"<ServerManagerConfiguration>\n"
"  <ProxyGroup name=\"filters\">\n"
"    <!-- ================================================================== -->\n"
"    <SourceProxy name=\"MyDijkstra\" class=\"vtkMyDGGP\" label=\"My Djikstra Graph Geodesic Path\">\n"
"      <Documentation\n"
"         long_help=\"uses Dijkstra Path Geodesic Path to return shortest path between\n"
"         user selected points\"\n"
"         short_help=\"Shortest path between string of points.\">\n"
"      </Documentation>\n"
"\n"
"      <InputProperty\n"
"         name=\"Input\"\n"
"	       port_index=\"0\"\n"
"         command=\"SetInputConnection\">\n"
"        <ProxyGroupDomain name=\"groups\">\n"
"          <Group name=\"sources\"/>\n"
"          <Group name=\"filters\"/>\n"
"        </ProxyGroupDomain>\n"
"        <DataTypeDomain name=\"input_type\">\n"
"          <DataType value=\"vtkPolyData\"/>\n"
"        </DataTypeDomain>\n"
"        <Documentation>\n"
"          Set the source data set.\n"
"        </Documentation>\n"
"      </InputProperty>\n"
"\n"
"      <InputProperty\n"
"         name=\"Selection\"\n"
"         port_index=\"1\"\n"
"         command=\"SetInputConnection\">\n"
"        <ProxyGroupDomain name=\"groups\">\n"
"          <Group name=\"sources\"/>\n"
"          <Group name=\"filters\"/>\n"
"        </ProxyGroupDomain>\n"
"        <DataTypeDomain name=\"input_type\">\n"
"          <DataType value=\"vtkUnstructuredGrid\"/>\n"
"        </DataTypeDomain>\n"
"        <Documentation>\n"
"          Set the extracted selection data set.\n"
"        </Documentation>\n"
"      </InputProperty>\n"
"\n"
"      <IntVectorProperty\n"
"         name=\"LineType\"\n"
"         command=\"SetLineType\"\n"
"         number_of_elements=\"1\"\n"
"         default_values=\"0\" >\n"
"       <EnumerationDomain name=\"enum\">\n"
"         <Entry value=\"0\" text=\"Geodesic\"/>\n"
"         <Entry value=\"1\" text=\"Sulcus\"/>\n"
"         <Entry value=\"2\" text=\"Gyrus\"/>\n"
"       </EnumerationDomain>\n"
"       <Documentation>\n"
"         This propery specifies which type of curvature to compute.\n"
"       </Documentation>\n"
"     </IntVectorProperty>\n"
"\n"
"\n"
"\n"
"    </SourceProxy>\n"
"  </ProxyGroup>\n"
"</ServerManagerConfiguration>\n"
"\n";
// Get single string
char* MyDijkstravtkMyDGGPInterfaces()
{

  const size_t len0 = strlen(MyDijkstravtkMyDGGPInterfaces0);
  size_t len = ( 0
    + len0 );
  char* res = new char[ len + 1];
  size_t offset = 0;
  std::copy(MyDijkstravtkMyDGGPInterfaces0, MyDijkstravtkMyDGGPInterfaces0 + len0, res + offset); offset += len0;
  assert(offset == len);
  res[offset] = 0;
  return res;
}



#endif