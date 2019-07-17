#include <KDApplication.h>
#include <Main.h>
#include <SurfPatches.h>
#include <TDL/ByuSurface.h>
#include <ADL/Array1D.h>
#include <ADL/Array2D.h>
#include <SurfaceDistance.h>

//
// Instantiate Patch dynamic arrays
//
#ifdef GNU_COMPILER
template class DynArray<Patch>;
#endif

#ifdef RS6K
#pragma define(DynArray<Patch>)
#endif

template<>
u_int DynArray<Patch>::m_increment = 10;

#define VERT_UNUSED    0U
#define VERT_PRIMARY   1U
#define VERT_SECONDARY 2U
#define VERT_EDGE      3U

SurfPatches::SurfPatches()
{
}

SurfPatches::~SurfPatches()
{
}

void SurfPatches::buildPatch(int v)
{
   int i,j,n,d,nnei;
   float sumarea;

   int   ns = m_patches.nextSlot();

   m_patches[ns].centerVert = v;
   m_patches[ns].verts.clear();
   m_patches[ns].verts.add(v);

   m_vertUsed[v] = VERT_PRIMARY;

   nnei = m_S->neighborhood()[v].numNeighbors();

   //
   // keep adding a neighborhood's worth of vertices to patch
   // until area is greater or equal to target area
   //
   sumarea = 0.f;
   for(j=1;j<=m_maxDepth && (sumarea < m_sqarea);j++) {
     for(i=0;i<nnei;i++) {
         if(m_S->neighborhood()[v].getDepth(i) != (u_char)j)
		continue;
         n = m_S->neighborhood()[v].getNeighbor(i);
         m_patches[ns].verts.addUnique(n);
         //
         // Set vertex as inside an assigned patch
         //
         m_vertUsed[n] = VERT_SECONDARY;
     }
     sumarea = calcPatchArea(m_patches[ns]);
   }
 
   if(sumarea < m_sqarea)  {
     cerr << "WARNING: SurfPatches:: patch " << ns;
     cerr << " did not meet required size" << endl;
   }

   m_patches[ns].area  = sumarea;
   m_patches[ns].depth = j;
   //
   // tag vertices at edge of neighborhood (depth = j)
   //
   for(i=0;i<nnei;i++) {
     if(m_S->neighborhood()[v].getDepth(i) != (u_char)j)
	 continue;
      n = m_S->neighborhood()[v].getNeighbor(i);
      m_vertUsed[n] = VERT_EDGE;
   }

}

float SurfPatches::calcPatchArea(Patch &P)
{
  static Array1D<u_char> gotVert(m_np);

  int i,a,b,c;
  float sumarea;

  gotVert = 0U;
  for(i=0;i<P.verts.size();i++)
      gotVert[P.verts[i]] = 1U;

  sumarea = 0.;
  for(i=0;i<m_nf;i++) {
    a = m_S->facets()[i][0];
    b = m_S->facets()[i][1];
    c = m_S->facets()[i][2];
    if(gotVert[a] && gotVert[b] && gotVert[c]) {
	sumarea += m_triAreas[i];
    }
  }
  return(sumarea);
}

void SurfPatches::calcAreas()
{
   int i,a,b,c;
   float area;
   for(i=0;i<m_nf;i++) {
      a = m_S->getFacet(i)[0];
      b = m_S->getFacet(i)[1];
      c = m_S->getFacet(i)[2];
      Point pA = m_S->getVert(a);
      Point pB = m_S->getVert(b);
      Point pC = m_S->getVert(c);

      // convert to mm
      pA.set(pA.x()*m_px,pA.y()*m_py,pA.z()*m_pz);
      pB.set(pB.x()*m_px,pB.y()*m_py,pB.z()*m_pz);
      pC.set(pC.x()*m_px,pC.y()*m_py,pC.z()*m_pz);

      Point p1 = pA - pC; // m_S->getVert(a) - m_S->getVert(c);
      Point p2 = pB - pC; // m_S->getVert(b) - m_S->getVert(c);
      area = (p1.cross(p2)).norm()/2.f;
      m_triAreas[i] = area;
   }
}


int SurfPatches::findVertex(Patch &P)
{
   int i,j,v,d,n,nn;

   d = P.depth;

   for(i=0;i<P.verts.size();i++) {
       v = P.verts[i];
       if(m_vertUsed[v] == VERT_EDGE) {
          nn = m_S->neighborhood()[v].numNeighbors();
          for(j=0;j<nn;j++) {
              if(m_S->neighborhood()[v].getDepth(j) != d)
                 continue;
              n = m_S->neighborhood()[v].getNeighbor(j);
              if(m_vertUsed[n] == VERT_UNUSED)
                 return(n);
          } // for j

          for(j=0;j<nn;j++) {
              if(m_S->neighborhood()[v].getDepth(j) != (d-1))
                 continue;
              n = m_S->neighborhood()[v].getNeighbor(j);
              if(m_vertUsed[n] == VERT_UNUSED)
                 return(n);
          } // for j
       } // VERT_EDGE
   }

   // now try rest of surface

   for(i=0;i<m_np;i++) {
      if(m_vertUsed[i] != VERT_UNUSED)
         continue;

      int ns = 0;

      nn = m_S->neighborhood()[i].numNeighbors();
      if(nn<1) continue;

      for(j=0;j<nn;j++) {
          n = m_S->neighborhood()[i].getNeighbor(j); 
          if(m_vertUsed[n] == VERT_PRIMARY) {
             j = nn;
          } else if(m_vertUsed[n] == VERT_SECONDARY)
             ns++;
      }
      //
      // if < 10% of its neighbors arent used, okay
      //
      if((((float)ns)/((float)nn)) < 0.1)
        return(i);
   }

   return(-1);
}


//
// Does vertex 'v' exist in patch 'idx' ?
//
bool SurfPatches::vertExists(int idx, int v)
{
  int i;
  for(i=0;i<m_patches[idx].verts.size();i++)
     if(m_patches[idx].verts[i] == v)
        return(true);
  return(false);
}

bool SurfPatches::calculate(Surface &S, float sq_area,
                            float px, float py, float pz)
{
   int i,j;
   char fname[1024];
   char spref[1024];
  

   if(!&S) {
      cerr << "ERROR: SurfPatches::calculate() invalid Vol/Surface"<<endl;
      return(false);
   }

   DynArray<int>::setIncrement(10);
   
   m_S = &S;
   m_sqarea  = sq_area;
   m_px      = px;
   m_py      = py;
   m_pz      = pz;

   sprintf(spref,"%s",m_S->Filename());
   j = strlen(spref);
   for(i=j-1;i>=0;i--)
      if(spref[i] == '.')
         break;
   if(i>=0) spref[i] = 0;

   m_np = m_S->getNumVert();
   m_nf = m_S->getNumPoly();

   m_patches.clear();

   // KWD : THIS IS TEMP!!!
   if(sq_area < 20.)
      m_maxDepth = 3;
   else
      m_maxDepth = 4;

   m_triAreas.setDim(m_nf);
   m_vertUsed.setDim(m_np);
   m_vertUsed = VERT_UNUSED;

   GenericWorkProc("Calculating neighborhoods...");
   m_S->setWorkingInactive();
   m_S->genNeighborhoods(m_maxDepth);  

   GenericWorkProc("Calculating triangle areas...");
   calcAreas();

   GenericWorkProc("Calculating patches...");

   // Start with first patch at vertex 0
   buildPatch(0);

   bool done = false;

   while(!done) {
     j = m_patches.size()-1;
     i = findVertex(m_patches[j]);
     if(i < 0) {
        done = true;
        break;
     }
     buildPatch(i);
   }

   int cnt=0;
   for(i=0;i<m_np;i++)
      if(m_vertUsed[i] == VERT_UNUSED) 
         cnt++;

   GenericWorkProc("Patching un-used vertices...");
   for(i=0;i<m_np;i++)
      if(m_vertUsed[i] == VERT_UNUSED) 
         buildPatch(i);
   GenericWorkProc(NULL);

   m_triAreas.setDim(0);
   m_vertUsed.setDim(0);

   sprintf(fname,"%s_patches.txt",spref);
   if(Browser->GetFilename("Patches output file",fname,1024,"*.txt")) {
      if(!save(fname)) 
         App->ShowMessage("ERROR: patches did not save. See stdout");
   }

   return(true);
}


bool SurfPatches::load(const char *fname)
{
   int i,j,idx,center,depth,num,n;
   float area;
   char line[512];

   ifstream fp(fname);
   if(!fp || fp.fail()) {
      cerr << "ERROR: SurfPatches::load() ";
      cerr << "could not open " << fname << " for output" << endl;
      return(false);
   }

   line[0] = '#';
   line[1] = 0;
   while((line[0]=='#')||(strlen(line)<1)||fp.fail())
     fp.getline(line,512);

   if(fp.fail())
      return(false);

   sscanf(line,"%d",&num);
   if((num<0)||(num>100000))
       goto badpatches;

   m_patches.setSize(num);
   for(i=0;i<num;i++) {
       fp.getline(line,512);
       if(fp.fail() || strlen(line)<2)
              goto badpatches;
       if(sscanf(line,"%d %d %d %f %d",&idx,
          &center,&depth,&area,&n) != 5)
              goto badpatches;

       m_patches[idx].centerVert = center;
       m_patches[idx].depth      = depth;
       m_patches[idx].area       = area;
       m_patches[idx].verts.setSize(n);
       for(j=0;j<n;j++)
         fp >> m_patches[idx].verts[j];
       fp.getline(line,512); // eat cr
   }

   fp.close();
   return(true);
 
badpatches:
   fp.close();
   cerr << "ERROR: SurfPatches::load() does not appear to be a";
   cerr << " valid patch file" << endl;
   return(false);

}


bool SurfPatches::save(const char *fname)
{
   int i,j;
   ofstream ofp(fname);
   if(!ofp || ofp.fail()) {
      cerr << "ERROR: SurfPatches::save() ";
      cerr << "could not open " << fname << " for output" << endl;
      return(false);
   }

   ofp << "#" << endl;
   ofp << "# BrainWorks SurfPatches file. DO NOT EDIT" << endl;
   ofp << "#" << endl;
   ofp << "# Surface: " << m_S->Filename() << endl;
   ofp << "# Patch Area (mm^2): " << m_sqarea << endl;
   ofp << "#" << endl;
   ofp << "# [FILE FORMAT]" << endl;
   ofp << "#   number_patches" << endl;
   ofp << "#   [then for EACH patch]" << endl;
   ofp << "#   patch_num center_vert nghbd_depth";
   ofp << "tot_area_mm2 nverts" << endl;
   ofp << "# vert1 vert2 vert3 ..." << endl;
   ofp << "# ... vertNNN-1 vertNNN" << endl;
   ofp << "#" << endl;
   ofp << "#" << endl;
   ofp << m_patches.size() << endl;
   for(i=0;i<m_patches.size();i++) {
      ofp << i << " " << m_patches[i].centerVert << " ";
      ofp << m_patches[i].depth << " ";
      ofp << m_patches[i].area << " ";
      ofp  << m_patches[i].verts.size() << endl;
      for(j=0;j<m_patches[i].verts.size();j++) {
         ofp << m_patches[i].verts[j] << " ";
         if(j % 15 == 14) ofp << endl;
      }
      if(j % 15 != 0)
         ofp << endl;
   }
   ofp.close();
   return(true);
}



