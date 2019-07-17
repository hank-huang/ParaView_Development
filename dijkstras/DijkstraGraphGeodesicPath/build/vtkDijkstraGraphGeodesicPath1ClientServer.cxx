// ClientServer wrapper for vtkDijkstraGraphGeodesicPath1 object
//
#define VTK_WRAPPING_CXX
#define VTK_STREAMS_FWD_ONLY
#include "vtkDijkstraGraphGeodesicPath1.h"
#include "vtkSystemIncludes.h"
#include "vtkClientServerInterpreter.h"
#include "vtkClientServerStream.h"


vtkObjectBase *vtkDijkstraGraphGeodesicPath1ClientServerNewCommand(void* /*ctx*/)
{
  return vtkDijkstraGraphGeodesicPath1::New();
}


int VTK_EXPORT vtkDijkstraGraphGeodesicPath1Command(vtkClientServerInterpreter *arlu, vtkObjectBase *ob, const char *method, const vtkClientServerStream& msg, vtkClientServerStream& resultStream, void* /*ctx*/)
{
  vtkDijkstraGraphGeodesicPath1 *op = vtkDijkstraGraphGeodesicPath1::SafeDownCast(ob);
  if(!op)
    {
    vtkOStrStreamWrapper vtkmsg;
    vtkmsg << "Cannot cast " << ob->GetClassName() << " object to vtkDijkstraGraphGeodesicPath1.  "
           << "This probably means the class specifies the incorrect superclass in vtkTypeMacro.";
    resultStream.Reset();
    resultStream << vtkClientServerStream::Error
                 << vtkmsg.str() << 0 << vtkClientServerStream::End;
    return 0;
    }
  (void)arlu;
  if (!strcmp("New",method) && msg.GetNumberOfArguments(0) == 2)
    {
    vtkDijkstraGraphGeodesicPath1  *temp20;
      {
      temp20 = (op)->New();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << (vtkObjectBase *)temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SafeDownCast",method) && msg.GetNumberOfArguments(0) == 3)
    {
    vtkObjectBase  *temp0;
    vtkDijkstraGraphGeodesicPath1  *temp20;
    if(vtkClientServerStreamGetArgumentObject(msg, 0, 2, &temp0, "vtkObjectBase"))
      {
      temp20 = (op)->SafeDownCast(temp0);
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << (vtkObjectBase *)temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("NewInstance",method) && msg.GetNumberOfArguments(0) == 2)
    {
    vtkDijkstraGraphGeodesicPath1  *temp20;
      {
      temp20 = (op)->NewInstance();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << (vtkObjectBase *)temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SetLineType",method) && msg.GetNumberOfArguments(0) == 3)
    {
    int      temp0;
    if(msg.GetArgument(0, 2, &temp0))
      {
      op->SetLineType(temp0);
      return 1;
      }
    }
  if (!strcmp("GetLineType",method) && msg.GetNumberOfArguments(0) == 2)
    {
    int      temp20;
      {
      temp20 = (op)->GetLineType();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SetLineTypeToGeodesic",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->SetLineTypeToGeodesic();
      return 1;
      }
    }
  if (!strcmp("SetLineTypeToSulcus",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->SetLineTypeToSulcus();
      return 1;
      }
    }
  if (!strcmp("SetLineTypeToGyrus",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->SetLineTypeToGyrus();
      return 1;
      }
    }
  if (!strcmp("GetIdList",method) && msg.GetNumberOfArguments(0) == 2)
    {
    vtkIdList  *temp20;
      {
      temp20 = (op)->GetIdList();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << (vtkObjectBase *)temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("StopWhenEndReachedOn",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->StopWhenEndReachedOn();
      return 1;
      }
    }
  if (!strcmp("StopWhenEndReachedOff",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->StopWhenEndReachedOff();
      return 1;
      }
    }
  if (!strcmp("UseScalarWeightsOn",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->UseScalarWeightsOn();
      return 1;
      }
    }
  if (!strcmp("UseScalarWeightsOff",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->UseScalarWeightsOff();
      return 1;
      }
    }
  if (!strcmp("RepelPathFromVerticesOn",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->RepelPathFromVerticesOn();
      return 1;
      }
    }
  if (!strcmp("RepelPathFromVerticesOff",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->RepelPathFromVerticesOff();
      return 1;
      }
    }
  if (!strcmp("SetRepelVertices",method) && msg.GetNumberOfArguments(0) == 3)
    {
    vtkPoints  *temp0;
    if(vtkClientServerStreamGetArgumentObject(msg, 0, 2, &temp0, "vtkPoints"))
      {
      op->SetRepelVertices(temp0);
      return 1;
      }
    }
  if (!strcmp("GetRepelVertices",method) && msg.GetNumberOfArguments(0) == 2)
    {
    vtkPoints  *temp20;
      {
      temp20 = (op)->GetRepelVertices();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << (vtkObjectBase *)temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("GetCumulativeWeights",method) && msg.GetNumberOfArguments(0) == 3)
    {
    vtkDoubleArray  *temp0;
    if(vtkClientServerStreamGetArgumentObject(msg, 0, 2, &temp0, "vtkDoubleArray"))
      {
      op->GetCumulativeWeights(temp0);
      return 1;
      }
    }

  {
    const char* commandName = "vtkGraphGeodesicPath";
    if (arlu->HasCommandFunction(commandName) &&
        arlu->CallCommandFunction(commandName, op, method, msg, resultStream)) { return 1; }
  }
  if(resultStream.GetNumberOfMessages() > 0 &&
     resultStream.GetCommand(0) == vtkClientServerStream::Error &&
     resultStream.GetNumberOfArguments(0) > 1)
    {
    /* A superclass wrapper prepared a special message. */
    return 0;
    }
  vtkOStrStreamWrapper vtkmsg;
  vtkmsg << "Object type: vtkDijkstraGraphGeodesicPath1, could not find requested method: \""
         << method << "\"\nor the method was called with incorrect arguments.\n";
  resultStream.Reset();
  resultStream << vtkClientServerStream::Error
               << vtkmsg.str() << vtkClientServerStream::End;
  vtkmsg.rdbuf()->freeze(0);
  return 0;
}


//-------------------------------------------------------------------------auto
void VTK_EXPORT vtkDijkstraGraphGeodesicPath1_Init(vtkClientServerInterpreter* csi)
{
  static vtkClientServerInterpreter* last = NULL;
  if(last != csi)
    {
    last = csi;
    csi->AddNewInstanceFunction("vtkDijkstraGraphGeodesicPath1", vtkDijkstraGraphGeodesicPath1ClientServerNewCommand);
    csi->AddCommandFunction("vtkDijkstraGraphGeodesicPath1", vtkDijkstraGraphGeodesicPath1Command);
    }
}
