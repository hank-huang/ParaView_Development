// ClientServer wrapper for vtkCurvatures object
//
#define VTK_WRAPPING_CXX
#define VTK_STREAMS_FWD_ONLY
#include "vtkCurvatures.h"
#include "vtkSystemIncludes.h"
#include "vtkClientServerInterpreter.h"
#include "vtkClientServerStream.h"


vtkObjectBase *vtkCurvaturesClientServerNewCommand(void* /*ctx*/)
{
  return vtkCurvatures::New();
}


int VTK_EXPORT vtkCurvaturesCommand(vtkClientServerInterpreter *arlu, vtkObjectBase *ob, const char *method, const vtkClientServerStream& msg, vtkClientServerStream& resultStream, void* /*ctx*/)
{
  vtkCurvatures *op = vtkCurvatures::SafeDownCast(ob);
  if(!op)
    {
    vtkOStrStreamWrapper vtkmsg;
    vtkmsg << "Cannot cast " << ob->GetClassName() << " object to vtkCurvatures.  "
           << "This probably means the class specifies the incorrect superclass in vtkTypeMacro.";
    resultStream.Reset();
    resultStream << vtkClientServerStream::Error
                 << vtkmsg.str() << 0 << vtkClientServerStream::End;
    return 0;
    }
  (void)arlu;
  if (!strcmp("SafeDownCast",method) && msg.GetNumberOfArguments(0) == 3)
    {
    vtkObjectBase  *temp0;
    vtkCurvatures  *temp20;
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
    vtkCurvatures  *temp20;
      {
      temp20 = (op)->NewInstance();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << (vtkObjectBase *)temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("New",method) && msg.GetNumberOfArguments(0) == 2)
    {
    vtkCurvatures  *temp20;
      {
      temp20 = (op)->New();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << (vtkObjectBase *)temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SetCurvatureType",method) && msg.GetNumberOfArguments(0) == 3)
    {
    int      temp0;
    if(msg.GetArgument(0, 2, &temp0))
      {
      op->SetCurvatureType(temp0);
      return 1;
      }
    }
  if (!strcmp("GetCurvatureType",method) && msg.GetNumberOfArguments(0) == 2)
    {
    int      temp20;
      {
      temp20 = (op)->GetCurvatureType();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  if (!strcmp("SetCurvatureTypeToGaussian",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->SetCurvatureTypeToGaussian();
      return 1;
      }
    }
  if (!strcmp("SetCurvatureTypeToMean",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->SetCurvatureTypeToMean();
      return 1;
      }
    }
  if (!strcmp("SetCurvatureTypeToMaximum",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->SetCurvatureTypeToMaximum();
      return 1;
      }
    }
  if (!strcmp("SetCurvatureTypeToMinimum",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->SetCurvatureTypeToMinimum();
      return 1;
      }
    }
  if (!strcmp("InvertMeanCurvatureOn",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->InvertMeanCurvatureOn();
      return 1;
      }
    }
  if (!strcmp("InvertMeanCurvatureOff",method) && msg.GetNumberOfArguments(0) == 2)
    {
      {
      op->InvertMeanCurvatureOff();
      return 1;
      }
    }

  {
    const char* commandName = "vtkPolyDataAlgorithm";
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
  vtkmsg << "Object type: vtkCurvatures, could not find requested method: \""
         << method << "\"\nor the method was called with incorrect arguments.\n";
  resultStream.Reset();
  resultStream << vtkClientServerStream::Error
               << vtkmsg.str() << vtkClientServerStream::End;
  vtkmsg.rdbuf()->freeze(0);
  return 0;
}


//-------------------------------------------------------------------------auto
void VTK_EXPORT vtkCurvatures_Init(vtkClientServerInterpreter* csi)
{
  static vtkClientServerInterpreter* last = NULL;
  if(last != csi)
    {
    last = csi;
    csi->AddNewInstanceFunction("vtkCurvatures", vtkCurvaturesClientServerNewCommand);
    csi->AddCommandFunction("vtkCurvatures", vtkCurvaturesCommand);
    }
}
