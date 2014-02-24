

// Qt includes
#include <QApplication>
#include <QDebug>
#include <QVBoxLayout>

// CTK includes
#include <ctkSliderWidget.h>
#include <ctkVTKSliceView.h>

// MRML includes
#include <vtkImageResliceMask.h>
#include <vtkMRMLColorTableNode.h>

// VTK includes
#include <vtkGeneralTransform.h>
#include <vtkImageAppendComponents.h>
#include <vtkImageBlend.h>
#include <vtkImageCast.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageData.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageLogic.h>
#include <vtkImageMapToColors.h>
#include <vtkImageMapToWindowLevelColors.h>
#include <vtkImageMathematics.h>
#include <vtkImageThreshold.h>
#include <vtkLookupTable.h>
#include <vtkMath.h>
#include <vtkMatrix4x4.h>
#include <vtkNew.h>

// vtkITK includes
#include <vtkITKArchetypeImageSeriesScalarReader.h>

// STD includes
#include <cstdlib>


namespace
{
  double Origin[3];
  double Spacing[3];
  double Dirs[3][3];

  vtkNew<vtkMatrix4x4> SliceToRAS;

  // vtkMRMLSliceNode
  vtkNew<vtkMatrix4x4> XYToSlice;
  vtkNew<vtkMatrix4x4> XYToRAS;
  double FieldOfView[3];
  int Dimensions[3];
  double XYZOrigin[3];
  double SliceSpacing[3];

  // Pipeline - vtkMRMLScalarVolumeDisplayNode
  vtkNew<vtkImageCast>                   ResliceAlphaCast;
  vtkNew<vtkImageLogic>                  AlphaLogic;
  vtkNew<vtkImageMapToColors>            MapToColors;
  vtkNew<vtkImageAppendComponents>       AppendComponents;

  vtkNew<vtkImageExtractComponents>      ExtractRGB;
  vtkNew<vtkImageExtractComponents>      ExtractAlpha;
  vtkNew<vtkImageMathematics>            MultiplyAlpha;
  vtkNew<vtkImageMapToWindowLevelColors> MapToWindowLevelColors;

  // Pipeline - vtkMRMLSliceLayerLogic
  vtkNew<vtkImageResliceMask> Reslice;
  vtkNew<vtkImageThreshold>   Threshold;
  vtkNew<vtkGeneralTransform> XYToIJKTransform;

  // Pipeline - vtkMRMLSliceLogic
  vtkNew<vtkImageBlend> Blend;
}

//----------------------------------------------------------------------------
void GetIJKToRASMatrix(vtkMatrix4x4* mat)
{
  // this is the full matrix including the spacing and origin
  mat->Identity();
  int row, col;
  for (row=0; row<3; row++)
    {
    for (col=0; col<3; col++)
      {
      mat->SetElement(row, col, ::Spacing[col] * ::Dirs[row][col]);
      }
    mat->SetElement(row, 3, ::Origin[row]);
    }
}

//---------------------------------------------------------------------------
void GetRASToIJKMatrix(vtkMatrix4x4* mat)
{
  GetIJKToRASMatrix(mat);
  mat->Invert();
}

//---------------------------------------------------------------------------
void SetIJKToRASMatrix(vtkMatrix4x4* argMat)
{
  if (argMat == NULL)
    {
    return;
    }
  vtkNew<vtkMatrix4x4> mat;
  mat->DeepCopy(argMat);

  // normalize direction vectors
  int col;
  for (col=0; col<3; col++)
    {
    double len =0;
    int row;
    for (row=0; row<3; row++)
      {
      len += mat->GetElement(row, col) * mat->GetElement(row, col);
      }
    len = sqrt(len);
    ::Spacing[col] = len;
    for (row=0; row<3; row++)
      {
      mat->SetElement(row, col,  mat->GetElement(row, col)/len);
      }
    }

  for (int row=0; row<3; row++)
    {
    for (int col=0; col<3; col++)
      {
      ::Dirs[row][col] = mat->GetElement(row, col);
      }
    ::Origin[row] = mat->GetElement(row, 3);
    }
}

//---------------------------------------------------------------------------
void GetRASBounds(vtkImageData *volumeImage, double bounds[6])
{
  vtkMath::UninitializeBounds(bounds);

  if (!volumeImage)
    {
    return;
    }

  //
  // Get the size of the volume in RAS space
  // - map the size of the volume in IJK into RAS
  // - map the middle of the volume to RAS for the center
  //   (IJK space always has origin at first pixel)
  //

  vtkNew<vtkGeneralTransform> transform;
  transform->PostMultiply();
  transform->Identity();

  vtkNew<vtkMatrix4x4> ijkToRAS;
  GetIJKToRASMatrix(ijkToRAS.GetPointer());

  transform->Concatenate(ijkToRAS.GetPointer());

  /*
  vtkMRMLTransformNode *transformNode = this->GetParentTransformNode();

  if ( transformNode )
    {
    vtkNew<vtkGeneralTransform> worldTransform;
    worldTransform->Identity();
    //transformNode->GetTransformFromWorld(worldTransform);
    transformNode->GetTransformToWorld(worldTransform.GetPointer());
    transform->Concatenate(worldTransform.GetPointer());
    }
  */

  int dimensions[3];
  int i,j,k;
  volumeImage->GetDimensions(dimensions);
  double doubleDimensions[4], *rasHDimensions;
  double minBounds[3], maxBounds[3];

  for ( i=0; i<3; i++)
    {
    minBounds[i] = 1.0e10;
    maxBounds[i] = -1.0e10;
    }
  for ( i=0; i<2; i++)
    {
    for ( j=0; j<2; j++)
      {
      for ( k=0; k<2; k++)
        {
        doubleDimensions[0] = i*(dimensions[0]) - 0.5;
        doubleDimensions[1] = j*(dimensions[1]) - 0.5 ;
        doubleDimensions[2] = k*(dimensions[2]) - 0.5;
        doubleDimensions[3] = 1;
        rasHDimensions = transform->TransformDoublePoint( doubleDimensions);
        for (int n=0; n<3; n++) {
          if (rasHDimensions[n] < minBounds[n])
            {
            minBounds[n] = rasHDimensions[n];
            }
          if (rasHDimensions[n] > maxBounds[n])
            {
            maxBounds[n] = rasHDimensions[n];
            }
          }
        }
      }
    }

   for ( i=0; i<3; i++)
    {
    bounds[2*i]   = minBounds[i];
    bounds[2*i+1] = maxBounds[i];
    }
}

// --------------------------------------------------------------------------
void GetVolumeRASBox(vtkImageData* volumeImage, double rasDimensions[3], double rasCenter[3])
{
  rasCenter[0] = rasDimensions[0] = 0.0;
  rasCenter[1] = rasDimensions[1] = 0.0;
  rasCenter[2] = rasDimensions[2] = 0.0;

  double bounds[6];
  GetRASBounds(volumeImage, bounds);

  for (int i=0; i<3; i++)
    {
    rasDimensions[i] = bounds[2*i+1] - bounds[2*i];
    rasCenter[i] = 0.5*(bounds[2*i+1] + bounds[2*i]);
  }
}

// --------------------------------------------------------------------------
void GetVolumeSliceBounds(vtkImageData* volumeImage, vtkMatrix4x4* sliceToRAS, double sliceBounds[6])
{
  sliceBounds[0] = sliceBounds[1] = 0.0;
  sliceBounds[2] = sliceBounds[3] = 0.0;
  sliceBounds[4] = sliceBounds[5] = 0.0;

  double rasDimensions[3], rasCenter[3];

  GetVolumeRASBox(volumeImage, rasDimensions, rasCenter);

  //
  // figure out how big that volume is on this particular slice plane
  //
  vtkNew<vtkMatrix4x4> rasToSlice;
  rasToSlice->DeepCopy(sliceToRAS);
  rasToSlice->SetElement(0, 3, 0.0);
  rasToSlice->SetElement(1, 3, 0.0);
  rasToSlice->SetElement(2, 3, 0.0);
  rasToSlice->Invert();

  double minBounds[3], maxBounds[3];
  double rasCorner[4], sliceCorner[4];
  int i,j,k;

  for ( i=0; i<3; i++)
    {
    minBounds[i] = 1.0e10;
    maxBounds[i] = -1.0e10;
    }
  for ( i=-1; i<=1; i=i+2)
    {
    for ( j=-1; j<=1; j=j+2)
      {
      for ( k=-1; k<=1; k=k+2)
        {
        rasCorner[0] = rasCenter[0] + i * rasDimensions[0] / 2.;
        rasCorner[1] = rasCenter[1] + j * rasDimensions[1] / 2.;
        rasCorner[2] = rasCenter[2] + k * rasDimensions[2] / 2.;
        rasCorner[3] = 1.;

        rasToSlice->MultiplyPoint( rasCorner, sliceCorner );

        for (int n=0; n<3; n++) {
          if (sliceCorner[n] < minBounds[n])
            {
            minBounds[n] = sliceCorner[n];
            }
          if (sliceCorner[n] > maxBounds[n])
            {
            maxBounds[n] = sliceCorner[n];
            }
          }
        }
      }
    }

  // ignore homogeneous coordinate
  sliceBounds[0] = minBounds[0];
  sliceBounds[1] = maxBounds[0];
  sliceBounds[2] = minBounds[1];
  sliceBounds[3] = maxBounds[1];
  sliceBounds[4] = minBounds[2];
  sliceBounds[5] = maxBounds[2];
}

// --------------------------------------------------------------------------
void ComputeAndSetVolumeSliceSpacing()
{
  // Based on vtkMRMLSliceLogic::GetVolumeSliceSpacing()

  vtkNew<vtkMatrix4x4> ijkToRAS;
  vtkNew<vtkMatrix4x4> rasToSlice;
  vtkNew<vtkMatrix4x4> ijkToSlice;

  GetIJKToRASMatrix(ijkToRAS.GetPointer());

  // Apply the transform, if it exists
  //vtkMRMLTransformNode *transformNode = volumeNode->GetParentTransformNode();

  rasToSlice->DeepCopy(::SliceToRAS.GetPointer());
  rasToSlice->Invert();

  ijkToSlice->Multiply4x4(rasToSlice.GetPointer(), ijkToRAS.GetPointer(), ijkToSlice.GetPointer());

  double invector[4] = {1., 1., 1., 0.};
  double spacing[4];
  ijkToSlice->MultiplyPoint(invector, spacing);
  for (int i = 0; i < 3; ++i)
    {
    ::SliceSpacing[i] = fabs(spacing[i]);
    }
}

// --------------------------------------------------------------------------
void UpdateMatrices()
{
  double spacing[3];
  unsigned int i;
  vtkNew<vtkMatrix4x4> xyToSlice;
  vtkNew<vtkMatrix4x4> xyToRAS;

  // the mapping from XY output slice pixels to Slice Plane coordinate
  xyToSlice->Identity();
  if (::Dimensions[0] > 0 &&
      ::Dimensions[1] > 0 &&
      ::Dimensions[2] > 0)
    {
    for (i = 0; i < 3; i++)
      {
      spacing[i] = ::FieldOfView[i] / ::Dimensions[i];
      xyToSlice->SetElement(i, i, spacing[i]);
      xyToSlice->SetElement(i, 3, -::FieldOfView[i] / 2. + ::XYZOrigin[i]);
      }
    //std::cout << "FieldOfView[2] = " << ::FieldOfView[2] << ", Dimensions[2] = " << ::Dimensions[2] << std::endl;
    //xyToSlice->SetElement(2, 2, 1.);

    xyToSlice->SetElement(2, 3, 0.);
    }

    // the mapping from slice plane coordinates to RAS
    // (the Orienation as in Axial, Sagittal, Coronal)
    //
    // The combined transform:
    //
    // | R | = [Slice to RAS ] [ XY to Slice ]  | X |
    // | A |                                    | Y |
    // | S |                                    | Z |
    // | 1 |                                    | 1 |
    //
    // or
    //
    // RAS = XYToRAS * XY
    //
    vtkMatrix4x4::Multiply4x4(::SliceToRAS.GetPointer(), xyToSlice.GetPointer(), xyToRAS.GetPointer());

    // check to see if the matrix actually changed
    //if ( !Matrix4x4AreEqual (xyToRAS.GetPointer(), ::XYToRAS.GetPointer()) )
    //  {
      ::XYToSlice->DeepCopy(xyToSlice.GetPointer());
      ::XYToRAS->DeepCopy(xyToRAS.GetPointer());
    //  }
}

//----------------------------------------------------------------------------
void UpdateTransforms()
{
  vtkNew<vtkMatrix4x4> xyToIJK;
  xyToIJK->Identity();

  ::XYToIJKTransform->Identity();
  ::XYToIJKTransform->PostMultiply();

  vtkMatrix4x4::Multiply4x4(::XYToRAS.GetPointer(), xyToIJK.GetPointer(), xyToIJK.GetPointer());
  ::XYToIJKTransform->Concatenate(xyToIJK.GetPointer());

  vtkNew<vtkMatrix4x4> rasToIJK;
  ::GetRASToIJKMatrix(rasToIJK.GetPointer());

  ::XYToIJKTransform->Concatenate(rasToIJK.GetPointer());

  ::Reslice->SetResliceTransform( ::XYToIJKTransform.GetPointer() );

  ::Reslice->SetOutputExtent( 0, ::Dimensions[0]-1,
                              0, ::Dimensions[1]-1,
                              0, ::Dimensions[2]-1);
}

// --------------------------------------------------------------------------
class SSPWidget : public QWidget
{
  Q_OBJECT
public:
  SSPWidget(QWidget* parent = 0):QWidget(parent)
    {
    this->OldOffset = 0.0;
    this->VLayout  = new QVBoxLayout(this);
    this->Slider = new ctkSliderWidget;
    this->VLayout->addWidget(this->Slider);
    this->SliceView = new ctkVTKSliceView;
    this->VLayout->addWidget(this->SliceView);

    connect(this->Slider, SIGNAL(valueChanged(double)), this, SLOT(onSliderValueChanged(double)));
    connect(this->SliceView, SIGNAL(resized(QSize)),
          this, SLOT(onWindowResized(QSize)));
    }
  virtual ~SSPWidget(){}

  void setSliceOffsetRange(double min, double max)
    {
    this->Slider->setRange(min, max);
    }

  void setImageData(vtkImageData* imageData)
    {
    this->SliceView->setImageData(imageData);
    }

public slots:

  void onWindowResized(const QSize& newSize)
    {
    double newWidth = newSize.width();
    double newHeight = newSize.height();

    // The following was previously in SliceSWidget.tcl
    double sliceStep = ::SliceSpacing[2];
    int oldDimensions[3];
    oldDimensions[0] = ::Dimensions[0];
    oldDimensions[1] = ::Dimensions[1];
    oldDimensions[2] = ::Dimensions[2];

    double oldFOV[3];
    oldFOV[0] = ::FieldOfView[0];
    oldFOV[1] = ::FieldOfView[1];
    oldFOV[2] = ::FieldOfView[2];

    double scalingX = (newWidth != 0 && oldDimensions[0] != 0 ? newWidth / oldDimensions[0] : 1.);
    double scalingY = (newHeight != 0 && oldDimensions[1] != 0 ? newHeight / oldDimensions[1] : 1.);

    double magnitudeX = (scalingX >= 1. ? scalingX : 1. / scalingX);
    double magnitudeY = (scalingY >= 1. ? scalingY : 1. / scalingY);

    double newFOV[3];
    if (magnitudeX < magnitudeY)
      {
      newFOV[0] = oldFOV[0];
      newFOV[1] = oldFOV[1] * scalingY / scalingX;
      }
    else
      {
      newFOV[0] = oldFOV[0] * scalingX / scalingY;
      newFOV[1] = oldFOV[1];
      }
    newFOV[2] = sliceStep * oldDimensions[2];
    double windowAspect = (newWidth != 0. ? newHeight / newWidth : 1.);
    double planeAspect = (newFOV[0] != 0. ? newFOV[1] / newFOV[0] : 1.);
    if (windowAspect != planeAspect)
      {
      newFOV[0] = (windowAspect != 0. ? newFOV[1] / windowAspect : newFOV[0]);
      }

    ::Dimensions[0] = newWidth;
    ::Dimensions[1] = newHeight;
    ::Dimensions[2] = oldDimensions[2];

    ::FieldOfView[0] = newFOV[0];
    ::FieldOfView[1] = newFOV[1];
    ::FieldOfView[2] = newFOV[2];

    UpdateMatrices();
    UpdateTransforms();

    this->SliceView->scheduleRender();
    }

  void onSliderValueChanged(double offset)
    {
    //
    // Set the Offset
    // - get the current translation in RAS space and convert it to Slice space
    //   by transforming it by the invers of the upper 3x3 of SliceToRAS
    // - replace the z value of the translation with the new value given by the slider
    // - this preserves whatever translation was already in place
    //

    double oldOffset = this->OldOffset;
    if (fabs(offset - oldOffset) <= 1.0e-6)
      {
      return;
      }

    vtkNew<vtkMatrix4x4> sliceToRAS;
    sliceToRAS->DeepCopy( ::SliceToRAS.GetPointer() );
    for (int i = 0; i < 3; i++)
      {
      sliceToRAS->SetElement( i, 3, 0.0 );  // Zero out the tranlation portion
      }
    vtkNew<vtkMatrix4x4> sliceToRASInverted; // inverse sliceToRAS
    sliceToRASInverted->DeepCopy(sliceToRAS.GetPointer());
    sliceToRASInverted->Invert();
    double v1[4], v2[4], v3[4];
    for (int i = 0; i < 4; i++)
      { // get the translation back as a vector
      v1[i] = ::SliceToRAS->GetElement( i, 3 );
      }
    // bring the translation into slice space
    // and overwrite the z part
    sliceToRASInverted->MultiplyPoint(v1, v2);

    v2[2] = offset;

    // Now bring the new translation vector back into RAS space
    sliceToRAS->MultiplyPoint(v2, v3);

    // if the translation has changed, update the rest of the matrices
    double eps=1.0e-6;
    if ( fabs(v1[0] - v3[0]) > eps ||
         fabs(v1[1] - v3[1]) > eps ||
         fabs(v1[2] - v3[2]) > eps )
      {
      // copy new translation into sliceToRAS
      for (int i = 0; i < 4; i++)
        {
        sliceToRAS->SetElement( i, 3, v3[i] );
        }
      ::SliceToRAS->DeepCopy(sliceToRAS.GetPointer());
      UpdateMatrices();
      UpdateTransforms();
      }
    this->OldOffset = offset;

    this->SliceView->scheduleRender();
    }

private:
  QVBoxLayout* VLayout;
  ctkSliderWidget* Slider;
  ctkVTKSliceView* SliceView;
  double OldOffset;
};

// --------------------------------------------------------------------------
int main(int argc, char*argv[])
{
  QApplication app(argc, argv);

  const char* fullName = "/tmp/Slicer/RemoteIO/MR-head.nrrd";
  // Hardcoded value - Usually computed in vtkMRMLScalarVolumeDisplayNode::CalculateAutoLevels()
  double window = 128;
  double level = 67;
  double lowerThreshold = -600;
  double higherThreshold = 600;

  vtkNew<vtkITKArchetypeImageSeriesScalarReader> reader;
  reader->SetUseOrientationFromFile(1);
  reader->SetSingleFile(1);
  reader->ResetFileNames();
  reader->SetArchetype(fullName);

  // Use native origin
  reader->SetOutputScalarTypeToNative();
  reader->SetDesiredCoordinateOrientationToNative();
  reader->SetUseNativeOriginOn();
  reader->Update();

  vtkNew<vtkImageChangeInformation> ici;
  ici->SetInput(reader->GetOutput());
  ici->SetOutputSpacing( 1, 1, 1 );
  ici->SetOutputOrigin( 0, 0, 0 );
  ici->Update();

  vtkNew<vtkMatrix4x4> ijkToRas;
  ijkToRas->DeepCopy(reader->GetRasToIjkMatrix());
  ijkToRas->Invert();
  SetIJKToRASMatrix(ijkToRas.GetPointer());

  if (ici->GetOutput() == NULL)
    {
    std::cerr << "Cannot read file: " << fullName << std::endl;
    return EXIT_FAILURE;
    }

  vtkNew<vtkImageData> iciOutputCopy;
  iciOutputCopy->ShallowCopy(ici->GetOutput());

  vtkImageData * imageData = iciOutputCopy.GetPointer();

  SSPWidget widget;

  ::FieldOfView[0] = 250.0;
  ::FieldOfView[1] = 250.0;
  ::FieldOfView[2] = 1.0;

  ::Dimensions[0] = 256;
  ::Dimensions[1] = 256;
  ::Dimensions[2] = 1;

  ::XYZOrigin[0] = ::XYZOrigin[1] = XYZOrigin[2] = 0;

  ::SliceSpacing[0] = ::SliceSpacing[1] = ::SliceSpacing[2] = 1;

  // vtkMRMLSliceLogic::UpdateSliceNodes -> vtkMRMLSliceLogic::UpdateSliceNode \
  //     -> vtkMRMLSliceLogic::UpdateSliceNodeFromLayout -> vtkMRMLSliceNode::SetOrientationToAxial
  ::SliceToRAS->Identity();
  // Px -> Patient Left
  ::SliceToRAS->SetElement(0, 0, -1.0);
  ::SliceToRAS->SetElement(1, 0,  0.0);
  ::SliceToRAS->SetElement(2, 0,  0.0);
  // Py -> Patient Anterior
  ::SliceToRAS->SetElement(0, 1,  0.0);
  ::SliceToRAS->SetElement(1, 1,  1.0);
  ::SliceToRAS->SetElement(2, 1,  0.0);
  // Pz -> Patient Inferior
  ::SliceToRAS->SetElement(0, 2,  0.0);
  ::SliceToRAS->SetElement(1, 2,  0.0);
  ::SliceToRAS->SetElement(2, 2,  1.0);

  // GetVolumeSliceBounds
  double sliceBounds[6] = {0, -1, 0, -1, 0, -1};
  GetVolumeSliceBounds(imageData, ::SliceToRAS.GetPointer(), sliceBounds);
  std::cout << "min:" << sliceBounds[4] << " / max: " << sliceBounds[5] << std::endl;
  widget.setSliceOffsetRange(sliceBounds[4], sliceBounds[5]);

  // Slice spacing
  ComputeAndSetVolumeSliceSpacing();

  // ------------------------------------------
  // Pipeline - vtkMRMLScalarVolumeDisplayNode
  // ------------------------------------------
  ::MapToWindowLevelColors->SetOutputFormatToLuminance();
  ::MapToWindowLevelColors->SetWindow(window); // default: 256.
  ::MapToWindowLevelColors->SetLevel(level); // default: 128.

  //      MapToWindowLevelColors -> MapToColors
  ::MapToColors->SetOutputFormatToRGBA();
  ::MapToColors->SetInputConnection(::MapToWindowLevelColors->GetOutputPort() );

  //      MapToColors -> ExtractRGB
  ::ExtractRGB->SetInputConnection(::MapToColors->GetOutputPort());
  ::ExtractRGB->SetComponents(0,1,2);

  //      MapToColors -> ExtractAlpha
  ::ExtractAlpha->SetInputConnection(::MapToColors->GetOutputPort());
  ::ExtractAlpha->SetComponents(3);

  ::Threshold->ReplaceInOn();
  ::Threshold->SetInValue(255);
  ::Threshold->ReplaceOutOn();
  ::Threshold->SetOutValue(255);
  ::Threshold->SetOutputScalarTypeToUnsignedChar();
  ::Threshold->ThresholdBetween(lowerThreshold, higherThreshold); // default: VTK_SHORT_MIN / VTK_SHORT_MAX
  ::ResliceAlphaCast->SetOutputScalarTypeToUnsignedChar();

  //      [ExtractAlpha, ResliceAlphaCast] -> MultiplyAlpha
  ::MultiplyAlpha->SetInputConnection(0, ::ExtractAlpha->GetOutputPort() );
  ::MultiplyAlpha->SetInputConnection(1, ::ResliceAlphaCast->GetOutputPort() );
  ::MultiplyAlpha->SetOperationToMultiply();

  //      [Threshold, MultiplyAlpha] -> AlphaLogic
  ::AlphaLogic->SetOperationToAnd();
  ::AlphaLogic->SetOutputTrueValue(255);
  ::AlphaLogic->SetInputConnection(0, ::Threshold->GetOutputPort() );
  //::AlphaLogic->SetInputConnection(1, ::Threshold->GetOutputPort() );
  ::AlphaLogic->SetInputConnection(1, ::MultiplyAlpha->GetOutputPort() );

  //      [ExtractRGB, AlphaLogic] -> AppendComponents
  ::AppendComponents->RemoveAllInputs();
  ::AppendComponents->AddInputConnection(0, ::ExtractRGB->GetOutputPort() );
  ::AppendComponents->AddInputConnection(0, ::AlphaLogic->GetOutputPort() );

  // ------------------------------------------
  // Pipeline - vtkMRMLSliceLayerLogic
  // ------------------------------------------
  ::Reslice->SetBackgroundColor(0, 0, 0, 0); // only first two are used
  ::Reslice->AutoCropOutputOff();
  ::Reslice->SetOptimization(1);
  ::Reslice->SetOutputOrigin( 0, 0, 0 );
  ::Reslice->SetOutputSpacing( 1, 1, 1 );
  ::Reslice->SetOutputDimensionality( 3 );

  // ------------------------------------------
  // Full Pipeline
  //      [[[InputImage]]] -> Reslice
  //      Reslice -> Threshold
  //      Reslice -> MapToWindowLevelColors
  //      MapToWindowLevelColors -> MapToColors
  //      MapToColors -> ExtractRGB
  //      MapToColors -> ExtractAlpha
  //      [ExtractAlpha, ResliceAlphaCast] -> MultiplyAlpha
  //      [Threshold, MultiplyAlpha] -> AlphaLogic
  //      [ExtractRGB, AlphaLogic] -> AppendComponents
  // ------------------------------------------

  // ------------------------------------------
  // Pipeline - vtkMRMLSliceLayerLogic - SetInput
  // ------------------------------------------
  ::Reslice->SetInterpolationModeToLinear();
  ::Reslice->SetInput(imageData);
  ::ResliceAlphaCast->SetInput(::Reslice->GetBackgroundMask());
  ::Reslice->GetBackgroundMask()->SetUpdateExtentToWholeExtent();

  vtkImageData * reslicedImage = ::Reslice->GetOutput();

  // ------------------------------------------
  // Pipeline - vtkMRMLScalarVolumeDisplayNode - SetInput
  // ------------------------------------------
  vtkNew<vtkMRMLColorTableNode> colorTableNode;
  colorTableNode->SetTypeToGrey();
  ::MapToColors->SetLookupTable(colorTableNode->GetLookupTable());
  ::Threshold->SetInput(reslicedImage);
  ::MapToWindowLevelColors->SetInput(reslicedImage);

  vtkImageData * foregroundImage = ::AppendComponents->GetOutput();

  // ------------------------------------------
  // Pipeline - vtkMRMLSliceLogic - SetInput
  // ------------------------------------------
  int layerIndex = 0;
  ::Blend->SetInput( layerIndex,  foregroundImage);
  ::Blend->SetOpacity( layerIndex++, 1.0 );

  widget.setImageData(::Blend->GetOutput());

  widget.show();

  return app.exec();
}

#include "moc_main.cpp"
