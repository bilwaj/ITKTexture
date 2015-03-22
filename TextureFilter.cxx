#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const int MAX_VAL=100000;
typedef short datatype;
typedef itk::Image<datatype, 3> ImageType;
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileWriter<ImageType> WriterType;
typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;


void Tex(ImageType::Pointer Inp, ImageType::Pointer Tex)
{ 
	IteratorType InpIt(Inp,Inp->GetRequestedRegion());
	IteratorType TexIt(Tex,Tex->GetRequestedRegion());
	//Setting up the neighborhood iterator
	 ImageType::SizeType radius;
     radius[0] = 1;
     radius[1] = 1;
     radius[2] = 1;
	 itk::NeighborhoodIterator<ImageType> InpNIt(radius, Inp, Inp->GetRequestedRegion());

	InpNIt.GoToBegin();
	TexIt.GoToBegin();
	while(!TexIt.IsAtEnd())
	{ 
	  TexIt.Set(0);
	  ++TexIt;
      ++InpNIt;
    }
}



//Expects smart pointers to scalar image objects. These need to be pointing to 2 separate memory blocks
void MeanTexture(ImageType::Pointer Inp, ImageType::Pointer MeanTexture)
{ 
	IteratorType TexIt(MeanTexture,MeanTexture->GetRequestedRegion());
	typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIterator;
	NeighborhoodIterator::RadiusType radius;
	for(int i=0;i<ImageType::ImageDimension;i++)
	{radius[i]=1;}
	NeighborhoodIterator InpIt(radius,Inp,Inp->GetRequestedRegion());
	TexIt.GoToBegin();
	while(!TexIt.IsAtEnd())
	{
		datatype accum=0.0;
		for(int j=0;j<InpIt.Size();j++)
		{
			accum=accum+InpIt.GetPixel(j);
		}
		TexIt.Set((datatype)accum/InpIt.Size());
	    TexIt++;
		InpIt++;
	}
 	std::cout<<"Traversing done"<<std::endl;
}

//Returns the estimated moment around the mean associated of the degree(th) order
//TO DO Ext defaults to 1; Ideally should be part of moments class
void MomentTexture(ImageType::Pointer Inp, ImageType::Pointer MomentTexture, int degree, int Ext=1)
{ 
	IteratorType TexIt(MomentTexture,MomentTexture->GetRequestedRegion());
	typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIterator;
	NeighborhoodIterator::RadiusType radius;
	for(int i=0;i<ImageType::ImageDimension;i++)
	{radius[i]=1;}
	NeighborhoodIterator InpIt(radius,Inp,Inp->GetRequestedRegion());
	TexIt.GoToBegin();
	while(!TexIt.IsAtEnd())
	{
		datatype accum=0.0;
		datatype min=0; //Store minimum voxel intensity value
		datatype max=0; //Store maximum voxel intensity value
	
		//first compute the mean
		//TO DO use mean texture function to not recompute the mean every time
		for(int j=0;j<InpIt.Size();j++)
		{
			accum=accum+InpIt.GetPixel(j);
			min=MIN(min,InpIt.GetPixel(j));
			max=MAX(max,InpIt.GetPixel(j));
		}
		datatype range=max-min;
		float mean=accum/InpIt.Size();
		double accumX=0.0;
		for(int j=0;j<InpIt.Size();j++)
		{   
			//Modified expression for moments to keep them range bound and manageable
			accumX=accumX+std::pow((InpIt.GetPixel(j)-mean)/range,degree);
		}
		//100 is usedas a multiplier to keep images produced to be visualized
		TexIt.Set((datatype)(1000*accumX/InpIt.Size()));
	    TexIt++;
		InpIt++;
	}
 	std::cout<<"Traversing done"<<std::endl;
}


int main(int argc, char *argv[])
{
  ReaderType::Pointer input = ReaderType::New();
  //Input file
  input->SetFileName("1.nii");
  input->Update();

  ImageType::Pointer Inp=ImageType::New();
  Inp=input->GetOutput();
  
  ImageType::Pointer Tex = ImageType::New();;
  Tex->CopyInformation(Inp);
  Tex->SetRequestedRegion( Inp->GetRequestedRegion() );
  Tex->SetBufferedRegion( Inp->GetBufferedRegion() );
  Tex->Allocate();
  Tex->FillBuffer(0);

  
  //Second moment neighborhood radius of 3
  MomentTexture(Inp,Tex,2,3);
  //The function is generic and multiple runs ,ay be used to create feature vectors as desired Inp and Tex are image types

  
  WriterType::Pointer output = WriterType::New();
  //Output file
  output->SetFileName("out.nii");
  output->SetInput(Tex);
  output->Update();

  //std::cout<<max_val<<std::endl;
  return EXIT_SUCCESS;
}
