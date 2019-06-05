#include "stdafx.h"
#include <math.h>  
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkPNGImageIOFactory.h"
#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkExtractImageFilter.h"
#include "itkNiftiImageIOFactory.h"
#include "itkImageToVtkImageFilter.h"
#include "itkPNGImageIOFactory.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include "itkImageDuplicator.h"

typedef double PixelType;
typedef itk::Image< PixelType, 2 > ImageType;
typedef itk::Image< PixelType, 2 > OutputImageType;
typedef itk::Image< unsigned char, 2 > CastFilterType;
typedef itk::ImageDuplicator< ImageType > DuplicatorType;
typedef itk::ImageToVTKImageFilter< OutputImageType > ConnectorType;
using namespace std;

typedef itk::ExtractImageFilter< ImageType, OutputImageType > FilterType;
typedef itk::ImageFileReader< ImageType > ReaderType;

class parametros {
	public:
		string idIndividuo;
		double melhorQ;	
		double maxPorcentagem;
	friend std::ostream & operator << (std::ostream &out, const parametros & obj){
			out << obj.idIndividuo << "\n" << obj.melhorQ << "\n" << obj.maxPorcentagem << std::endl;
			return out;
	}
};

void qsigmoide(OutputImageType::Pointer image, int numSlice, double q, double a, double b, double threshold) {

	int x = image->GetLargestPossibleRegion().GetSize()[0];
	int y = image->GetLargestPossibleRegion().GetSize()[1];

	OutputImageType::SizeType regionSize;
	regionSize[0] = x;
	regionSize[1] = y;

	OutputImageType::IndexType regionIndex;
	regionIndex[0] = 0;
	regionIndex[1] = 0;

	OutputImageType::RegionType region;
	region.SetSize(regionSize);
	region.SetIndex(regionIndex);

	itk::ImageRegionIterator<OutputImageType> imageIterator1(image, region);
	double min = 0, max = 0;
	while (!(imageIterator1.IsAtEnd())) {
		if (imageIterator1.Get() < min) min = imageIterator1.Get();
		if (imageIterator1.Get() > max) max = imageIterator1.Get();
		++imageIterator1;
	}


	itk::ImageRegionIterator<OutputImageType> imageIterator3(image, region);
	while (!(imageIterator3.IsAtEnd())) {
		imageIterator3.Set((imageIterator3.Get() - min) / (max - min));
		++imageIterator3;
	}

	double expq = 0;
	int i = 0;
	int scale = 1;
	itk::ImageRegionIterator<OutputImageType> imageIterator(image, region);
	while (!(imageIterator.IsAtEnd())) {
		double D = abs((imageIterator.Get() - b) / a);

		if (q < 1.0) {

			expq = scale / (pow((1 + (1 - q)*D), (1 / (1 - q))));
		}
		else if (q == 1) {
			expq = scale / (1 + exp(-D));
		}
		else {
			D = -1 / D;
			expq = scale / (pow((1 + (1 - q) * D), (1 / (1 - q))));
		}
		if (expq >= numeric_limits<int>::max()) {
			expq = 1;
		}
		imageIterator.Set(expq);

		++imageIterator;
	}
	itk::ImageRegionIterator<OutputImageType> imageIterator2(image, region);
	while (!(imageIterator2.IsAtEnd())) {
		if (imageIterator2.Get() < threshold) {
			imageIterator2.Set(255);
		}
		else {
			imageIterator2.Set(0);
		}
		++imageIterator2;
	}
}


void readOriginal(ReaderType::Pointer reader, int i, string idIndividuo) {
	stringstream ss;
	ss << i;
	string numSlice;
	ss >> numSlice;
	reader->SetFileName("C:\\Users\\Windows 10Pro\\Desktop\\Emotion\\" + idIndividuo + "\\Original\\Original"
						+ numSlice + ".png");
	try {
		reader->Update();
	}
	catch (itk::ExceptionObject & error) {
		std::cerr << "Exception in reading file! " << std::endl;
		std::cerr << error << std::endl;
		return;
	}
	reader->UpdateOutputInformation();
}
void readActivation(ReaderType::Pointer reader, int i) {
	stringstream ss;
	ss << i;
	string numSlice;
	ss >> numSlice;
	reader->SetFileName("C:\\Users\\Windows 10Pro\\Desktop\\Emotion\\Ativação Emotion\\ativacao" + numSlice + ".png");
	try {
		reader->Update();
	}
	catch (itk::ExceptionObject & error) {
		std::cerr << "Exception in reading file! " << std::endl;
		std::cerr << error << std::endl;
		return;
	}

	int x = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
	int y = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
	ImageType::SizeType regionSize;
	regionSize[0] = x;
	regionSize[1] = y;

	ImageType::IndexType regionIndex;
	regionIndex[0] = 0;
	regionIndex[1] = 0;

	ImageType::RegionType region;
	region.SetSize(regionSize);
	region.SetIndex(regionIndex);
	itk::ImageRegionIterator<ImageType> activationIterator(reader->GetOutput(), region);
	while (!activationIterator.IsAtEnd()) {
		if (activationIterator.Get() > 0) {
			activationIterator.Set(255);
		}
		++activationIterator;
	}


}

double calculoErro(ImageType::Pointer sigmoidImage, ImageType::Pointer originalImage) {
	int x = sigmoidImage->GetLargestPossibleRegion().GetSize()[0];
	int y = sigmoidImage->GetLargestPossibleRegion().GetSize()[1];

	int numPixels = 0;

	ImageType::SizeType regionSize;
	regionSize[0] = x;
	regionSize[1] = y;

	ImageType::IndexType regionIndex;
	regionIndex[0] = 0;
	regionIndex[1] = 0;

	ImageType::RegionType region;
	region.SetSize(regionSize);
	region.SetIndex(regionIndex);

	itk::ImageRegionIterator<ImageType> sigmoidIterator(sigmoidImage, region);
	itk::ImageRegionIterator<ImageType> originalIterator(originalImage, region);
	double brancosMascara = 0, pretosMascara = 0, brancosMatch = 0, pretosMatch = 0;
	while (!originalIterator.IsAtEnd()) {
		if (originalIterator.Get() > 0) {
			brancosMascara++;
			if (sigmoidIterator.Get() > 0) {
				brancosMatch++;
			}
		}
		else {
			pretosMascara++;
			if (sigmoidIterator.Get() == 0) {
				pretosMatch++;
			}
		}
		++originalIterator;
		++sigmoidIterator;
	}

	double porcentagemBrancos = brancosMatch / brancosMascara;
	if (brancosMascara == 0) porcentagemBrancos = 0;
	double porcentagemPretos = pretosMatch / pretosMascara;
	return ((porcentagemBrancos + porcentagemPretos) / 2);
}
vector<string> lerArquivoIndividuos() {
	vector<string> individuos;
	ifstream arquivoIndividuos("C:\\Users\\Windows 10Pro\\Desktop\\IndividuosRestantes.txt");
	string strIndividuo;
	while (std::getline(arquivoIndividuos, strIndividuo))
	{
		// Line contains string of length > 0 then save it in vector
		if (strIndividuo.size() > 0)
			individuos.push_back(strIndividuo);
	}
	return individuos;
}

vector<parametros> lerArquivoParametros() {
	vector<parametros> parametrosSalvos;
	ifstream arquivoParametros("C:\\Users\\Windows 10Pro\\Desktop\\ParametrosOtimizadosEmotion2.txt");
	int idIndividuo;
	double melhorQ, maxPorcentagem;
	while (arquivoParametros >> idIndividuo >> melhorQ >> maxPorcentagem)
	{
		parametros parametroIndividual;
		parametroIndividual.melhorQ = melhorQ;
		parametroIndividual.maxPorcentagem = maxPorcentagem;
		parametroIndividual.idIndividuo = idIndividuo;
		// Line contains string of length > 0 then save it in vector
		parametrosSalvos.push_back(parametroIndividual);
	}
	return parametrosSalvos;
}
int main(int argc, char ** argv) {
	itk::PNGImageIOFactory::RegisterOneFactory();
	vector<string> individuos = lerArquivoIndividuos();
	vector<parametros> parametrosSalvos = lerArquivoParametros();

	
	vector<ReaderType::Pointer> activationReaderVector;
	double z = 91;
	double somaPorcentagemAcerto = 0;
	double porcentagemAcerto = 0;
	itk::PNGImageIOFactory::RegisterOneFactory();
	for (int i = 0; i <= 90; i++) {
		ReaderType::Pointer activationReader = ReaderType::New();
		readActivation(activationReader, i);
		activationReaderVector.push_back(activationReader);
	}

	for (int numIndividuo = 0; numIndividuo < individuos.size(); numIndividuo++) {
		cout << "individuo: " << individuos[numIndividuo] << endl;
		vector<ReaderType::Pointer> originalReaderVector;
		for (int i = 0; i <= 90; i++) {
			ReaderType::Pointer originalReader = ReaderType::New();
			readOriginal(originalReader, i, individuos[numIndividuo]);
			originalReaderVector.push_back(originalReader);
		}
		double maxPorcentagem = 0, maiorQ = 0, q = 0;
		double b = .2, a = .1, t = 0.15;
		for (q = 0.01; q <= 2; q += 0.01) {

			if (q == 1.0000000000000006) {
				q = q + 0.01;
			}
			cout << "q = " << q << endl;
			cout << endl;
			for (int i = 0; i <= 90; i++) {
				ImageType::Pointer originalImage;
				DuplicatorType::Pointer duplicator = DuplicatorType::New();
				duplicator->SetInputImage(originalReaderVector[i]->GetOutput());
				duplicator->Update();
				originalImage = duplicator->GetOutput();
				ImageType::Pointer activationImage = activationReaderVector[i]->GetOutput();
				qsigmoide(originalImage, 0, q, a, b, t);
				somaPorcentagemAcerto += calculoErro(originalImage, activationImage);
			}
			porcentagemAcerto = (somaPorcentagemAcerto / z) * 100;
			
			if (porcentagemAcerto > maxPorcentagem) {
				maxPorcentagem = porcentagemAcerto;
				maiorQ = q;
			}

			porcentagemAcerto = 0;
			somaPorcentagemAcerto = 0;
		}


		cout << endl << endl << endl;
		cout << "*************MELHOR VALORES*************" << endl;
		cout << "individuo: " << individuos[numIndividuo] << endl;
		cout << "q: " << maiorQ << endl;
		cout << "porcentagem: " << maxPorcentagem << endl;
		cout << endl;
		
		parametros novoIndividuo;
		novoIndividuo.melhorQ = maiorQ;
		novoIndividuo.maxPorcentagem = maxPorcentagem;
		novoIndividuo.idIndividuo = individuos[numIndividuo];
		parametrosSalvos.push_back(novoIndividuo);
		ofstream saida("C:\\Users\\Windows 10Pro\\Desktop\\ParametrosOtimizadosEmotion2.txt");
		for (int qtdIndividuosSalvos = 0; qtdIndividuosSalvos < parametrosSalvos.size(); qtdIndividuosSalvos++) {
			saida << parametrosSalvos[qtdIndividuosSalvos];
		}
		saida.close();
		maxPorcentagem = 0;
		maiorQ = 0;

	}
	return EXIT_SUCCESS;
}
