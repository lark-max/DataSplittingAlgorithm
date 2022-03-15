// SOM algorithm
// Written by Chen J. Y. at 2020, Zhejiang UNV
#include<iostream>
#include<fstream>
#include<algorithm>
#include<ctime>
#include<cmath>
#include<numeric>
#include"som.h"
using som::matrix;
using som::Inner_struct;
using som::SOMProcess;

void SOMProcess::Train()
{
	Train(orderingRate, orderingNeiborSize, orderingEpochs);
	Train(tuningRate, tuningNeiborSize, tuningEpochs);

	// Cluster
	for (int i = 0; i < totalData; ++i)
		AllocateProcess(InputData[i], i);
}

void SOMProcess::Train(const double& learningRate, const int neighborSize, const int Epochs)
{
	int epoch = 0;
	int iteration = 0;
	maxIterations = Epochs * totalData;
	
	for (; epoch < Epochs; ++epoch)
	{
		// Randomise ordering of samples
		vector<int> order;
		method::Sequence(order, totalData);
		std::random_shuffle(order.begin(), order.end());

		// Train process of each epoch
		for (int i=0; i < totalData; ++i)
		{
			int index_ = order[i];
			double decayFunction = 1.0 - (double)iteration / (double)maxIterations;
			learningRate_ = decayFunction * learningRate;
			neighbourhoodSize_ = (int)(decayFunction * (double)neighborSize + 0.5);
			GetActiveNeuron(InputData[index_]);
			GetNeighbor();
			UpdateBias();
			UpdateWeight(InputData[index_]);
			++iteration;
		}
	}
}

void  SOMProcess::BeforeTrainProcess(matrix<double>& Data)
{
	InputData = Data;
	totalData = InputData.rows();
	dimensions_ = InputData.columns();

	neuronNumber_ = (int)(2 * sqrt(totalData) + 0.5);
	neuronColumn = (int)(sqrt(neuronNumber_ / 1.6) + 0.5);
	neuronRow = (int)(1.6 * neuronColumn + 0.5);

	orderingRate = 0.8;
	orderingNeiborSize = (int)(0.5 * (double)neuronRow + 0.5);
	orderingEpochs = 3;

	tuningRate = 0.01;
	tuningNeiborSize = (int)(0.5 * (double)neuronColumn + 0.5);
	tuningEpochs = 4;

	consicience = 10.0;
	bias = 0.0001;

	SetNeuronSize();
	
	// Standardization or normalization should be carried out before SOM 
	// in order to prevent "unfair" input
	Standardised(InputData);
	
	// Initially setting weights
	RandomlyInitialize();		
}

double SOMProcess::GetMean(const vector<double>& _vector)
{
	return accumulate(_vector.begin(), _vector.end(), 0.0) / _vector.size();
}

double SOMProcess::GetStdev(const vector<double>& _vector)
{
	double mean_ = GetMean(_vector);
	double var_ = 0.0;
	for (vector<double>::const_iterator iter = _vector.begin();
		iter != _vector.end(); ++iter)
	{
		var_ += pow(*iter - mean_, 2.0);
	}
	var_ /= _vector.size();
	return sqrt(var_);
}

double	SOMProcess::GaussRand(double E, double V)
{
	static	double V1, V2, s;
	static int phase = 0;
	double	x;

	if (phase == 0)
	{
		do
		{
			double	U1 = (double)rand() / RAND_MAX;
			double	U2 = (double)rand() / RAND_MAX;
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			s = V1 * V1 + V2 * V2;
		} while (s >= 1 || s == 0);
		x = V1 * sqrt(-2 * log(s) / s);
	}
	else
		x = V2 * sqrt(-2 * log(s) / s);

	phase = 1 - phase;
	return	(x * V + E);
}

// data standardised for one dimension data(only output)
void	SOMProcess::Standardised(matrix<double>& data)
{
	//if (data.columns() == 1)
	//	return;

	vector<double> tempvec(data.rows());
	for (int j = 0; j < data.columns(); ++j)
	{
		for (int i = 0; i < tempvec.size(); ++i)
			tempvec[i] = data[i][j];
		double mean = GetMean(tempvec);
		double stdev = GetStdev(tempvec);
		for (int i = 0; i < data.rows(); ++i)
			data[i][j] = (data[i][j] - mean) / stdev;
	}
}

void	SOMProcess::Normalize(matrix<double>& data)
{
	for (int j = 0; j < data.columns(); ++j)
	{
		double max = data.Max(j);
		for (int i = 0; i < data.rows(); ++i)
			data[i][j] /= max;
	}
}

// The weights are randomly initialized (multdimensional)
void	SOMProcess::RandomlyInitialize()
{
	int temp_ = 0;
	double bias_ = 1.0 / neuronNumber_;
	for (int i = 0; i < neuronRow; ++i)
		for (int j = 0; j < neuronColumn; ++j)
		{
			(*this)[i][j]._bias = bias_;
			while (temp_++ < dimensions_)
			{
				(*this)[i][j].W.push_back((double)rand() / RAND_MAX - 0.5);
				//(*this)[i][j].W.push_back(GaussRand());
			}
				
			temp_ = 0;
		}
}


//void	SOMProcess::VectorNormalize(vector<double>& _vector)
//{
//	if (_vector.size() < 2)
//		return;
//	double sum_ = accumulate(_vector.begin(), _vector.end(), 0.0);
//	for (vector<double>::iterator it = _vector.begin(); it != _vector.end(); ++it)
//		*it /= sum_;
//}

void	SOMProcess::SetNeuronSize()
{
	(*this).resize(neuronRow);
	ClusterSet.resize(neuronRow);
	for (int i = 0; i < neuronRow; ++i)
	{
		(*this)[i].resize(neuronColumn);
		ClusterSet[i].resize(neuronColumn);
	}

}

// By calculating Euclidean distance, 
// the index of winning neuron corresponding to each input layer data is obtained 
void	SOMProcess::GetActiveNeuron(const vector<double>& Data_I)
{
	activeRow = 0, activeColumn = 0;
	double Active_Distance = Getdistance_((*this)[0][0].W, Data_I);
	for (int i = 0; i < neuronRow; ++i)
		for (int j = 0; j < neuronColumn; ++j)
		{
			double distance_ = Getdistance_((*this)[i][j].W, Data_I) - consicience * (1.0 / neuronNumber_ - (*this)[i][j]._bias);
			if (distance_ < Active_Distance)
			{
				activeRow = i;
				activeColumn = j;
				Active_Distance = distance_;
			}
		}
}

// Design neighborhood size 
void	SOMProcess::GetNeighbor()
{
	if (activeRow < neighbourhoodSize_)
		neighbourhoodLowerRowIndex_ = 0;
	else
		neighbourhoodLowerRowIndex_ = activeRow - neighbourhoodSize_;

	if (neuronRow - 1 < activeRow + neighbourhoodSize_)
		neighbourhoodUpperRowIndex_ = neuronRow - 1;
	else
		neighbourhoodUpperRowIndex_ = activeRow + neighbourhoodSize_;

	if (activeColumn < neighbourhoodSize_)
		neighbourhoodLowerColumnIndex_ = 0;
	else
		neighbourhoodLowerColumnIndex_ = activeColumn - neighbourhoodSize_;

	if (neuronColumn - 1 < activeColumn + neighbourhoodSize_)
		neighbourhoodUpperColumnIndex_ = neuronColumn - 1;
	else
		neighbourhoodUpperColumnIndex_ = activeColumn + neighbourhoodSize_;
}

void	SOMProcess::UpdateBias()
{
	for (int i = 0; i < neuronRow; ++i)
		for (int j = 0; j < neuronColumn; ++j)
		{
			if (i == activeRow && j == activeColumn)
				(*this)[i][j]._bias += bias * (1 - (*this)[i][j]._bias);
			else
				(*this)[i][j]._bias -= bias * (*this)[i][j]._bias;
		}			
}

// The closer to the winning neuron, the larger the adjustment range 
void	SOMProcess::UpdateWeight(const vector<double>& Data_I)
{
	for (int i = neighbourhoodLowerRowIndex_; i < neighbourhoodUpperRowIndex_; ++i)
		for (int j = neighbourhoodLowerColumnIndex_; j < neighbourhoodUpperColumnIndex_; ++j)
		{
			int dx = i - activeRow, dy = j - activeColumn;
			double d = sqrt((double)(dx * dx + dy * dy));
			double h = exp(-(d * d / (2.0 * neighbourhoodSize_ * neighbourhoodSize_)));
			for (int i1 = 0; i1 < dimensions_; ++i1)
				(*this)[i][j].W[i1] += h * learningRate_ * (Data_I[i1] - (*this)[i][j].W[i1]);
			//// The weights of neurons need to be normalized after adjustment
			//VectorNormalize((*this)[i][j].W);	
		}
}

// Calculating Euclidean distance 
double	SOMProcess::Getdistance_(const vector<double>& m, const vector<double>& n)
{
	double _euclidean_distance_ = 0;
	if (m.size() != n.size())
	{
		std::cout << "Fatal error: Dimenson mismatch!" << std::endl;
		system("pause");
		exit(1);	// error
	}

	else
	{
		for (int i = 0; i < m.size(); ++i)
			_euclidean_distance_ += pow(m[i] - n[i], 2.0);
		return sqrt(_euclidean_distance_);
	}
}

// Gather the data points of input layer into the trained network 
void	SOMProcess::AllocateProcess(const vector<double>& Data_I, int& i)
{
	GetActiveNeuron(Data_I);
	ClusterSet[activeRow][activeColumn]._index.push_back(i);
}

void	SOMProcess::OutputClusterSet()
{
	std::ofstream outfile;
	outfile.open("cluster.txt");
	for (int i = 0; i < neuronRow; ++i)
	{
		for (int j = 0; j < neuronColumn; ++j)
			outfile << ClusterSet[i][j]._index.size() << "\t";
		outfile << "\n";
	}
	outfile.close();
}