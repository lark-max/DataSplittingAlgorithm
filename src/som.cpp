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
using som::Inner_Cell;
using som::SOMProcess;

void SOMProcess::Train()
{
	Train(m_orderingRate, m_orderingNeiborSize, m_orderingEpochs);
	Train(m_tuningRate, m_tuningNeiborSize, m_tuningEpochs);

	// Cluster
	for (int i = 0; i < m_totalData; ++i)
		AllocateProcess(m_InputData[i], i);
}

void SOMProcess::Train(const double& learnRate, const int neighborSize, const int Epochs)
{
	int epoch = 0;
	int iteration = 0;
	m_maxIterations = Epochs * m_totalData;
	
	for (; epoch < Epochs; ++epoch)
	{
		// Randomise ordering of samples
		vector<int> order;
		method::Sequence(order, m_totalData);
		std::random_shuffle(order.begin(), order.end());

		// Train process of each epoch
		for (int i=0; i < m_totalData; ++i)
		{
			int index_ = order[i];
			double decayFunction = 1.0 - (double)iteration / (double)m_maxIterations;
			m_learningRate = decayFunction * learnRate;
			m_neighbourhoodSize = (int)(decayFunction * (double)neighborSize + 0.5);
			GetActiveNeuron(m_InputData[index_]);
			GetNeighbor();
			UpdateBias();
			UpdateWeight(m_InputData[index_]);
			++iteration;
		}
	}
}

void  SOMProcess::BeforeTrainProcess(matrix<double>& Data)
{
	m_InputData.swap(Data);
	m_totalData = m_InputData.rows();
	m_dimensions = m_InputData.columns();

	m_neuronNumber = (int)(2 * sqrt(m_totalData) + 0.5);
	m_neuronColumn = (int)(sqrt(m_neuronNumber / 1.6) + 0.5);
	m_neuronRow = (int)(1.6 * m_neuronColumn + 0.5);

	m_orderingRate = 0.8;
	m_orderingNeiborSize = (int)(0.5 * (double)m_neuronRow + 0.5);
	m_orderingEpochs = 3;

	m_tuningRate = 0.01;
	m_tuningNeiborSize = (int)(0.5 * (double)m_neuronColumn + 0.5);
	m_tuningEpochs = 4;

	m_consicience = 10.0;
	m_bias = 0.0001;

	SetNeuronSize();
	
	// Standardization or normalization should be carried out before SOM 
	// in order to prevent "unfair" input
	Standardised(m_InputData);
	
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
	double bias_ = 1.0 / m_neuronNumber;
	for (int i = 0; i < m_neuronRow; ++i)
		for (int j = 0; j < m_neuronColumn; ++j)
		{
			(*this)[i][j].m_bias = bias_;
			while (temp_++ < m_dimensions)
			{
				(*this)[i][j].m_weight.push_back((double)rand() / RAND_MAX - 0.5);
				//(*this)[i][j].m_weight.push_back(GaussRand());
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
	(*this).resize(m_neuronRow);
	m_ClusterSet.resize(m_neuronRow);
	for (int i = 0; i < m_neuronRow; ++i)
	{
		(*this)[i].resize(m_neuronColumn);
		m_ClusterSet[i].resize(m_neuronColumn);
	}

}

// By calculating Euclidean distance, 
// the index of winning neuron corresponding to each input layer data is obtained 
void	SOMProcess::GetActiveNeuron(const vector<double>& Data_I)
{
	m_activeRow = 0, m_activeColumn = 0;
	double Active_Distance = Getdistance_((*this)[0][0].m_weight, Data_I);
	for (int i = 0; i < m_neuronRow; ++i)
		for (int j = 0; j < m_neuronColumn; ++j)
		{
			double distance_ = Getdistance_((*this)[i][j].m_weight, Data_I) - m_consicience * (1.0 / m_neuronNumber - (*this)[i][j].m_bias);
			if (distance_ < Active_Distance)
			{
				m_activeRow = i;
				m_activeColumn = j;
				Active_Distance = distance_;
			}
		}
}

// Design neighborhood size 
void	SOMProcess::GetNeighbor()
{
	if (m_activeRow < m_neighbourhoodSize)
		m_neighbourhoodLowerRowIndex = 0;
	else
		m_neighbourhoodLowerRowIndex = m_activeRow - m_neighbourhoodSize;

	if (m_neuronRow - 1 < m_activeRow + m_neighbourhoodSize)
		m_neighbourhoodUpperRowIndex = m_neuronRow - 1;
	else
		m_neighbourhoodUpperRowIndex = m_activeRow + m_neighbourhoodSize;

	if (m_activeColumn < m_neighbourhoodSize)
		m_neighbourhoodLowerColumnIndex = 0;
	else
		m_neighbourhoodLowerColumnIndex = m_activeColumn - m_neighbourhoodSize;

	if (m_neuronColumn - 1 < m_activeColumn + m_neighbourhoodSize)
		m_neighbourhoodUpperColumnIndex = m_neuronColumn - 1;
	else
		m_neighbourhoodUpperColumnIndex = m_activeColumn + m_neighbourhoodSize;
}

void	SOMProcess::UpdateBias()
{
	for (int i = 0; i < m_neuronRow; ++i)
		for (int j = 0; j < m_neuronColumn; ++j)
		{
			if (i == m_activeRow && j == m_activeColumn)
				(*this)[i][j].m_bias += m_bias * (1 - (*this)[i][j].m_bias);
			else
				(*this)[i][j].m_bias -= m_bias * (*this)[i][j].m_bias;
		}			
}

// The closer to the winning neuron, the larger the adjustment range 
void	SOMProcess::UpdateWeight(const vector<double>& Data_I)
{
	for (int i = m_neighbourhoodLowerRowIndex; i < m_neighbourhoodUpperRowIndex; ++i)
		for (int j = m_neighbourhoodLowerColumnIndex; j < m_neighbourhoodUpperColumnIndex; ++j)
		{
			int dx = i - m_activeRow, dy = j - m_activeColumn;
			double d = sqrt((double)(dx * dx + dy * dy));
			double h = exp(-(d * d / (2.0 * m_neighbourhoodSize * m_neighbourhoodSize)));
			for (int i1 = 0; i1 < m_dimensions; ++i1)
				(*this)[i][j].m_weight[i1] += h * m_learningRate * (Data_I[i1] - (*this)[i][j].m_weight[i1]);
			//// The weights of neurons need to be normalized after adjustment
			//VectorNormalize((*this)[i][j].m_weight);	
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
	m_ClusterSet[m_activeRow][m_activeColumn].m_index.push_back(i);
}

void	SOMProcess::OutputClusterSet()
{
	std::ofstream outfile;
	outfile.open("cluster.txt");
	for (int i = 0; i < m_neuronRow; ++i)
	{
		for (int j = 0; j < m_neuronColumn; ++j)
			outfile << m_ClusterSet[i][j].m_index.size() << "\t";
		outfile << "\n";
	}
	outfile.close();
}