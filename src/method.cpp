// This file includes several methods for data splitting
// Written by Chen J. Y. at 2020, Zhejiang UNV
#include<iostream>
#include<iomanip>
#include<fstream>
#include<numeric>
#include<algorithm>
#include<cmath>
#include<memory>
#include"method.h"
using method::Method;
using som::SOMProcess;
using som::matrix;

Method::Method(const double& a, const double& b, matrix<double>& data,matrix<double>& iniData)
	:m_Ratio_train(a),m_Ratio_test(b), m_data(data), m_Initial_Data(iniData)
{}

// Used to import data from the .txt file under the current path 
void	method::GetTXTData(matrix<double>& data,string& str,vector<int>& Idex)
{
	vector< vector<double> > temp_data;
	std::ifstream	inpFile;
	std::string  index_;
	double temp_input;

	inpFile.open(str);
	if (!inpFile)
	{
		std::ofstream outCerr;
		outCerr.open("log.out", std::ios::app);
		outCerr << "There is no such file named " << str << "\n"
			<< "Please make sure the file is in the current folder!" << std::endl;
		outCerr.close();
		exit(1);
	}
	std::getline(inpFile, index_);
	int column_size = std::count(index_.begin(), index_.end(), '\t') + 1;
	int temp_ = 0;
	temp_data.push_back(vector<double>());
	while (inpFile >> temp_input)
	{
		if (inpFile.eof())
			break;
		if (temp_ % column_size == 0 && temp_ / column_size != 0)
			temp_data.push_back(vector<double>());
		temp_data.back().push_back(temp_input);
		++temp_;
	}
	inpFile.close();
	int row = (int)temp_data.size();
	int col = (int)temp_data[0].size() - 1;
	data.resize(row, vector<double>(col));
	for (int i = 0; i < row; ++i)
	{
		Idex.push_back(round(temp_data[i][0]));
		for (int j = 0; j < col; ++j)
			data[i][j] = temp_data[i][j + 1];
	}
}

void	method::Sequence(vector<int>& Idex, int max, int min)
{
	for (int i = min; i < max; ++i)
		Idex.push_back(i);
}

void	Method::SetVector()
{
	m_NumberOfData_eachNeuron.resize(m_afterSOM->m_neuronRow, vector<int>(m_afterSOM->m_neuronColumn));
	m_TrainNumber_eachNeuron.resize(m_afterSOM->m_neuronRow, vector<int>(m_afterSOM->m_neuronColumn));
	m_TestNumber_eachNeuron.resize(m_afterSOM->m_neuronRow, vector<int>(m_afterSOM->m_neuronColumn));
	m_ValidNumber_EachNeuron.resize(m_afterSOM->m_neuronRow, vector<int>(m_afterSOM->m_neuronColumn));
	m_theta_numerator.resize(m_afterSOM->m_neuronRow, vector<double>(m_afterSOM->m_neuronColumn));
}

// The number of data in each stratum 
void	Method::GetNumberOfData_eachNeuron()
{
	for (int i = 0; i < m_afterSOM->m_neuronRow; ++i)
		for (int j = 0; j < m_afterSOM->m_neuronColumn; ++j)
			m_NumberOfData_eachNeuron[i][j] = m_afterSOM->m_ClusterSet[i][j].m_index.size();
}

// The amount of data allocated to train, test, and validate at each stratum
// According to equal proportion 
void	Method::GetSampleNumber_eachNeuron_proportion()
{
	for (int i = 0; i < m_afterSOM->m_neuronRow; ++i)
		for (int j = 0; j < m_afterSOM->m_neuronColumn; ++j)
		{
			m_TrainNumber_eachNeuron[i][j] = int(m_Ratio_train / (m_Ratio_train + m_Ratio_test)
				* m_NumberOfData_eachNeuron[i][j] * (m_Ratio_train + m_Ratio_test) + 0.5);
			m_TestNumber_eachNeuron[i][j] = int(m_Ratio_test / (m_Ratio_train + m_Ratio_test)
				* m_NumberOfData_eachNeuron[i][j] * (m_Ratio_train + m_Ratio_test) + 0.5);
			m_ValidNumber_EachNeuron[i][j] = m_NumberOfData_eachNeuron[i][j] - m_TrainNumber_eachNeuron[i][j]
				- m_TestNumber_eachNeuron[i][j];
		}
}

void	Method::GetSampleNumber_eachNeuron_neyman()
{
	for (int i = 0; i < m_afterSOM->m_neuronRow; ++i)
		for (int j = 0; j < m_afterSOM->m_neuronColumn; ++j)
		{
			m_TrainNumber_eachNeuron[i][j] = int(m_TrainNumber * m_theta_numerator[i][j] / m_theta_denominator + 0.5);
			m_TestNumber_eachNeuron[i][j] = int(m_TestNumber * m_theta_numerator[i][j] / m_theta_denominator + 0.5);

			// Prevent the number of samples from being greater than the number of data in this stratum
			// Sometimes Neyman allocation has this problem!
			if (m_TrainNumber_eachNeuron[i][j] > m_NumberOfData_eachNeuron[i][j])
			{
				m_TrainNumber_eachNeuron[i][j] = m_NumberOfData_eachNeuron[i][j];
				m_TestNumber_eachNeuron[i][j] = 0;
			}
			else if (m_TrainNumber_eachNeuron[i][j] + m_TestNumber_eachNeuron[i][j] > m_NumberOfData_eachNeuron[i][j])
				m_TestNumber_eachNeuron[i][j] = m_NumberOfData_eachNeuron[i][j] - m_TrainNumber_eachNeuron[i][j];
			else
			{
				;	// do nothing
			}			

		}
}

// Random sampling in each stratum 
void	Method::Random_sample()
{
	vector<int>::iterator it;
	for (int i = 0; i < m_afterSOM->m_neuronRow; ++i)
		for (int j = 0; j < m_afterSOM->m_neuronColumn; ++j)
		{
			int tmp = 0;
			int sampleTimes = 0;
			// m_trainSet sampling
			if (m_NumberOfData_eachNeuron[i][j] != 0 && m_TrainNumber_eachNeuron[i][j] != 0)
			{
				int trainSet_initial_index = 0;
				while (sampleTimes < m_TrainNumber_eachNeuron[i][j])
				{
					trainSet_initial_index = m_trainSet.size();
					if (trainSet_initial_index >= m_TrainNumber)
						return;
					tmp = rand() % m_afterSOM->m_ClusterSet[i][j].m_index.size();
					m_trainSet.push_back(m_afterSOM->m_ClusterSet[i][j].m_index[tmp]);
					it = m_afterSOM->m_ClusterSet[i][j].m_index.begin() + tmp;
					m_afterSOM->m_ClusterSet[i][j].m_index.erase(it);
					++sampleTimes;
				}
			}
			// m_testSet sampling
			if (m_NumberOfData_eachNeuron[i][j] != 0 && m_TestNumber_eachNeuron[i][j] != 0)
			{
				int testSet_initial_index = 0;
				sampleTimes = 0;
				while (sampleTimes < m_TestNumber_eachNeuron[i][j])
				{
					testSet_initial_index = m_testSet.size();
					if (testSet_initial_index >= m_TestNumber)
						return;
					tmp = rand() % m_afterSOM->m_ClusterSet[i][j].m_index.size();
					m_testSet.push_back(m_afterSOM->m_ClusterSet[i][j].m_index[tmp]);
					it = m_afterSOM->m_ClusterSet[i][j].m_index.begin() + tmp;
					m_afterSOM->m_ClusterSet[i][j].m_index.erase(it);
					++sampleTimes;
				}
			}
		}
}

void	Method::GetTheta()
{
	m_theta_denominator = 0;
	double temp;
	double average;
	vector<double> var(m_data[0].size());


	vector<int>::iterator iter;
	for (int i = 0; i < m_afterSOM->m_neuronRow; ++i)
		for (int j = 0; j < m_afterSOM->m_neuronColumn; ++j)
		{
			if (m_NumberOfData_eachNeuron[i][j] == 0)
			{
				m_theta_numerator[i][j] = 0;
				continue;
			}
			// Calculation formula of standard deviation for each stratum:
			// theta = sqrt(theta_x1^2+theta_x2^2+theta_x3^2+...+theta_xp^2+theta_y^2), p = dimension of data
			// where theta_xi^2 is the intra-stratum variance of component xi,
			// and theta_y^2 is the variance of the output variable, y.
			for (int k = 0; k < m_data.columns(); ++k)
			{
				temp = 0;
				average = 0;
				var[k] = 0;

				for (iter = m_afterSOM->m_ClusterSet[i][j].m_index.begin();
					iter != m_afterSOM->m_ClusterSet[i][j].m_index.end(); ++iter)
					temp += m_data[*iter][k];

				average = temp / m_NumberOfData_eachNeuron[i][j];

				for (iter = m_afterSOM->m_ClusterSet[i][j].m_index.begin();
					iter != m_afterSOM->m_ClusterSet[i][j].m_index.end(); ++iter)
					var[k] += pow(m_data[*iter][k] - average, 2)
					/ m_NumberOfData_eachNeuron[i][j];
			}
			double _temp = accumulate(var.begin(), var.end(), 0.0);
			m_theta_numerator[i][j] = m_NumberOfData_eachNeuron[i][j] * sqrt(_temp);
			m_theta_denominator += m_theta_numerator[i][j];
		}

}

// After sampling, the remaining data of each stratum is allocated to the valid set 
void	Method::AllocateRemainDataInto_validSet()
{
	vector<int>::iterator it;
	for (int i = 0; i < m_afterSOM->m_neuronRow; ++i)
		for (int j = 0; j < m_afterSOM->m_neuronColumn; ++j)
			for (it = m_afterSOM->m_ClusterSet[i][j].m_index.begin(); it != m_afterSOM->m_ClusterSet[i][j].m_index.end(); ++it)
				m_validSet.push_back(*it);

}

void	Method::GetThisColumnData(vector<double>& _column_,int& j, const vector<int>& CurrentDataIndex)
{
	for (int i = 0; i < _column_.size(); ++i)
		_column_[i] = m_data[CurrentDataIndex[i]][j];
}

// Calculate the mean and stdev  of each column
// Then datum = (datum - mean)/stdev;
void	Method::DP_Standardise(const vector<int>& CurrentDataIndex)
{
	if (m_data.columns()<2 || CurrentDataIndex.size() < 2)
		return;
	vector<double> _column_(CurrentDataIndex.size());
	for (int j = 0; j < m_data.columns(); ++j)
	{
		GetThisColumnData(_column_, j, CurrentDataIndex);
		double mean_ = SOMProcess::GetMean(_column_);
		double stdev_ =SOMProcess::GetStdev(_column_);
		for (int i = 0; i < CurrentDataIndex.size(); ++i)
			m_data[CurrentDataIndex[i]][j] = (m_data[CurrentDataIndex[i]][j] - mean_) / stdev_;
	}
}

// The first step of DUPLEX method: find two data points with the largest Euclidean distance 
void	Method::DP_InitialSample(vector<int>& ini_Data, vector<int>& Key)
{
	if (ini_Data.size() < 2)
		return;
	vector<int>::iterator a = ini_Data.begin();
	vector<int>::iterator b = ini_Data.begin() + 1;
	double	maxDistance = 0;
	for (vector<int>::iterator t = ini_Data.begin(); t != ini_Data.end() - 1; ++t)
		for (vector<int>::iterator tt = t + 1; tt != ini_Data.end(); ++tt)
		{
			//double _distance = m_Distance_Matrix[*t][*tt];
			//double _distance = abs(m_data[*t][0] - m_data[*tt][0]);
			double _distance = Distance_matrix(*t, *tt);
			if (_distance > maxDistance)
			{
				maxDistance = _distance;
				a = t;
				b = tt;
			}
		}
	Key.push_back(*a);
	Key.push_back(*b);
	ini_Data.erase(b);
	ini_Data.erase(a);

}

// Find the next sampled pair xi and xj, such that they maximize the minimum single-linkage distance ||x − s||
// for each previously sampled s ∈ T and allocate to T.
void	Method::DP_Resample(vector<int>& ini_Data, vector<int>& Key)
{
	if (ini_Data.size() <= 2)
	{
		while (ini_Data.size() > 0)
		{
			Key.push_back(ini_Data.back());
			ini_Data.pop_back();
		}
		return;
	}

	vector<int>::iterator a = ini_Data.begin();
	vector<int>::iterator b = ini_Data.begin() + 1;
	double distA = 0;
	double distB = 0;

	for (vector<int>::iterator t = ini_Data.begin(); t != ini_Data.end(); ++t)
	{
		double _distance = 1.0e+6;		
		for (vector<int>::iterator s = Key.begin(); s != Key.end(); ++s)
		{
			double temp_dist = Distance_matrix(*t, *s);
			if (temp_dist < _distance)
				_distance = temp_dist;
		}

		// First update a, then update b
		// a points to the farthest node, b points to the second farthest node
		if (_distance >= distA)
		{
			b = a;
			distB = distA;
			a = t;
			distA = _distance;
		}
		else if (_distance >= distB)
		{
			b = t;
			distB = _distance;
		}
		else
			continue;
	}
	Key.push_back(*a);
	Key.push_back(*b);

	if (a > b)
	{
		ini_Data.erase(a);
		ini_Data.erase(b);
	}
	else
	{
		ini_Data.erase(b);
		ini_Data.erase(a);
	}

}

// In each stratum, use the DUPLEX method to allocate data
void	Method::BasedSOM_DPSample()
{
	for (int i = 0; i < m_afterSOM->m_neuronRow; ++i)
		for (int j = 0; j < m_afterSOM->m_neuronColumn; ++j)
		{
			m_trainKey.clear();
			m_testKey.clear();
			m_validKey.clear();

			if (m_TrainNumber_eachNeuron[i][j] > 0)
				DP_InitialSample(m_afterSOM->m_ClusterSet[i][j].m_index, m_trainKey);
			if (m_TestNumber_eachNeuron[i][j] > 0)
				DP_InitialSample(m_afterSOM->m_ClusterSet[i][j].m_index, m_testKey);
			if (m_ValidNumber_EachNeuron[i][j] > 0)
				DP_InitialSample(m_afterSOM->m_ClusterSet[i][j].m_index, m_validKey);
			while (m_afterSOM->m_ClusterSet[i][j].m_index.size() > 0)
			{
				if (m_trainKey.size() < m_TrainNumber_eachNeuron[i][j])
					DP_Resample(m_afterSOM->m_ClusterSet[i][j].m_index, m_trainKey);
				if (m_testKey.size() < m_TestNumber_eachNeuron[i][j])
					DP_Resample(m_afterSOM->m_ClusterSet[i][j].m_index, m_testKey);
				if (m_validKey.size() < m_ValidNumber_EachNeuron[i][j])
					DP_Resample(m_afterSOM->m_ClusterSet[i][j].m_index, m_validKey);

				if (CheckFull(m_afterSOM->m_ClusterSet[i][j].m_index, m_trainKey, m_testKey, m_validKey,
					m_TrainNumber_eachNeuron[i][j], m_TestNumber_eachNeuron[i][j], m_ValidNumber_EachNeuron[i][j]))
					break;
			}
			std::copy(m_trainKey.begin(), m_trainKey.end(), std::back_inserter(m_trainSet));
			std::copy(m_testKey.begin(), m_testKey.end(), std::back_inserter(m_testSet));
			std::copy(m_validKey.begin(), m_validKey.end(), std::back_inserter(m_validSet));
		}
}


void	Method::BasedSOM_MDP_Sample()
{
	for (int i = 0; i < m_afterSOM->m_neuronRow; ++i)
		for (int j = 0; j < m_afterSOM->m_neuronColumn; ++j)
		{
			m_trainKey.clear();
			m_testKey.clear();
			m_validKey.clear();

			if (m_TrainNumber_eachNeuron[i][j] > 0)
				DP_InitialSample(m_afterSOM->m_ClusterSet[i][j].m_index, m_trainKey);
			if (m_TestNumber_eachNeuron[i][j] > 0)
				DP_InitialSample(m_afterSOM->m_ClusterSet[i][j].m_index, m_testKey);
			if (m_ValidNumber_EachNeuron[i][j] > 0)
				DP_InitialSample(m_afterSOM->m_ClusterSet[i][j].m_index, m_validKey);

			// obviously: group_size >= 3
			int group_size = (int)(1 / minInThree(m_Ratio_train, m_Ratio_test, 1 - m_Ratio_train - m_Ratio_test) + 0.5);
			int train_size = (int)(group_size * m_Ratio_train + 0.5);
			int test_size = (int)(group_size * m_Ratio_test + 0.5);
			int valid_size = group_size - train_size - test_size;

			int groupTMP, trainTMP, testTMP, validTMP;
			while (m_afterSOM->m_ClusterSet[i][j].m_index.size() > 0)
			{
				groupTMP = group_size; trainTMP = train_size; testTMP = test_size; validTMP = valid_size;
				if (m_trainKey.size() < m_TrainNumber_eachNeuron[i][j])
				{
					DP_Resample(m_afterSOM->m_ClusterSet[i][j].m_index, m_trainKey);
					--trainTMP;
					--groupTMP;
				}
				if (m_validKey.size() < m_ValidNumber_EachNeuron[i][j])
				{
					DP_Resample(m_afterSOM->m_ClusterSet[i][j].m_index, m_validKey);
					--validTMP;
					--groupTMP;
				}
				if (m_testKey.size() < m_TestNumber_eachNeuron[i][j])
				{
					DP_Resample(m_afterSOM->m_ClusterSet[i][j].m_index, m_testKey);
					--testTMP;
					--groupTMP;
				}

				while (groupTMP--)
				{
					if (trainTMP != 0 && m_trainKey.size() < m_TrainNumber_eachNeuron[i][j])
					{
						DP_Resample(m_afterSOM->m_ClusterSet[i][j].m_index, m_trainKey);
						--trainTMP;
					}
					if (validTMP != 0 && m_validKey.size() < m_ValidNumber_EachNeuron[i][j])
					{
						DP_Resample(m_afterSOM->m_ClusterSet[i][j].m_index, m_validKey);
						--validTMP;
					}
					if (testTMP != 0 && m_testKey.size() < m_TestNumber_eachNeuron[i][j])
					{
						DP_Resample(m_afterSOM->m_ClusterSet[i][j].m_index, m_testKey);
						--testTMP;
					}
				}
			}
			std::copy(m_trainKey.begin(), m_trainKey.end(), std::back_inserter(m_trainSet));
			std::copy(m_testKey.begin(), m_testKey.end(), std::back_inserter(m_testSet));
			std::copy(m_validKey.begin(), m_validKey.end(), std::back_inserter(m_validSet));
		}
}


void	Method::GetOutputList_Y(vector<int>& Idex)
{
	m_output_Y.clear();
	//Get the last number of each container,which is the output variable 
	for (vector<int>::iterator it = Idex.begin(); it != Idex.end(); ++it)
		m_output_Y.push_back(m_data[*it].back());
}

void	Method::output_Y_bubblesort(vector<double>& data, vector<int>& index)
{
	double tempdata;
	int tempindex;
	for (int i = 0; i < (int)data.size()-1; i++)
		for (int j = i + 1; j < data.size(); j++)
		{
			if (data[i] > data[j])
			{
				tempdata = data[i];
				data[i] = data[j];
				data[j] = tempdata;
				tempindex = index[i];
				index[i] = index[j];
				index[j] = tempindex;
			}
		}
}

std::vector<int>	Method::SSsample(const vector<int>& A,const double& _ratio)
{
	vector<int>	sample;
	double interval = 1.0 / _ratio;
	int m = ceil(interval);
	double Location = rand() % m;
	for (; Location < (double)A.size() - 1; Location += interval)
	{
		// For non-integer locations, round up to next integer value k
		int k = ceil(Location);
		sample.push_back(A[k]);
	}
	return sample;
}

std::vector<int>	Method::remainUnsample(const vector<int>& X, const vector<int>& Y)
{
	vector<int> remainData;
	vector<int>::const_iterator xi = X.begin(), yi = Y.begin();
	for (; xi != X.end(); ++xi)
	{
		if (yi == Y.end())		// Prevent dereference from being out of vector size 
			break;
		if (*xi == *yi)			// if yi equal to xi is found in the Y set
			++yi;
		else
			remainData.push_back(*xi);
	}
	for (; xi != X.end(); ++xi)
		remainData.push_back(*xi);

	return remainData;
}

void	Method::GetDistanceMatrix()
{
	try 
	{
		if (m_data.columns() < 2)
		{
			m_Distance_Matrix.resize(m_data.rows());
			for (int i = 0; i < m_data.rows(); ++i)
				for (int j = i + 1; j < m_data.rows(); ++j)
					m_Distance_Matrix[i].push_back(abs(m_data[i][0] - m_data[j][0]));
		}
		else
		{
			m_Distance_Matrix.resize(m_data.rows());
			for (int i = 0; i < m_data.rows(); ++i)
				for (int j = i + 1; j < m_data.rows(); ++j)
					m_Distance_Matrix[i].push_back(SOMProcess::Getdistance_(m_data[i], m_data[j]));
		}
	}
	catch (std::bad_alloc& ba)
	{
		std::cout << "Bad memory allocation in the process of getting DistanceMatrix!" << std::endl;
		system("pause");
		exit(1);
	}

}

bool	Method::CheckFull(vector<int>& Idexkey,vector<int>& key1, vector<int>& key2, vector<int>& key3,
			const int& key1_lim, const int& key2_lim, const int& key3_lim)
{
	if (!(key1.size() < key1_lim) && !(key2.size() < key2_lim))
	{
		std::copy(Idexkey.begin(), Idexkey.end(), std::back_inserter(key3));
		return true;
	}
	else if (!(key1.size() < key1_lim) && !(key3.size() < key3_lim))
	{
		std::copy(Idexkey.begin(), Idexkey.end(), std::back_inserter(key2));
		return true;
	}	
	else if (!(key2.size() < key2_lim) && !(key3.size() < key3_lim))
	{
		std::copy(Idexkey.begin(), Idexkey.end(), std::back_inserter(key1));
		return true;
	}
	else
		return false;
}

void	Method::SetSamplingParameters()
{
	m_TotalNumber = m_data.size();
	m_SampleNumber = (m_Ratio_train + m_Ratio_test) * m_TotalNumber;
	m_TrainNumber = m_Ratio_train * m_TotalNumber;
	m_TestNumber = m_Ratio_test * m_TotalNumber;
	m_ValidNumber = m_TotalNumber - (m_TrainNumber + m_TestNumber);
}

void	Method::SBSS_P_Process()								// SBSS-P
{
	SetSamplingParameters();
	SetVector();
	GetNumberOfData_eachNeuron();
	GetSampleNumber_eachNeuron_proportion();
	Random_sample();
	AllocateRemainDataInto_validSet();
}

void	Method::SBSS_N_Process()								// SBSS-N
{
	SetSamplingParameters();
	SetVector();
	GetNumberOfData_eachNeuron();
	GetTheta();
	GetSampleNumber_eachNeuron_neyman();
	Random_sample();
	AllocateRemainDataInto_validSet();
}

void	Method::SOMPLEX_Process()								// SOMPLEX
{
	SetSamplingParameters();
	SetVector();
	GetNumberOfData_eachNeuron();
	GetSampleNumber_eachNeuron_proportion();
	GetDistanceMatrix();
	BasedSOM_DPSample();
}

void	Method::SOMPLEX_M_Process()
{
	SetSamplingParameters();
	SetVector();
	GetNumberOfData_eachNeuron();
	GetSampleNumber_eachNeuron_proportion();
	GetDistanceMatrix();
	BasedSOM_MDP_Sample();
}

void	Method::DUPLEX_Process()								// DUPLEX
{
	SetSamplingParameters();
	vector<int>	DUPLEX_tempIndex;
	Sequence(DUPLEX_tempIndex, m_TotalNumber);		
	/* 
	* According to the source code from NNDK,
	* DUPLEX needs initially standardising!
	*/
	DP_Standardise(DUPLEX_tempIndex);
	GetDistanceMatrix();

	// DUPLEX method, step1:
	if (m_TrainNumber > 0)
		DP_InitialSample(DUPLEX_tempIndex, m_trainSet);
	if (m_TestNumber > 0)
		DP_InitialSample(DUPLEX_tempIndex, m_testSet);
	if (m_ValidNumber > 0)
		DP_InitialSample(DUPLEX_tempIndex, m_validSet);

	// DUPLEX method, step2:(resample)
	while (DUPLEX_tempIndex.size() > 0)
	{
		if (m_trainSet.size() < m_TrainNumber)
			DP_Resample(DUPLEX_tempIndex, m_trainSet);
		if (m_testSet.size() < m_TestNumber)
			DP_Resample(DUPLEX_tempIndex, m_testSet);
		if (m_validSet.size() < m_ValidNumber)
			DP_Resample(DUPLEX_tempIndex, m_validSet);

		if (CheckFull(DUPLEX_tempIndex, m_trainSet, m_testSet, m_validSet,
			m_TrainNumber, m_TestNumber, m_ValidNumber))
			break;
	}
}

void	Method::DPModified_Process()
{
	SetSamplingParameters();
	vector<int>	DUPLEX_tempIndex;
	Sequence(DUPLEX_tempIndex, m_TotalNumber);
	DP_Standardise(DUPLEX_tempIndex);
	GetDistanceMatrix();
	
	// obviously: group_size >= 3
	int group_size = (int)(1 / minInThree(m_Ratio_train,m_Ratio_test, 1 - m_Ratio_train - m_Ratio_test) + 0.5);
	int train_size = (int)(group_size * m_Ratio_train + 0.5);
	int test_size = (int)(group_size * m_Ratio_test + 0.5);
	int valid_size = group_size - train_size - test_size;

	int groupTMP, trainTMP, testTMP, validTMP;
	// step1:initially sample
	if (m_TrainNumber > 0)
		DP_InitialSample(DUPLEX_tempIndex, m_trainSet);
	if (m_TestNumber > 0)
		DP_InitialSample(DUPLEX_tempIndex, m_testSet);	
	if (m_ValidNumber > 0)
		DP_InitialSample(DUPLEX_tempIndex, m_validSet);
	// in each group,validation is in the second sample
	while (DUPLEX_tempIndex.size() > 0)
	{
		groupTMP = group_size; trainTMP = train_size; testTMP = test_size; validTMP = valid_size;

		//each data set must be sampled at least once 
		if (m_trainSet.size() < m_TrainNumber)
			DP_Resample(DUPLEX_tempIndex, m_trainSet);
		--trainTMP;
		if(m_testSet.size()<m_TestNumber)
			DP_Resample(DUPLEX_tempIndex, m_testSet);
		--testTMP;
		if (m_validSet.size() < m_ValidNumber)
			DP_Resample(DUPLEX_tempIndex, m_validSet);
		--validTMP;

		groupTMP -= 3;

		// sample in each group
		while (groupTMP--)
		{
			if (trainTMP != 0 && m_trainSet.size() < m_TrainNumber)
			{
				DP_Resample(DUPLEX_tempIndex, m_trainSet);
				--trainTMP;
			}
			if (testTMP != 0 && m_testSet.size() < m_TestNumber)
			{
				DP_Resample(DUPLEX_tempIndex, m_testSet);
				--testTMP;
			}
			if (validTMP != 0 && m_validSet.size() < m_ValidNumber)
			{
				DP_Resample(DUPLEX_tempIndex, m_validSet);
				--validTMP;
			}
		}
	}
}

void	Method::SS_Process()							//SS
{
	SetSamplingParameters();
	vector<int> SS_tempIndex;
	Sequence(SS_tempIndex, m_TotalNumber);

	// Firstly get the output variable list
	GetOutputList_Y(SS_tempIndex);

	// The data are first ordered along the output variable dimension in increasing order
	output_Y_bubblesort(m_output_Y, SS_tempIndex);

	double	calibratingRatio = m_Ratio_train + m_Ratio_test;	

	// First of all, the data is divided into two parts,
	// which need to be sampled and those that do not need to be sampled
	m_sampleKey = SSsample(SS_tempIndex, calibratingRatio);

	// Unsampling data are allocated to validating subset
	m_validSet = remainUnsample(SS_tempIndex, m_sampleKey);

	// Then split sample into systematic testing and training sets
	m_testSet = SSsample(m_sampleKey, m_Ratio_test / calibratingRatio);
	m_trainSet = remainUnsample(m_sampleKey, m_testSet);
}

void	Method::SOM_SS_Process()
{
	SetSamplingParameters();
	SetVector();
	double	calibratingRatio = m_Ratio_train + m_Ratio_test;

	for(int i=0;i<m_afterSOM->m_neuronRow;++i)
		for (int j = 0; j < m_afterSOM->m_neuronColumn; ++j)
		{
			m_trainKey.clear();
			m_testKey.clear();
			m_validKey.clear();

			GetOutputList_Y(m_afterSOM->m_ClusterSet[i][j].m_index);
			output_Y_bubblesort(m_output_Y, m_afterSOM->m_ClusterSet[i][j].m_index);
			m_sampleKey = SSsample(m_afterSOM->m_ClusterSet[i][j].m_index, calibratingRatio);
			m_validKey = remainUnsample(m_afterSOM->m_ClusterSet[i][j].m_index, m_sampleKey);
			m_testKey = SSsample(m_sampleKey, m_Ratio_test / calibratingRatio);
			m_trainKey = remainUnsample(m_sampleKey, m_testKey);

			std::copy(m_trainKey.begin(), m_trainKey.end(), std::back_inserter(m_trainSet));
			std::copy(m_testKey.begin(), m_testKey.end(), std::back_inserter(m_testSet));
			std::copy(m_validKey.begin(), m_testKey.end(), std::back_inserter(m_validSet));
		}
	
}

void	Method::Output_Result(vector<int>& Idex,string& trainFile,string& testFile,string& validFile) const
{
	std::ofstream outfile;
	const char* tab = "\t";
	outfile.open(trainFile);
	vector<int>::const_iterator iter;
	vector<double>::const_iterator iter2;

	outfile << "\"Idex\"";
	for (int i = 1; i < m_Initial_Data.columns(); ++i)
		outfile << "\t\"I\"";
	outfile << "\t\"O\"" << std::endl;
	for (iter = this->m_trainSet.begin(); iter != this->m_trainSet.end(); ++iter)
	{
		outfile << Idex[*iter];
		for (iter2 = m_Initial_Data[*iter].begin(); iter2 != m_Initial_Data[*iter].end(); ++iter2)
			outfile << *tab << *iter2;
		outfile << std::endl;
	}
	outfile.close();

	outfile.open(testFile);
	
	outfile << "\"Idex\"";
	for (int i = 1; i < m_Initial_Data.columns(); ++i)
		outfile << "\t\"I\"";
	outfile << "\t\"O\"" << std::endl;
	for (iter = this->m_testSet.begin(); iter != this->m_testSet.end(); ++iter)
	{
		outfile << Idex[*iter];
		for (iter2 = m_Initial_Data[*iter].begin(); iter2 != m_Initial_Data[*iter].end(); ++iter2)
			outfile << *tab << *iter2;
		outfile << std::endl;
	}
	outfile.close();

	outfile.open(validFile);
	
	outfile << "\"Idex\"";
	for (int i = 1; i < m_Initial_Data.columns(); ++i)
		outfile << "\t\"I\"";
	outfile << "\t\"O\"" << std::endl;
	for (iter = this->m_validSet.begin(); iter != this->m_validSet.end(); ++iter)
	{
		outfile << Idex[*iter];
		for (iter2 = m_Initial_Data[*iter].begin(); iter2 != m_Initial_Data[*iter].end(); ++iter2)
			outfile << *tab << *iter2;
		outfile << std::endl;
	}
	outfile.close();
}