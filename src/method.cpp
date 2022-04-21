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

Method::Method()
{}

Method::Method(const double& a, const double& b, matrix<double>& data,matrix<double>& inidata)
	:Ratio_train(a),Ratio_test(b), _data(data), Initial_Data(inidata)
{
}
Method::~Method()
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
	data.resize(row,vector<double>(col));
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

void	Method::Set_vector()
{

	NumberOfData_eachNeuron.resize(afterSOM->neuronRow, vector<int>(afterSOM->neuronColumn));
	TrainNumber_eachNeuron.resize(afterSOM->neuronRow, vector<int>(afterSOM->neuronColumn));
	TestNumber_eachNeuron.resize(afterSOM->neuronRow, vector<int>(afterSOM->neuronColumn));
	ValidNumber_EachNeuron.resize(afterSOM->neuronRow, vector<int>(afterSOM->neuronColumn));
	theta_numerator.resize(afterSOM->neuronRow, vector<double>(afterSOM->neuronColumn));
}

// The number of data in each stratum 
void	Method::GetNumberOfData_eachNeuron()
{
	for (int i = 0; i < afterSOM->neuronRow; ++i)
		for (int j = 0; j < afterSOM->neuronColumn; ++j)
			NumberOfData_eachNeuron[i][j] = afterSOM->ClusterSet[i][j]._index.size();
}

// The amount of data allocated to train, test, and validate at each stratum
// According to equal proportion 
void	Method::GetSampleNumber_eachNeuron_proportion()
{
	for (int i = 0; i < afterSOM->neuronRow; ++i)
		for (int j = 0; j < afterSOM->neuronColumn; ++j)
		{
			TrainNumber_eachNeuron[i][j] = int(Ratio_train / (Ratio_train + Ratio_test)
				* NumberOfData_eachNeuron[i][j] * (Ratio_train + Ratio_test) + 0.5);
			TestNumber_eachNeuron[i][j] = int(Ratio_test / (Ratio_train + Ratio_test)
				* NumberOfData_eachNeuron[i][j] * (Ratio_train + Ratio_test) + 0.5);
			ValidNumber_EachNeuron[i][j] = NumberOfData_eachNeuron[i][j] - TrainNumber_eachNeuron[i][j]
				- TestNumber_eachNeuron[i][j];
		}
}

void	Method::GetSampleNumber_eachNeuron_neyman()
{
	for (int i = 0; i < afterSOM->neuronRow; ++i)
		for (int j = 0; j < afterSOM->neuronColumn; ++j)
		{
			TrainNumber_eachNeuron[i][j] = int(TrainNumber * theta_numerator[i][j] / theta_denominator + 0.5);
			TestNumber_eachNeuron[i][j] = int(TestNumber * theta_numerator[i][j] / theta_denominator + 0.5);

			// Prevent the number of samples from being greater than the number of data in this stratum
			// Sometimes Neyman allocation has this problem!
			if (TrainNumber_eachNeuron[i][j] > NumberOfData_eachNeuron[i][j])
			{
				TrainNumber_eachNeuron[i][j] = NumberOfData_eachNeuron[i][j];
				TestNumber_eachNeuron[i][j] = 0;
			}
			else if (TrainNumber_eachNeuron[i][j] + TestNumber_eachNeuron[i][j] > NumberOfData_eachNeuron[i][j])
				TestNumber_eachNeuron[i][j] = NumberOfData_eachNeuron[i][j] - TrainNumber_eachNeuron[i][j];
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
	for (int i = 0; i < afterSOM->neuronRow; ++i)
		for (int j = 0; j < afterSOM->neuronColumn; ++j)
		{
			int tmp = 0;
			int sampleTimes = 0;
			// trainSet sampling
			if (NumberOfData_eachNeuron[i][j] != 0 && TrainNumber_eachNeuron[i][j] != 0)
			{
				int trainSet_initial_index = 0;
				while (sampleTimes < TrainNumber_eachNeuron[i][j])
				{
					trainSet_initial_index = trainSet.size();
					if (trainSet_initial_index >= TrainNumber)
						return;
					tmp = rand() % afterSOM->ClusterSet[i][j]._index.size();
					trainSet.push_back(afterSOM->ClusterSet[i][j]._index[tmp]);
					it = afterSOM->ClusterSet[i][j]._index.begin() + tmp;
					afterSOM->ClusterSet[i][j]._index.erase(it);
					++sampleTimes;
				}
			}
			// testSet sampling
			if (NumberOfData_eachNeuron[i][j] != 0 && TestNumber_eachNeuron[i][j] != 0)
			{
				int testSet_initial_index = 0;
				sampleTimes = 0;
				while (sampleTimes < TestNumber_eachNeuron[i][j])
				{
					testSet_initial_index = testSet.size();
					if (testSet_initial_index >= TestNumber)
						return;
					tmp = rand() % afterSOM->ClusterSet[i][j]._index.size();
					testSet.push_back(afterSOM->ClusterSet[i][j]._index[tmp]);
					it = afterSOM->ClusterSet[i][j]._index.begin() + tmp;
					afterSOM->ClusterSet[i][j]._index.erase(it);
					++sampleTimes;
				}
			}
		}
}

void	Method::Get_theta()
{
	theta_denominator = 0;
	double temp;
	double average;
	vector<double> var(_data[0].size());


	vector<int>::iterator iter;
	for (int i = 0; i < afterSOM->neuronRow; ++i)
		for (int j = 0; j < afterSOM->neuronColumn; ++j)
		{
			if (NumberOfData_eachNeuron[i][j] == 0)
			{
				theta_numerator[i][j] = 0;
				continue;
			}
			// Calculation formula of standard deviation for each stratum:
			// theta = sqrt(theta_x1^2+theta_x2^2+theta_x3^2+...+theta_xp^2+theta_y^2), p = dimension of data
			// where theta_xi^2 is the intra-stratum variance of component xi,
			// and theta_y^2 is the variance of the output variable, y.
			for (int k = 0; k < _data.columns(); ++k)
			{
				temp = 0;
				average = 0;
				var[k] = 0;

				for (iter = afterSOM->ClusterSet[i][j]._index.begin();
					iter != afterSOM->ClusterSet[i][j]._index.end(); ++iter)
					temp += _data[*iter][k];

				average = temp / NumberOfData_eachNeuron[i][j];

				for (iter = afterSOM->ClusterSet[i][j]._index.begin();
					iter != afterSOM->ClusterSet[i][j]._index.end(); ++iter)
					var[k] += pow(_data[*iter][k] - average, 2)
					/ NumberOfData_eachNeuron[i][j];
			}
			double _temp = accumulate(var.begin(), var.end(), 0.0);
			theta_numerator[i][j] = NumberOfData_eachNeuron[i][j] * sqrt(_temp);
			theta_denominator += theta_numerator[i][j];
		}

}

// After sampling, the remaining data of each stratum is allocated to the valid set 
void	Method::AllocateRemainDataInto_validSet()
{
	vector<int>::iterator it;
	for (int i = 0; i < afterSOM->neuronRow; ++i)
		for (int j = 0; j < afterSOM->neuronColumn; ++j)
			for (it = afterSOM->ClusterSet[i][j]._index.begin(); it != afterSOM->ClusterSet[i][j]._index.end(); ++it)
				validSet.push_back(*it);

}

void	Method::GetThisColumnData(vector<double>& _column_,int& j, const vector<int>& CurrentDataIndex)
{
	for (int i = 0; i < _column_.size(); ++i)
		_column_[i] = _data[CurrentDataIndex[i]][j];
}

// Calculate the mean and stdev  of each column
// Then datum = (datum - mean)/stdev;
void	Method::DP_Standardise(const vector<int>& CurrentDataIndex)
{
	if (_data.columns()<2 || CurrentDataIndex.size() < 2)
		return;
	vector<double> _column_(CurrentDataIndex.size());
	for (int j = 0; j < _data.columns(); ++j)
	{
		GetThisColumnData(_column_, j, CurrentDataIndex);
		double mean_ = SOMProcess::GetMean(_column_);
		double stdev_ =SOMProcess::GetStdev(_column_);
		for (int i = 0; i < CurrentDataIndex.size(); ++i)
			_data[CurrentDataIndex[i]][j] = (_data[CurrentDataIndex[i]][j] - mean_) / stdev_;
	}
}

// The first step of DUPLEX method: find two data points with the largest Euclidean distance 
void	Method::DP_initialSample(vector<int>& ini_Data, vector<int>& Key)
{
	if (ini_Data.size() < 2)
		return;
	vector<int>::iterator a = ini_Data.begin();
	vector<int>::iterator b = ini_Data.begin() + 1;
	double	maxDistance = 0;
	for (vector<int>::iterator t = ini_Data.begin(); t != ini_Data.end() - 1; ++t)
		for (vector<int>::iterator tt = t + 1; tt != ini_Data.end(); ++tt)
		{
			//double _distance = Distance_Matrix[*t][*tt];
			//double _distance = abs(_data[*t][0] - _data[*tt][0]);
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
void	Method::DP_resample(vector<int>& ini_Data, vector<int>& Key)
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
			//double temp_dist = Distance_Matrix[*t][*s];
			//double temp_dist = abs(_data[*t][0] - _data[*s][0]);
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
		{
		}
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
	for (int i = 0; i < afterSOM->neuronRow; ++i)
		for (int j = 0; j < afterSOM->neuronColumn; ++j)
		{
			trainKey.clear();
			testKey.clear();
			validKey.clear();

			if (TrainNumber_eachNeuron[i][j] > 0)
				DP_initialSample(afterSOM->ClusterSet[i][j]._index, trainKey);
			if (TestNumber_eachNeuron[i][j] > 0)
				DP_initialSample(afterSOM->ClusterSet[i][j]._index, testKey);
			if (ValidNumber_EachNeuron[i][j] > 0)
				DP_initialSample(afterSOM->ClusterSet[i][j]._index, validKey);
			while (afterSOM->ClusterSet[i][j]._index.size() > 0)
			{
				if (trainKey.size() < TrainNumber_eachNeuron[i][j])
					DP_resample(afterSOM->ClusterSet[i][j]._index, trainKey);
				if (testKey.size() < TestNumber_eachNeuron[i][j])
					DP_resample(afterSOM->ClusterSet[i][j]._index, testKey);
				if (validKey.size() < ValidNumber_EachNeuron[i][j])
					DP_resample(afterSOM->ClusterSet[i][j]._index, validKey);

				if (CheckFull(afterSOM->ClusterSet[i][j]._index, trainKey, testKey, validKey,
					TrainNumber_eachNeuron[i][j], TestNumber_eachNeuron[i][j], ValidNumber_EachNeuron[i][j]))
					break;
			}
			std::copy(trainKey.begin(), trainKey.end(), std::back_inserter(trainSet));
			std::copy(testKey.begin(), testKey.end(), std::back_inserter(testSet));
			std::copy(validKey.begin(), validKey.end(), std::back_inserter(validSet));
		}
}


void	Method::BasedSOM_DP_M_Sample()
{
	for (int i = 0; i < afterSOM->neuronRow; ++i)
		for (int j = 0; j < afterSOM->neuronColumn; ++j)
		{
			trainKey.clear();
			testKey.clear();
			validKey.clear();

			if (TrainNumber_eachNeuron[i][j] > 0)
				DP_initialSample(afterSOM->ClusterSet[i][j]._index, trainKey);
			if (TestNumber_eachNeuron[i][j] > 0)
				DP_initialSample(afterSOM->ClusterSet[i][j]._index, testKey);
			if (ValidNumber_EachNeuron[i][j] > 0)
				DP_initialSample(afterSOM->ClusterSet[i][j]._index, validKey);

			// obviously: group_size >= 3
			int group_size = (int)(1 / minInThree(Ratio_train, Ratio_test, 1 - Ratio_train - Ratio_test) + 0.5);
			int train_size = (int)(group_size * Ratio_train + 0.5);
			int test_size = (int)(group_size * Ratio_test + 0.5);
			int valid_size = group_size - train_size - test_size;

			int groupTMP, trainTMP, testTMP, validTMP;
			while (afterSOM->ClusterSet[i][j]._index.size() > 0)
			{
				groupTMP = group_size; trainTMP = train_size; testTMP = test_size; validTMP = valid_size;
				if (trainKey.size() < TrainNumber_eachNeuron[i][j])
				{
					DP_resample(afterSOM->ClusterSet[i][j]._index, trainKey);
					--trainTMP;
					--groupTMP;
				}
				if (validKey.size() < ValidNumber_EachNeuron[i][j])
				{
					DP_resample(afterSOM->ClusterSet[i][j]._index, validKey);
					--validTMP;
					--groupTMP;
				}
				if (testKey.size() < TestNumber_eachNeuron[i][j])
				{
					DP_resample(afterSOM->ClusterSet[i][j]._index, testKey);
					--testTMP;
					--groupTMP;
				}

				while (groupTMP--)
				{
					if (trainTMP != 0 && trainKey.size() < TrainNumber_eachNeuron[i][j])
					{
						DP_resample(afterSOM->ClusterSet[i][j]._index, trainKey);
						--trainTMP;
					}
					if (validTMP != 0 && validKey.size() < ValidNumber_EachNeuron[i][j])
					{
						DP_resample(afterSOM->ClusterSet[i][j]._index, validKey);
						--validTMP;
					}
					if (testTMP != 0 && testKey.size() < TestNumber_eachNeuron[i][j])
					{
						DP_resample(afterSOM->ClusterSet[i][j]._index, testKey);
						--testTMP;
					}
				}
			}
			std::copy(trainKey.begin(), trainKey.end(), std::back_inserter(trainSet));
			std::copy(testKey.begin(), testKey.end(), std::back_inserter(testSet));
			std::copy(validKey.begin(), validKey.end(), std::back_inserter(validSet));
		}
}


void	Method::GetOutputList_Y(vector<int>& Idex)
{
	output_Y.clear();
	//Get the last number of each container,which is the output variable 
	for (vector<int>::iterator it = Idex.begin(); it != Idex.end(); ++it)
		output_Y.push_back(_data[*it].back());
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
		if (_data.columns() < 2)
		{
			Distance_Matrix.resize(_data.rows());
			for (int i = 0; i < _data.rows(); ++i)
				for (int j = i + 1; j < _data.rows(); ++j)
					Distance_Matrix[i].push_back(abs(_data[i][0] - _data[j][0]));
		}
		else
		{
			Distance_Matrix.resize(_data.rows());
			for (int i = 0; i < _data.rows(); ++i)
				for (int j = i + 1; j < _data.rows(); ++j)
					Distance_Matrix[i].push_back(SOMProcess::Getdistance_(_data[i], _data[j]));
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
	TotalNumber = _data.size();
	SampleNumber = (Ratio_train + Ratio_test) * TotalNumber;
	TrainNumber = Ratio_train * TotalNumber;
	TestNumber = Ratio_test * TotalNumber;
	ValidNumber = TotalNumber - (TrainNumber + TestNumber);
}

void	Method::SBSS_P_Process()								// SBSS-P
{
	SetSamplingParameters();
	Set_vector();
	GetNumberOfData_eachNeuron();
	GetSampleNumber_eachNeuron_proportion();
	Random_sample();
	AllocateRemainDataInto_validSet();
}

void	Method::SBSS_N_Process()								// SBSS-N
{
	SetSamplingParameters();
	Set_vector();
	GetNumberOfData_eachNeuron();
	Get_theta();
	GetSampleNumber_eachNeuron_neyman();
	Random_sample();
	AllocateRemainDataInto_validSet();
}

void	Method::SOMPLEX_Process()								// SOMPLEX
{
	SetSamplingParameters();
	Set_vector();
	GetNumberOfData_eachNeuron();
	GetSampleNumber_eachNeuron_proportion();
	GetDistanceMatrix();
	BasedSOM_DPSample();
}

void	Method::SOMPLEX_M_Process()
{
	SetSamplingParameters();
	Set_vector();
	GetNumberOfData_eachNeuron();
	GetSampleNumber_eachNeuron_proportion();
	GetDistanceMatrix();
	BasedSOM_DP_M_Sample();
}

void	Method::DUPLEX_Process()								// DUPLEX
{
	SetSamplingParameters();
	vector<int>	DUPLEX_tempIndex;
	Sequence(DUPLEX_tempIndex, TotalNumber);		
	/* 
	* According to the source code from NNDK,
	* DUPLEX needs initially standardising!
	*/
	DP_Standardise(DUPLEX_tempIndex);
	GetDistanceMatrix();

	// DUPLEX method, step1:
	if (TrainNumber > 0)
		DP_initialSample(DUPLEX_tempIndex, trainSet);
	if (TestNumber > 0)
		DP_initialSample(DUPLEX_tempIndex, testSet);
	if (ValidNumber > 0)
		DP_initialSample(DUPLEX_tempIndex, validSet);

	// DUPLEX method, step2:(resample)
	while (DUPLEX_tempIndex.size() > 0)
	{
		if (trainSet.size() < TrainNumber)
			DP_resample(DUPLEX_tempIndex, trainSet);
		if (testSet.size() < TestNumber)
			DP_resample(DUPLEX_tempIndex, testSet);
		if (validSet.size() < ValidNumber)
			DP_resample(DUPLEX_tempIndex, validSet);

		if (CheckFull(DUPLEX_tempIndex, trainSet, testSet, validSet,
			TrainNumber, TestNumber, ValidNumber))
			break;
	}
}

void	Method::DPModified_Process()
{
	SetSamplingParameters();
	vector<int>	DUPLEX_tempIndex;
	Sequence(DUPLEX_tempIndex, TotalNumber);
	DP_Standardise(DUPLEX_tempIndex);
	GetDistanceMatrix();
	
	// obviously: group_size >= 3
	int group_size = (int)(1 / minInThree(Ratio_train,Ratio_test, 1 - Ratio_train - Ratio_test) + 0.5);
	int train_size = (int)(group_size * Ratio_train + 0.5);
	int test_size = (int)(group_size * Ratio_test + 0.5);
	int valid_size = group_size - train_size - test_size;

	int groupTMP, trainTMP, testTMP, validTMP;
	// step1:initially sample
	if (TrainNumber > 0)
		DP_initialSample(DUPLEX_tempIndex, trainSet);
	if (TestNumber > 0)
		DP_initialSample(DUPLEX_tempIndex, testSet);	
	if (ValidNumber > 0)
		DP_initialSample(DUPLEX_tempIndex, validSet);
	// in each group,validation is in the second sample
	while (DUPLEX_tempIndex.size() > 0)
	{
		groupTMP = group_size; trainTMP = train_size; testTMP = test_size; validTMP = valid_size;

		//each data set must be sampled at least once 
		if (trainSet.size() < TrainNumber)
			DP_resample(DUPLEX_tempIndex, trainSet);
		--trainTMP;
		if(testSet.size()<TestNumber)
			DP_resample(DUPLEX_tempIndex, testSet);
		--testTMP;
		if (validSet.size() < ValidNumber)
			DP_resample(DUPLEX_tempIndex, validSet);
		--validTMP;

		groupTMP -= 3;

		// sample in each group
		while (groupTMP--)
		{
			if (trainTMP != 0 && trainSet.size() < TrainNumber)
			{
				DP_resample(DUPLEX_tempIndex, trainSet);
				--trainTMP;
			}
			if (testTMP != 0 && testSet.size() < TestNumber)
			{
				DP_resample(DUPLEX_tempIndex, testSet);
				--testTMP;
			}
			if (validTMP != 0 && validSet.size() < ValidNumber)
			{
				DP_resample(DUPLEX_tempIndex, validSet);
				--validTMP;
			}
		}
	}
}

void	Method::SS_Process()							//SS Sampler
{
	SetSamplingParameters();
	vector<int> SS_tempIndex;
	Sequence(SS_tempIndex, TotalNumber);

	// Firstly get the output variable list
	GetOutputList_Y(SS_tempIndex);

	// The data are first ordered along the output variable dimension in increasing order
	output_Y_bubblesort(output_Y, SS_tempIndex);

	double	calibratingRatio = Ratio_train + Ratio_test;	

	// First of all, the data is divided into two parts,
	// which need to be sampled and those that do not need to be sampled
	sampleKey = SSsample(SS_tempIndex, calibratingRatio);

	// Unsampling data are allocated to validating subset
	validSet = remainUnsample(SS_tempIndex, sampleKey);

	// Then split sample into systematic testing and training sets
	testSet = SSsample(sampleKey, Ratio_test / calibratingRatio);
	trainSet = remainUnsample(sampleKey, testSet);
}

void	Method::SOM_SS_Process()
{
	SetSamplingParameters();
	Set_vector();
	double	calibratingRatio = Ratio_train + Ratio_test;

	for(int i=0;i<afterSOM->neuronRow;++i)
		for (int j = 0; j < afterSOM->neuronColumn; ++j)
		{
			trainKey.clear();
			testKey.clear();
			validKey.clear();

			GetOutputList_Y(afterSOM->ClusterSet[i][j]._index);
			output_Y_bubblesort(output_Y, afterSOM->ClusterSet[i][j]._index);
			sampleKey = SSsample(afterSOM->ClusterSet[i][j]._index, calibratingRatio);
			validKey = remainUnsample(afterSOM->ClusterSet[i][j]._index, sampleKey);
			testKey = SSsample(sampleKey, Ratio_test / calibratingRatio);
			trainKey = remainUnsample(sampleKey, testKey);

			std::copy(trainKey.begin(), trainKey.end(), std::back_inserter(trainSet));
			std::copy(testKey.begin(), testKey.end(), std::back_inserter(testSet));
			std::copy(validKey.begin(), testKey.end(), std::back_inserter(validSet));
		}
	
}

void	Method::Output_result(vector<int>& Idex,string& trainFile,string& testFile,string& validFile) const
{
	std::ofstream outfile;
	const char* tab = "\t";
	outfile.open(trainFile);
	vector<int>::const_iterator iter;
	vector<double>::const_iterator iter2;

	outfile << "\"Idex\"";
	for (int i = 1; i < Initial_Data.columns(); ++i)
		outfile << "\t\"I\"";
	outfile << "\t\"O\"" << std::endl;
	for (iter = this->trainSet.begin(); iter != this->trainSet.end(); ++iter)
	{
		outfile << Idex[*iter];
		for (iter2 = Initial_Data[*iter].begin(); iter2 != Initial_Data[*iter].end(); ++iter2)
			outfile << *tab << *iter2;
		outfile << std::endl;
	}
	outfile.close();

	outfile.open(testFile);
	
	outfile << "\"Idex\"";
	for (int i = 1; i < Initial_Data.columns(); ++i)
		outfile << "\t\"I\"";
	outfile << "\t\"O\"" << std::endl;
	for (iter = this->testSet.begin(); iter != this->testSet.end(); ++iter)
	{
		outfile << Idex[*iter];
		for (iter2 = Initial_Data[*iter].begin(); iter2 != Initial_Data[*iter].end(); ++iter2)
			outfile << *tab << *iter2;
		outfile << std::endl;
	}
	outfile.close();

	outfile.open(validFile);
	
	outfile << "\"Idex\"";
	for (int i = 1; i < Initial_Data.columns(); ++i)
		outfile << "\t\"I\"";
	outfile << "\t\"O\"" << std::endl;
	for (iter = this->validSet.begin(); iter != this->validSet.end(); ++iter)
	{
		outfile << Idex[*iter];
		for (iter2 = Initial_Data[*iter].begin(); iter2 != Initial_Data[*iter].end(); ++iter2)
			outfile << *tab << *iter2;
		outfile << std::endl;
	}
	outfile.close();
}
