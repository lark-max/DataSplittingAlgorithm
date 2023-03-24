/*
* This file contains the implementation classes for the SOM neural network
* Written by Chen J. Y. at 2020, Zhejiang UNV
*/
#ifndef _SOMMethod
#define	_SOMMethod
#include<vector>
#include<string>

namespace method {
	class Method;
	void	Sequence(std::vector<int>& Idex, int max, int min = 0);
}

namespace som {
	using std::vector;

	// class matrix,basic template class
	template <typename T>
	class matrix : public vector< vector<T> >
	{
	public:
		matrix() = default;
		matrix(int row, int col);
		~matrix() = default;
	
		int	rows()	const;
		int	columns()	const;
		T	Max()	const;
		T	Min()	const;
		T	Max(int j)	const;
		T	Min(int j)	const;
	};

	template<typename T>
	matrix<T>::matrix(int row, int col)
	{
		this->resize(row);
		for (int i = 0; i < row; ++i)
			(*this)[i].resize(col);		
	}

	template <typename T>
	inline int matrix<T>::rows() const
	{
		return (*this).size();
	}

	template <typename T>
	inline int matrix<T>::columns() const
	{
		if (this->rows() == 0)
			return 0;
		else
			return (*this)[0].size();
	}

	template <typename T>
	T	 matrix<T>::Max()	const
	{
		T tmpMAX = (*this)[0][0];
		for (int i = 0; i < (*this).rows(); ++i)
			for (int j = 0; j < (*this).columns(); ++j)
				if ((*this)[i][j] > tmpMAX)
					tmpMAX = (*this)[i][j];
		return tmpMAX;
	}

	template <typename T>
	T	matrix<T>::Min() const
	{
		T tmpMIN = (*this)[0][0];
		for (int i = 0; i < (*this).rows(); ++i)
			for (int j = 0; j < (*this).columns(); ++j)
				if ((*this)[i][j] < tmpMIN)
					tmpMIN = (*this)[i][j];
		return tmpMIN;
	}

	template <typename T>
	T	 matrix<T>::Max(int j)	const
	{
		T tmpMAX = (*this)[0][j];
		for (int i = 0; i < (*this).rows(); ++i)
				if ((*this)[i][j] > tmpMAX)
					tmpMAX = (*this)[i][j];
		return tmpMAX;
	}

	template <typename T>
	T	 matrix<T>::Min(int j)	const
	{
		T tmpMIN = (*this)[0][j];
		for (int i = 0; i < (*this).rows(); ++i)
				if ((*this)[i][j] < tmpMIN)
					tmpMIN = (*this)[i][j];
		return tmpMIN;
	}
	// class Inner_Cell
	// this struct is using for setting each neuron(belongs to the vector "SOMProcess") as a struct,
	// in order to use the member functions and objects
	class Inner_Cell
	{
	public:
		Inner_Cell():m_bias(0.0){}
		~Inner_Cell() = default;

		double	m_bias;
		vector<int>	m_index;			// Represents the information stored in each neuron 
		vector<double>	m_weight;				// m_weight represents the weight value of each neuron, which is variable	
	};

	// main class
	class SOMProcess : public matrix<Inner_Cell>
	{
		friend class method::Method;
		
	public:
		// Basic functions
		void	RandomlyInitialize();
		void	OutputClusterSet();
		void	SetNeuronSize();
		void	Normalize(matrix<double>& data);
		void	GetNeighbor();
		void	AllocateProcess(const vector<double>& Data_I, int& i);
		void	BeforeTrainProcess(matrix<double>& m_InputData);
		void	Train();
		void	Train(const double& m_learningRate, const int neighborSize, const int Epochs);
		void	GetActiveNeuron(const vector<double>& Data_I);
		void	UpdateBias();
		void	UpdateWeight(const vector<double>& Data_I);
		void	Standardised(matrix<double>& data);
		static double	GaussRand(double E = 0.0, double V = 1.0);
		static double	GetMean(const vector<double>& _vector);
		static double	GetStdev(const vector<double>& _vector);
		static double	Getdistance(const vector<double>& m, const vector<double>& n);

	private:
		// basic parameters
		int	m_totalData;
		int	m_dimensions;
		int	m_neuronNumber;
		int	m_neuronRow;
		int	m_neuronColumn;
		int	m_maxIterations;
		int	m_activeRow;
		int	m_activeColumn;
		int	m_neighbourhoodLowerRowIndex;
		int	m_neighbourhoodLowerColumnIndex;
		int	m_neighbourhoodUpperRowIndex;
		int	m_neighbourhoodUpperColumnIndex;
		int	m_neighbourhoodSize;
		double	m_learningRate;

		double	m_orderingRate;
		int	m_orderingNeiborSize;
		int	m_orderingEpochs;

		double	m_tuningRate;
		int	m_tuningNeiborSize;
		int	m_tuningEpochs;

		double	m_consicience;
		double	m_bias;

		matrix<double>	m_InputData;
		matrix<Inner_Cell>	m_ClusterSet;
	};

}
#endif