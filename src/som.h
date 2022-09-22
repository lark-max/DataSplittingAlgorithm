// SOM algorithm
// Written by Chen J. Y. at 2020, Zhejiang UNV
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
		Inner_Cell():_bias(0.0){}
		~Inner_Cell() = default;

		double	_bias;
		vector<int>	_index;			// Represents the information stored in each neuron 
		vector<double>	W;				// W represents the weight value of each neuron, which is variable	
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
		void	BeforeTrainProcess(matrix<double>& InputData);
		void	Train();
		void	Train(const double& learningRate, const int neighborSize, const int Epochs);
		void	GetActiveNeuron(const vector<double>& Data_I);
		void	UpdateBias();
		void	UpdateWeight(const vector<double>& Data_I);
		void	Standardised(matrix<double>& data);
		static double	GaussRand(double E = 0.0, double V = 1.0);
		static double	GetMean(const vector<double>& _vector);
		static double	GetStdev(const vector<double>& _vector);
		static double	Getdistance_(const vector<double>& m, const vector<double>& n);

	private:
		// basic parameters
		int	totalData;
		int	dimensions;
		int	neuronNumber;
		int	neuronRow;
		int	neuronColumn;
		int	maxIterations;
		int	activeRow;
		int	activeColumn;
		int	neighbourhoodLowerRowIndex_;
		int	neighbourhoodLowerColumnIndex_;
		int	neighbourhoodUpperRowIndex_;
		int	neighbourhoodUpperColumnIndex_;
		int	neighbourhoodSize;
		double	learningRate;

		double	orderingRate;
		int	orderingNeiborSize;
		int	orderingEpochs;

		double	tuningRate;
		int	tuningNeiborSize;
		int	tuningEpochs;

		double	consicience;
		double	bias;

		matrix<double>	InputData;
		matrix<Inner_Cell>	ClusterSet;
	};

}
#endif