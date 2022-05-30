// This file includes several methods for data splitting
// Written by Chen J. Y. at 2020, Zhejiang UNV
#ifndef METHOD_
#define METHOD_
#include"som.h"
#define minInThree(x,y,z)   ((x) < (y) ? ((x) < (z) ?  (x): (z)) : ((y) < (z) ? (y) : (z)))
#define maxInThree(x,y,z)	((x) > (y) ? ((x) > (z) ?  (x): (z)) : ((y) > (z) ? (y) : (z)))
#define Distance_matrix(i,j)  (((i) < (j)) ? (Distance_Matrix[(i)][(j)-(i)-1]) : (Distance_Matrix[(j)][(i)-(j)-1]))

namespace method{
	using std::vector;
	using std::string;

	void	GetTXTData(som::matrix<double>& data,string& str, vector<int>& Idex);

	class Method
	{
	public:

		Method() = default;
		Method(const double& a, const double& b, som::matrix<double>& data, 
					som::matrix<double>& iniData);
		~Method() = default;
		// basic methods	
		void	MoveData(std::unique_ptr<som::SOMProcess>& ptr){afterSOM = std::move(ptr);}
		void	Set_vector();
		void	GetNumberOfData_eachNeuron();
		void	GetSampleNumber_eachNeuron_proportion();
		void	GetSampleNumber_eachNeuron_neyman();
		void	SetSamplingParameters();
		void	GetDistanceMatrix();
		

		// SBSS-P
		void	SBSS_P_Process();
		void	Random_sample();
		void	AllocateRemainDataInto_validSet();

		// SBSS-N
		void	SBSS_N_Process();
		void	Get_theta();

		// SOM-DUPLEX(SOMPLEX)
		void	SOMPLEX_Process();
		void	BasedSOM_DPSample();
		
		// SOMPLEX-M
		void	SOMPLEX_M_Process();
		void	BasedSOM_DP_M_Sample();

		// DUPLEX
		void	DUPLEX_Process();
		void	DP_Standardise(const vector<int>& CurrentDataIndex);
		void	DP_initialSample(vector<int>& ini_Data, vector<int>& Key);
		void	DP_resample(vector<int>& ini_Data, vector<int>& Key);
		void	GetThisColumnData(vector<double>& _column_, int& j, const vector<int>& CurrentDataIndex);
		bool	CheckFull(vector<int>& Idexkey,vector<int>& key1,vector<int>& key2,vector<int>& key3,
							const int& key1_lim, const int& key2_lim, const int& key3_lim);

		// DUPLEX modified(MDUPLEX)
		void	DPModified_Process();

		// SS Sampler
		void	SS_Process();
		void	GetOutputList_Y(vector<int>& Idex);
		void	output_Y_bubblesort(vector<double>& data, vector<int>& index);
		vector<int>	SSsample(const vector<int>& A, const double& _ratio);
		vector<int>	remainUnsample(const vector<int>& X, const vector<int>& Y);

		// SOM-SS
		void	SOM_SS_Process();


		// Output as txt
		void	Output_result(vector<int>& Idex,string& trainFile, string& testFile, string& validFile)	const;

	private:
		int	TotalNumber;
		int	SampleNumber;
		int	TrainNumber;
		int	TestNumber;
		int	ValidNumber;
		double	Ratio_train;
		double	Ratio_test;
		double	theta_denominator;
		vector<int>	trainSet;
		vector<int>	testSet;
		vector<int>	validSet;
		vector<vector<double>>	theta_numerator;
		vector<int>	trainKey;
		vector<int>	testKey;
		vector<int>	validKey;
		vector<int>	sampleKey;
		vector<double>	output_Y;
		som::matrix<double>	_data;		
		const som::matrix<double>	Initial_Data;
		vector<vector<int>>	NumberOfData_eachNeuron;
		vector<vector<int>>	TrainNumber_eachNeuron;
		vector<vector<int>>	TestNumber_eachNeuron;
		vector<vector<int>>	ValidNumber_EachNeuron;
		std::unique_ptr<som::SOMProcess>	afterSOM;
		som::matrix<double>	Distance_Matrix;
	};
	
}

#endif // !METHOD_
