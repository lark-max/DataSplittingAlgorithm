// This file includes several methods for data splitting
// Written by Chen J. Y. at 2020, Zhejiang UNV
#ifndef METHOD_
#define METHOD_
#include"som.h"
#define minInThree(x,y,z)   ((x) < (y) ? ((x) < (z) ?  (x): (z)) : ((y) < (z) ? (y) : (z)))
#define maxInThree(x,y,z)	((x) > (y) ? ((x) > (z) ?  (x): (z)) : ((y) > (z) ? (y) : (z)))
#define Distance_matrix(i,j)  (((i) < (j)) ? (m_Distance_Matrix[(i)][(j)-(i)-1]) : (m_Distance_Matrix[(j)][(i)-(j)-1]))

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
		void	MoveData(std::unique_ptr<som::SOMProcess>& ptr){m_afterSOM = std::move(ptr);}
		void	SetVector();
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
		void	GetTheta();

		// SOM-DUPLEX(SOMPLEX)
		void	SOMPLEX_Process();
		void	BasedSOM_DPSample();
		
		// SOMPLEX-M
		void	SOMPLEX_M_Process();
		void	BasedSOM_MDP_Sample();

		// DUPLEX
		void	DUPLEX_Process();
		void	DP_Standardise(const vector<int>& CurrentDataIndex);
		void	DP_InitialSample(vector<int>& ini_Data, vector<int>& Key);
		void	DP_Resample(vector<int>& ini_Data, vector<int>& Key);
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
		void	Output_Result(vector<int>& Idex,string& trainFile, string& testFile, string& validFile)	const;

	private:
		int	m_TotalNumber;
		int	m_SampleNumber;
		int	m_TrainNumber;
		int	m_TestNumber;
		int	m_ValidNumber;
		double	m_Ratio_train;
		double	m_Ratio_test;
		double	m_theta_denominator;
		vector<int>	m_trainSet;
		vector<int>	m_testSet;
		vector<int>	m_validSet;
		vector<vector<double>>	m_theta_numerator;
		vector<int>	m_trainKey;
		vector<int>	m_testKey;
		vector<int>	m_validKey;
		vector<int>	m_sampleKey;
		vector<double>	m_output_Y;
		som::matrix<double>	m_data;		
		const som::matrix<double>	m_Initial_Data;
		vector<vector<int>>	m_NumberOfData_eachNeuron;
		vector<vector<int>>	m_TrainNumber_eachNeuron;
		vector<vector<int>>	m_TestNumber_eachNeuron;
		vector<vector<int>>	m_ValidNumber_EachNeuron;
		std::unique_ptr<som::SOMProcess>	m_afterSOM;
		som::matrix<double>	m_Distance_Matrix;
	};
	
}

#endif // !METHOD_
