#include<iostream>
#include<stdlib.h>
#include<cstring>
#include<fstream>
#include<iomanip>
#include<memory>
#include"utils.h"
#include"som.h"
#include"method.h"
//#define DEBUG

int main(int argc,const char* argv[])
{
#ifdef DEBUG
	argc = 10;
	argv[1] = "data.txt";argv[2] = "False"; argv[3] = "SOMPLEX"; 
	argv[4] = "1095"; argv[5] = "0.6"; argv[6] = "0.2"; argv[7] = "train.txt"; argv[8] = "test.txt"; argv[9] = "valid.txt";
#endif // DEBUG
	Console console(argc, argv);
	console.check(); //Check the input parameter

	// Get args from console
	string inputFileName = console.getInputFileName();
	int choice = console.getChoice();
	methodType method = console.getMethod();
	int seed = console.getSeed();
	double trainFraction = console.getTrainFrac();
	double testFraction = console.getTestFrac();
	string trainFile = console.getOutTrain();
	string testFile = console.getOutTest();
	string validFile = console.getOutVaild();

	// Get raw data from ...txt
	som::matrix<double> iniData;
	std::vector<int> Idex;
	method::GetTXTData(iniData,inputFileName,Idex);

	som::matrix<double> data;
	if (choice || iniData.columns() == 1)
		data = iniData;
	else{
		data.resize(iniData.rows());
		for (int i = 0; i < iniData.rows(); ++i)
			data[i].push_back(iniData[i].back());
	}

	// Basic preparation 
	std::unique_ptr<method::Method> mainProcess(new method::Method(trainFraction, testFraction, data, iniData));
	srand(seed);
	som::matrix<double>().swap(iniData);

	// Preparation for SOM clustering
	std::unique_ptr<som::SOMProcess> somCluster(new som::SOMProcess());
	if (method <= 2) {
		somCluster->BeforeTrainProcess(data);
		somCluster->Train();
		mainProcess->MoveData(somCluster);
	}
	else 
		som::matrix<double>().swap(data); // Release unused memory
	
	// Select splitting method
	switch (method)
	{
		case SOMPLEX:
			mainProcess->SOMPLEX_Process();
			break;
		case SBSS_P:
			mainProcess->SBSS_P_Process();
			break;
		case SBSS_N:
			mainProcess->SBSS_N_Process();
			break;
		case MDUPLEX:
			mainProcess->DPModified_Process();
			break;
		case DUPLEX:
			mainProcess->DUPLEX_Process();
			break;
		case SS:
			mainProcess->SS_Process();
			break;
		default:
			break;
	}
	mainProcess->Output_Result(Idex, trainFile, testFile, validFile);

	return 0;
}
