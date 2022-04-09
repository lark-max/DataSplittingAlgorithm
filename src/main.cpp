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
	argv[1] = "data.txt";argv[2] = "false"; argv[3] = "SOMPLEX"; 
	argv[4] = "1000"; argv[5] = "0.6"; argv[6] = "0.2"; argv[7] = "train.txt"; argv[8] = "test.txt"; argv[9] = "valid.txt";
#endif // DEBUG
	Console console(argc, argv);
	console.check(); //Check the input parameter

	string inputFileName = console.getInputFileName();
	int choice = console.getChoice();
	methodType method = console.getMethod();
	int seed = console.getSeed();
	float trainFraction = console.getTrainFrac();
	float testFraction = console.getTestFrac();
	string trainFile = console.getOutTrain();
	string testFile = console.getOutTest();
	string validFile = console.getOutVaild();

	som::matrix<double> inidata;
	std::vector<int> Idex;
	method::GetTXTData(inidata,inputFileName,Idex);

	som::matrix<double> data;
	if (choice)
		data = inidata;
	else{
		data.resize(inidata.rows());
		for (int i = 0; i < inidata.rows(); ++i)
			data[i].push_back(inidata[i].back());
	}
	
	// Basic preparation 
	std::unique_ptr<method::Method> mainProcess(new method::Method(trainFraction, testFraction, data, inidata));
	srand(seed);
	std::unique_ptr<som::SOMProcess> somCluster(new som::SOMProcess());
	if (method <= 2) {
		somCluster->BeforeTrainProcess(data);
		somCluster->Train();
		mainProcess->StructCopy(somCluster);
	}
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
	mainProcess->Output_result(Idex, trainFile, testFile, validFile);

	return 0;
}
