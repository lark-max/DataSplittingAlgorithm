#include<iostream>
#include<stdlib.h>
#include<cstring>
#include<fstream>
#include<iomanip>
#include<memory>
#include"som.h"
#include"method.h"
#define DEBUG

int main(int argc,const char* argv[])
{
#ifdef DEBUG
	argc = 10;
	argv[1] = "data.txt";argv[2] = "only_output"; argv[3] = "SS"; argv[4] = "1001";
	argv[5] = "0.6"; argv[6] = "0.2"; argv[7] = "train.txt"; argv[8] = "test.txt"; argv[9] = "valid.txt";
#endif // DEBUG

	if (argc < 10)
	{
		std::cerr<< "The program lacks some necessary parameters!" << std::endl;
		system("pause");
		exit(1);
	}
	std::string file_name_ = argv[1], trainFile = argv[7], testFile = argv[8], validFile = argv[9];
	int seed = atoi(argv[4]);
	double trainFraction, testFraction;

	som::matrix<double> inidata;
	std::vector<int> Idex;
	try
	{
		method::GetTXTData(inidata,file_name_,Idex);
	}
	catch (std::bad_alloc& ba)
	{
		std::cerr << "Bad memory allocation in the process of input data."
			<< "Probably the data are too large" << std::endl;
		system("pause");
		exit(1);
	}

	som::matrix<double> data;
	if (!strcmp(argv[2],"include_input"))
		data = inidata;
	else if (!strcmp(argv[2],"only_output"))
	{
		data.resize(inidata.rows());
		for (int i = 0; i < inidata.rows(); ++i)
			data[i].push_back(inidata[i].back());
	}
	else
	{
		std::cerr << "Error: Incorrect schema name!" << std::endl;
		system("pause");
		exit(1);
	}

	trainFraction = atof(argv[5]);
	testFraction = atof(argv[6]);
	if (trainFraction + testFraction >= 1.0)
	{
		std::cerr << "Error: The splitting ratio is greater than 1!" << std::endl;
		system("pause");
		exit(1);
	}
	try
	{
		if (trainFraction + testFraction >= 1.0)
			throw trainFraction + testFraction;
	}
	catch (double)
	{
		std::cerr << "Error: The splitting ratio is greater than 1!" << std::endl;
		system("pause");
		exit(1);
	}
	// Basic preparation 
	std::unique_ptr<method::Method> test(new method::Method(trainFraction, testFraction, data, inidata));
	srand(seed);

	if (!strcmp(argv[3],"SS"))		// SS
	{
		test->SS_Process();
		test->Output_result(Idex, trainFile, testFile, validFile);
	}
	else if (!strcmp(argv[3], "DUPLEX"))	// DUPLEX
	{
		test->DUPLEX_Process();
		test->Output_result(Idex, trainFile, testFile, validFile);
	}
	else if (!strcmp(argv[3], "DUPLEX-M")) // modified DUPLEX
	{
		test->DPModified_Process();
		test->Output_result(Idex, trainFile, testFile, validFile);
	}
	// Based on SOM
	else if (!strcmp(argv[3], "SOMPLEX") || !strcmp(argv[3], "SBSS-P") || !strcmp(argv[3], "SBSS-N") || !strcmp(argv[3], "SOM-SS") || !strcmp(argv[3], "SOMPLEX-M"))	
	{
		std::unique_ptr<som::SOMProcess> test2(new som::SOMProcess());
		if (!strcmp(argv[3], "SOMPLEX"))
		{
			// SOMPLEX
			test2->BeforeTrainProcess(data);
			test2->Train();
			test2->OutputClusterSet();
			test->StructCopy(test2);
			test->SOMPLEX_Process();
			test->Output_result(Idex, trainFile, testFile, validFile);
		}
		else if (!strcmp(argv[3], "SBSS-P"))
		{
			// SBSS-P
			test2->BeforeTrainProcess(data);
			test2->Train();
			test->StructCopy(test2);
			test->SBSS_P_Process();
			test->Output_result(Idex, trainFile, testFile, validFile);
		}
		else if(!strcmp(argv[3], "SBSS-N"))
		{
			// SBSS-N
			test2->BeforeTrainProcess(data);
			test2->Train();
			test->StructCopy(test2);
			test->SBSS_N_Process();
			test->Output_result(Idex, trainFile, testFile, validFile);
		}
		else if (!strcmp(argv[3], "SOM-SS"))
		{
			//SOM-SS
			test2->BeforeTrainProcess(data);
			test2->Train();
			test->StructCopy(test2);
			test->SOM_SS_Process();
			test->Output_result(Idex, trainFile, testFile, validFile);
		}
		else if (!strcmp(argv[3], "SOMPLEX-M"))
		{
			//SOMPLEX-M
			test2->BeforeTrainProcess(data);
			test2->Train();
			test->StructCopy(test2);
			test->SOMPLEX_M_Process();
			test->Output_result(Idex, trainFile, testFile, validFile);
		}
		else
		{}
	}
	else
	{
		std::cerr << "Error: Incorrcet splitting method name!";
		system("pause");
		exit(1);
	}
	return 0;
}
