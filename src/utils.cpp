#include<cstdlib>
#include<fstream>
#include"utils.h"

Console::Console(){
}

Console::Console(int argc, const char* argv[]):seed(1000), trainFrac(0.6),
		testFrac(0.2), outTrain("train.txt"), outTest("test.txt"), outValid("valid.txt") {
	if (argc < 4) // Lack som necessary parameter: inputFileName; choice; method
		inputFileName = "";
	else  // Basic parameters are existed
	{
		inputFileName = argv[1]; 

		string cchoice = argv[2];
		if (argv[2] == "true")
		choice = 1;
		else if (argv[2] == "false")
			choice = 0;
		else
			choice = -1;

		string mmethod = argv[3];
		if (mmethod == "SOMPLEX")
			method = SOMPLEX;
		else if (mmethod == "MDUPLEX")
			method = MDUPLEX;
		else if (mmethod == "SBSS-P")
			method = SBSS_P;
		else if (mmethod == "SBSS-N")
			method = SBSS_N;
		else if (mmethod == "DUPLEX")
			method = DUPLEX;
		else if (mmethod == "SS")
			method = SS;
		else
			method = Wrong;
		
	}
	if (argc >= 10) { //If other arguments are specified, the default arguments are overridden
		seed = atoi(argv[4]);
		trainFrac = atof(argv[5]);
		testFrac = atof(argv[6]);
		outTrain = argv[7];
		outTest = argv[8];
		outValid = argv[9];
	}
}

bool Console::checkValid() {
	if (inputFileName.size()==0)
		throw inputFileName;
	if (choice == -1)
		throw choice;
	if (method == Wrong)
		throw method;
	if (trainFrac + testFrac >= 1.0)
		throw trainFrac + testFrac;
	return true;
}

void Console::check() {
	try {
		checkValid();
	}
	catch (string) {
		std::ofstream outCerr;
		outCerr.open("log.out",std::ios::app);
		outCerr << "Lack some basic parameter: inputFilename!" << std::endl;
		outCerr.close();
		exit(1);
	}
	catch (int) {
		std::ofstream outCerr;
		outCerr.open("log.out", std::ios::app);
		outCerr << "Lack some basic parameter: [true] or [false] for IO!" << std::endl;
		outCerr.close();
		exit(1);
	}
	catch (methodType) {
		std::ofstream outCerr;
		outCerr.open("log.out", std::ios::app);
		outCerr << "Lack some basic parameter: splitting method!" << std::endl;
		outCerr.close();
		exit(1);
	}
	catch (float) {
		std::ofstream outCerr;
		outCerr.open("log.out", std::ios::app);
		outCerr << "Invalid fraction!: trainFrac + testFrac > 1.0!" << std::endl;
		outCerr.close();
		exit(1);
	}
}
