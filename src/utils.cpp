#include<cstdlib>
#include<fstream>
#include"utils.h"

Console::Console(int argc, const char* argv[]):m_seed(1000), m_trainFrac(0.6),
		m_testFrac(0.2), m_outTrain("train.txt"), m_outTest("test.txt"), m_outValid("valid.txt") {
	if (argc < 4) // Lack som necessary parameter: inputFileName; choice; m_method
		m_inputFileName = "";
	else  // Basic parameters are existed
	{
		m_inputFileName = argv[1]; 

		string cchoice = argv[2];
		if (cchoice == "T" || cchoice == "t" || cchoice == "TRUE" || cchoice == "True" || cchoice == "true")
			m_choice = 1;
		else if (cchoice == "F" || cchoice == "f" || cchoice == "FALSE" || cchoice == "False" || cchoice == "false")
			m_choice = 0;
		else
			m_choice = -1;

		string mmethod = argv[3];
		if (mmethod == "SOMPLEX")
			m_method = SOMPLEX;
		else if (mmethod == "MDUPLEX")
			m_method = MDUPLEX;
		else if (mmethod == "SBSS-P")
			m_method = SBSS_P;
		else if (mmethod == "SBSS-N")
			m_method = SBSS_N;
		else if (mmethod == "DUPLEX")
			m_method = DUPLEX;
		else if (mmethod == "SS")
			m_method = SS;
		else
			m_method = Wrong;
		
	}
	if (argc >= 10) { //If other arguments are specified, the default arguments are overridden
		m_seed = atoi(argv[4]);
		m_trainFrac = atof(argv[5]);
		m_testFrac = atof(argv[6]);
		m_outTrain = argv[7];
		m_outTest = argv[8];
		m_outValid = argv[9];
	}
}

void Console::checkValid() {
	if (m_inputFileName.size()==0)
		throw m_inputFileName;
	if (m_choice == -1)
		throw m_choice;
	if (m_method == Wrong)
		throw m_method;
	if (m_trainFrac + m_testFrac >= 1.0)
		throw m_trainFrac + m_testFrac;
}

void Console::check() {
	try {
		checkValid();
	}
	catch (string err) {
		std::ofstream outCerr;
		outCerr.open("log.out",std::ios::app);
		outCerr << "Lack some basic parameter: [inputFilename]: " << err << std::endl;
		outCerr.close();
		exit(1);
	}
	catch (int err) {
		std::ofstream outCerr;
		outCerr.open("log.out", std::ios::app);
		outCerr << "Lack some basic parameter: set [T,t,TRUE,True,true] or [F,f,FALSE,False,false] for whether include input vectors!" << std::endl;
		outCerr.close();
		exit(1);
	}
	catch (methodType err) {
		std::ofstream outCerr;
		outCerr.open("log.out", std::ios::app);
		outCerr << "Lack some basic parameter: [splitting method name]: " << methodTypeName[err] << std::endl;
		outCerr.close();
		exit(1);
	}
	catch (double err) {
		std::ofstream outCerr;
		outCerr.open("log.out", std::ios::app);
		outCerr << "Invalid fraction!: trainFrac + testFrac > 1.0: " << err << std::endl;
		outCerr.close();
		exit(1);
	}
}
