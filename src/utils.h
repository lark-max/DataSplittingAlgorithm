/*
* This file implements terminal interaction and exception handling classes
* Written by Chen J. Y. at 2020, Zhejiang UNV
*/
#ifndef _UTIL_S
#define _UTIL_S
#include<string>
using std::string;

enum methodType {
	SOMPLEX,
	SBSS_P,
	SBSS_N,
	MDUPLEX,
	DUPLEX,
	SS,
	Wrong
};

static const string methodTypeName[] = { "SOMPLEX", "SBSS-P", "SBSS-N", "MDUPLEX", "DUPLEX", "SS", "Wrong" };

class Console {
public:
	Console() = default;
	Console(int argc, const char* argv[]);
	void check();
	void checkValid();

	// interface
	string getInputFileName() { return m_inputFileName; }
	int getChoice() { return m_choice; }
	methodType getMethod() { return m_method; }
	int	getSeed() { return m_seed; }
	double getTrainFrac() { return m_trainFrac; }
	double getTestFrac() { return m_testFrac; }
	string getOutTrain() { return m_outTrain; }
	string getOutTest() { return m_outTest; }
	string getOutVaild() { return m_outValid; }

private:
	string m_inputFileName;	// The name of input file, must be given
	int m_choice;	// Whether the input in included in the process of data splitting, must be given
	methodType m_method;	// Choose data splitting method, must be given
	int m_seed;	// Set seed, default = 1000
	double m_trainFrac;	// Set fraction for training set, default = 0.6
	double m_testFrac;	// Set fraction for test set, default = 0.2
	string m_outTrain;	// The name of output file, default = train.txt
	string m_outTest;		// default = test.txt
	string m_outValid;	// default = valid.txt
};

#endif // !_UTIL_S

