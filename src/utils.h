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
class Console {
public:
	Console();
	Console(int argc, const char* argv[]);
	void check();
	bool checkValid();

	// interface
	string getInputFileName() { return inputFileName; }
	int getChoice() { return choice; }
	methodType getMethod() { return method; }
	int	getSeed() { return seed; }
	float getTrainFrac() { return trainFrac; }
	float getTestFrac() { return testFrac; }
	string getOutTrain() { return outTrain; }
	string getOutTest() { return outTest; }
	string getOutVaild() { return outValid; }

private:
	string inputFileName;	// The name of input file, must be given
	int choice;	// Whether the input in included in the process of data splitting, must be given
	methodType method;	// Choose data splitting method, must be given
	int seed;	// Set seed, default = 1000
	float trainFrac;	// Set fraction for training set, default = 0.6
	float testFrac;	// Set fraction for test set, default = 0.2
	string outTrain;	// The name of output file, default = train.txt
	string outTest;		// default = test.txt
	string outValid;	// default = valid.txt
};

#endif // !_UTIL_S

