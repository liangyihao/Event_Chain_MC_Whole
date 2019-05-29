#ifndef RANDOM_NUMBER_ECMC
#define RANDOM_NUMBER_ECMC
#include <vector>
using namespace std;
void rand_init(int seed);
double Uniform_Random();
double Exponential_Random(double rate);


typedef struct{
    int integer1,integer2;
    double weight1;
}int2_double;

class Frequency_Generator{
    private:
        vector<double>Frequency_Hist;
        vector<int2_double>Bin;
    public:
        Frequency_Generator(vector<double>Frequency_List_Unnormalized);
        int Gen();
};
#endif
