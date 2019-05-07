#include <random>
using namespace std;
mt19937*generator;
uniform_real_distribution<double> UD(0.0,1.0);
void rand_init(int seed){
    generator=new mt19937(seed);
}

double Uniform_Random(){
    return UD(*generator);
}
