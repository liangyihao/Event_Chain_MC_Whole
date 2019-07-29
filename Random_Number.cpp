/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
#include <random>
#include <iostream>
#include <limits>
#include "Random_Number.hpp"
using namespace std;
mt19937*generator;
uniform_real_distribution<double> UD(0.0,1.0);

void rand_init(int seed){
    generator=new mt19937(seed);
}

double Uniform_Random(){
    return UD(*generator);
}

double Exponential_Random(double rate){
	if(rate==0)return (numeric_limits<double>::max());
    double r=Uniform_Random();
    return -log(1-r)/rate;
}

Frequency_Generator::Frequency_Generator(vector<double>Frequency_List_Unnormalized){
    Frequency_Hist.clear();
    double Tot=0;
    for(int k=0;k<Frequency_List_Unnormalized.size();k++){
        Frequency_Hist.push_back(Frequency_List_Unnormalized[k]);
        Tot+=Frequency_List_Unnormalized[k];
    }

    ////////////////////////////////
    int2_double temp1,temp2;
    vector<int2_double>G;
    vector<int2_double>L;
    G.clear();L.clear();Bin.clear();
    for(int k=0;k<Frequency_Hist.size();k++){
        Frequency_Hist[k]/=Tot;
        temp1.integer1=k;
        temp1.integer2=-1;
        temp1.weight1=Frequency_Hist[k]*Frequency_Hist.size();
        if(temp1.weight1>=1){
            G.push_back(temp1);
        }
        else{
            L.push_back(temp1);
        }
    }
    while((!L.empty())&&(!G.empty())){
        temp1=L[L.size()-1];
        temp2=G[G.size()-1];
        temp1.integer2=temp2.integer1;
        temp2.weight1-=(1-temp1.weight1);
        Bin.push_back(temp1);
        L.pop_back();
        if(temp2.weight1<1){
            G.pop_back();
            L.push_back(temp2);
        }else{
            G[G.size()-1]=temp2;
        }
    }
    while(!G.empty()){
        G[G.size()-1].integer2=G[G.size()-1].integer1;
        if(abs(G[G.size()-1].weight1-1)>1E-15)cout<<"Warning"<<G[G.size()-1].weight1<<endl;
        G[G.size()-1].weight1=1;
        Bin.push_back(G[G.size()-1]);
        G.pop_back();
    }
    while(!L.empty()){
        L[L.size()-1].integer2=L[L.size()-1].integer1;
        if(abs(L[L.size()-1].weight1-1)>1E-15)cout<<"Warning"<<L[L.size()-1].weight1<<endl;
        L[L.size()-1].weight1=1;
        Bin.push_back(L[L.size()-1]);
        L.pop_back();
    }
	/*
    cout<<"Frequency Generator constructed"<<endl;
    vector<double>Prob(Frequency_Hist.size(),0);
    for(int k=0;k<Bin.size();k++){
        cout<<Bin[k].integer1<<'('<<Bin[k].weight1<<") "<<Bin[k].integer2<<'('<<(1.0-Bin[k].weight1)<<')'<<endl;
        Prob[Bin[k].integer1]+=Bin[k].weight1;
        Prob[Bin[k].integer2]+=(1.0-Bin[k].weight1);
    }
    cout<<endl<<"Expect Probability"<<endl;
    for(int k=0;k<Prob.size();k++)cout<<k<<" "<<Prob[k]/Bin.size()<<"(contructed)  "<<Frequency_Hist[k]<<"(Input)"<<endl;
    cout<<endl;
    */
}

int Frequency_Generator::Gen(){
    double r;
    r=Uniform_Random()*Bin.size();
    int bin_id;
    bin_id=(int)r;
    r-=bin_id;
    if(r>Bin[bin_id].weight1)return Bin[bin_id].integer2;else return Bin[bin_id].integer1;
    return -1;
}