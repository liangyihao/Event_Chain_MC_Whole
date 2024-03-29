#ifndef INPUT_FILE_PARSER_HPP
#define INPUT_FILE_PARSER_HPP
#include <vector>
#include <string>
#include "public.hpp"
using namespace std;

extern vector<Instruction>Instruction_list;
extern int loop_times;
extern vector<int2>IO_Order_of_Beads;//Record the order of beads in the input file, output in this order
int Input_File_Parser(const char* fileName);
void xml_write(const char* fileName, const int timestep);
void next_input_file_writer(const char* fileName);
void Hard_Repulsion_Checker();
#endif
