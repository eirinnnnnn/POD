#ifndef _LIB_PARSER_H_
#define _LIB_PARSER_H_
#include <string>
#include <vector>
#include <map>

char parseBinaryMatrix(std::string &matrixPath, std::vector<std::vector<char> > &matrix);

char parseConfig(std::string config_path, std::map<std::string, std::map<std::string, std::string> > &config);
char writeBackConfig(std::string config_path, std::map<std::string, std::map<std::string, std::string> > &config);

#endif