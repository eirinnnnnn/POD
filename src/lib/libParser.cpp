#include "libParser.h"
#include "libDebug.h"
#include <cassert>
#include <string>
#include <vector>


char parseBinaryMatrix(std::string &matrixPath, std::vector<std::vector<char> > &matrix){
	FILE *matrix_data = fopen(matrixPath.c_str(),"r");
	unsigned int row, col, tmpInt;
	if(!matrix_data){
		INFO("%s is not exist.", matrixPath.c_str());
		return(1);
	}
	fscanf(matrix_data,"%d",&row);
	fscanf(matrix_data,"%d",&col);

	matrix.resize(row);
	for(unsigned int i_idx=0;i_idx<row;i_idx++){
		matrix[i_idx].resize(col);
		for(unsigned int j_idx=0;j_idx<col;j_idx++){
			assert(!feof(matrix_data));
			fscanf(matrix_data,"%1d",&(tmpInt));
			matrix[i_idx][j_idx] = tmpInt;
		}
	}
	fclose(matrix_data);	
	return 0;
}

char parseConfig(std::string config_path, std::map<std::string, std::map<std::string, std::string> > &config){
	FILE *config_data = fopen(config_path.c_str(),"r");
	if(!config_data){
		INFO("%s is not exist.", config_path.c_str());
		return(1);
	}
	INFO("set config with %s\n",config_path.c_str());
	char line[20000];
	char key_name[10000];
	char value_name[10000];
	char section_name[10000];
	std::string target_name;
	while (fgets(line, sizeof(line), config_data)) {
        target_name = "key";
        for(int line_idx=0,name_idx=0; line_idx<20000; line_idx++){
        	if      (line[line_idx] == '['){
        		target_name = "section";
        		continue;
        	}else if(line[line_idx] == ']'){
        		section_name[name_idx] = '\0';
        		continue;
        	}else if(line[line_idx] == '='){
        		target_name = "value";
        		key_name[name_idx] = '\0';
        		name_idx = 0;
        		continue;
        	}else if(line[line_idx] == ' ' || line[line_idx] == '\n'){
        		continue;
        	}else{
        		if(target_name == "section")	section_name[name_idx] = line[line_idx];
        		if(target_name == "key")        key_name[name_idx] = line[line_idx];
        		if(target_name == "value")      value_name[name_idx] = line[line_idx];
        		name_idx++;
        		if(line[line_idx] == '\0')      break;
        	}
        }
        if(target_name == "value"){
        	if(config.find(section_name) != config.end())
        		if(config.find(section_name)->second.find(key_name) != config.find(section_name)->second.end())
        			config.find(section_name)->second.find(key_name)->second = value_name;
        		else
        			ERROR("key_name '%s' fail.",key_name);
        	else{
        		printf("%s\n", section_name);
        		ERROR("section_name '%s' fail.",section_name);
        	}
        }
    }	
	fclose(config_data);
	return 0;
}


char writeBackConfig(std::string config_path, std::map<std::string, std::map<std::string, std::string> > &config){
	FILE *config_data = fopen(config_path.c_str(),"w");
	assert(config_data);
	for(std::map<std::string, std::map<std::string, std::string> >::iterator pair_idx = config.begin();
	    pair_idx != config.end(); pair_idx++){
		fprintf(config_data,"[%s]\n",pair_idx->first.c_str());
		for(std::map<std::string, std::string>::iterator sub_pair_idx = pair_idx->second.begin(); 
			sub_pair_idx != pair_idx->second.end(); sub_pair_idx++){
			fprintf(config_data,"%-30s = %s\n",sub_pair_idx->first.c_str(),sub_pair_idx->second.c_str());
		}
	}
	fclose(config_data);
	return 0;
}