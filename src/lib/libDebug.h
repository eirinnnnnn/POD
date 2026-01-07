#ifndef _LIB_DEBUG_H_
#define _LIB_DEBUG_H_
#include <vector>
#include <iostream>
#define ERROR(format, args...) 											\
	do{ 																\
		printf("ERROR[%s:%d] " format"\n", __FILE__, __LINE__, ##args); \
		exit(1);														\
	} while(0)
#	ifndef NDEBUG
#		define DEBUG_MODE 1
#		define INFO(format, args...) printf("[INFO] " format"\n", ##args)
#		define WARRING(format, args...) printf("WARRING[%s:%d] " format"\n", __FILE__, __LINE__, ##args)
#	else
#		define DEBUG_MODE 0
#		define INFO(args...)
#		define WARRING(args...)
#	endif
#endif // end of _LIB_DEBUG_H_