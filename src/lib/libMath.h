#ifndef _LIB_MATH_H_
#define _LIB_MATH_H_
#include <vector>
#include <cmath>

#define ABS(a) (a>=0 ? a:-a)
#define SIGN(a) (a>=0 ? 1:-1)
#define MIN(a,b) (a>b ? b:a)
#define MAX(a,b) (a>b ? a:b)

// probability 
double ran0(long *idum);
double NormalDistribution(long *seed_ptr);
double NormalDistributionPDF(double x,int mean,double var);
double NormalDistributionCDF(double x,int mean,double var);

// linear algebra
bool printMatrix(std::vector<std::vector<char> > &matrix);
bool matrixMultiplication(std::vector<char> &in, std::vector<std::vector<char> > &matrix, std::vector<char> &out);
bool matrixMultiplication(std::vector<std::vector<char> > &in, std::vector<std::vector<char> > &matrix, std::vector<std::vector<char> > &out);
unsigned int GaussianJordanElimination(std::vector<std::vector<char> > &matrix, int corner);
void findNullSpace(std::vector<std::vector<char> > inMatrix, std::vector<std::vector<char> > &outMatrix);
bool transposeMatrix(std::vector<std::vector<char> > &in, std::vector<std::vector<char> > &out);
bool copyMatrix(std::vector<std::vector<char> > &in, std::vector<std::vector<char> > &out);
bool permutationMatrix(std::vector<std::vector<char> > &matrix, unsigned int size, long *random_seed);
bool permutationMatrix(std::vector<std::vector<char> > &matrix, std::vector<unsigned int> &permutation_array);
   
// finite field, field poly x^3 + 2x^2 = [1, 2, 0, 0]
// subfunction about field contruct
bool isPrime(unsigned int in);
bool num2FieldPoly(unsigned int in, unsigned int prime_num, std::vector<unsigned int> &out);
unsigned int fieldPoly2Num(std::vector<unsigned int> &in, unsigned int prime_num);
bool creatFieldPoly(unsigned int prime_num, unsigned int dimension, std::vector<unsigned int> &out);
// field is prime
bool fieldPolyMod(std::vector<unsigned int> poly_A, std::vector<unsigned int> &poly_B, unsigned int prime_num, std::vector<unsigned int> &out);
// field is power fo prime
unsigned int fieldAdd(unsigned int A, unsigned int B);
unsigned int fieldMulti(unsigned int A, unsigned int B, std::vector<unsigned int> &AtoD, std::vector<unsigned int> &DtoA);
unsigned int fieldDivide(unsigned int A, unsigned int B, std::vector<unsigned int> &AtoD, std::vector<unsigned int> &DtoA);
bool fieldPolyMulti(std::vector<unsigned int> &poly_A, std::vector<unsigned int> &poly_B, std::vector<unsigned int> &out, std::vector<unsigned int> &AtoD, std::vector<unsigned int> &DtoA);
bool fieldPolyMod(std::vector<unsigned int> poly_A, std::vector<unsigned int> &poly_B, std::vector<unsigned int> &out, std::vector<unsigned int> &AtoD, std::vector<unsigned int> &DtoA);


unsigned int oneCount(std::vector<char> &in);
double EuclideanDistance(std::vector<char> &codeword, std::vector<double> &received);
#endif
