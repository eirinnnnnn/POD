#include "libMath.h"
#include "libDebug.h"
#include <cstdlib>
#include <cassert>
#include <cstdio>
#include <list>

// probability 
double ran0(long *idum){
    // Copy from Numerical Recipes in C  2nd edt  page 279
    // set idum to be any integer value
    // Input a seed , this sub-function will return a random number between [0 ,1]
    long IA = 16807;
    long IM = 2147483647;
    double AM = (1.0/IM);
    long IQ = 127773;
    long IR = 2836;
    long MASK = 123459876;

    long k ;
    double ans;
    *idum ^= MASK;
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum <0)  *idum+=IM;
    ans=(double)(AM*(*idum));
    *idum ^= MASK;
    return ans;
}

double NormalDistribution(long *seed_ptr){
    // Normal(0,1), mean = 0, var = 1
    double noise,s,v[2];
    do{
        v[0] = 2*ran0(seed_ptr) - 1;
        v[1] = 2*ran0(seed_ptr) - 1;
        s = v[0]*v[0] + v[1]*v[1];
    }while(s>=1);
    noise = v[0]*sqrt(-2*log(s)/s);
    return noise;
}

double NormalDistributionPDF(double x,int mean,double var){
    /* http://mathworld.wolfram.com/NormalDistribution.html */
    // printf("state : NormalDistributionPDF\n");
    return exp(pow(mean-x,2)/-2/var)/sqrt(var*2*3.1415926);
}

double NormalDistributionCDF(double x,int mean,double var){
    /* http://mathworld.wolfram.com/NormalDistribution.html */
    // printf("state : NormalDistributionCDF\n");
    return 0.5*(1+erf((x-mean)/sqrt(var*2)));
}

// linear algebra
bool printMatrix(std::vector<std::vector<char> > &matrix){
    for(unsigned int idx=0; idx<matrix.size(); idx++){
        for(unsigned int jdx=0; jdx<matrix[idx].size(); jdx++)
            printf(" %d",matrix[idx][jdx]);
        printf("\n");
    }
    return 0;
}

bool matrixMultiplication(std::vector<char> &in, 
                          std::vector<std::vector<char> > &matrix, 
                          std::vector<char> &out){
     
    assert(in.size() > 0);
    if(in.size() != matrix.size()){
        ERROR("in.size(%ld) != matrix.size(%ld)",in.size(), matrix.size()); exit(1);
    }
    assert(matrix[0].size() > 0);
    out.resize(matrix[0].size());

    for(unsigned int i_idx=0; i_idx<out.size(); i_idx++){
        out[i_idx] = 0;
        for(unsigned int j_idx=0; j_idx<in.size(); j_idx++)
            out[i_idx] ^= (in[j_idx] & matrix[j_idx][i_idx] );
    }
     
    return 0;
}

bool matrixMultiplication(std::vector<std::vector<char> > &in, 
                          std::vector<std::vector<char> > &matrix, 
                          std::vector<std::vector<char> > &out){

    // assert(in.size() > 0);
    // assert(in[0].size() > 0);
    // assert(in[0].size() == matrix.size());
    // assert(matrix[0].size() > 0);
    out.resize(in.size());

    for(unsigned int i_idx=0; i_idx<out.size(); i_idx++){
        out[i_idx].resize(matrix[0].size());
        for(unsigned int j_idx=0; j_idx<out[i_idx].size(); j_idx++){
            out[i_idx][j_idx] = 0;
            for(unsigned int m_idx=0; m_idx<in[0].size(); m_idx++)
                out[i_idx][j_idx] ^= (in[i_idx][m_idx] & matrix[m_idx][j_idx]);
        }
    }

    return 0;
}

unsigned int GaussianJordanElimination(std::vector<std::vector<char> > &matrix, int corner){

    unsigned int row_start, row_step, row_end;
    unsigned int col_start, col_step, col_end;

    if(matrix.size()==0){
        ERROR("matrix.size()==0");  exit(1);}
    // if(matrix[0].size() < matrix.size()){
        // ERROR("matrix[0].size(%d) < matrix.size(%d)",matrix[0].size(),matrix.size());  exit(1);}

    switch(corner){
        case 0:        // left top
            row_start=0; row_step=1; row_end=matrix.size()-1;
            col_start=0; col_step=1; col_end=matrix[0].size()-1;
            break;
        case 1:        // right top
            row_start=0; row_step=1; row_end=matrix.size()-1;
            col_start=matrix[0].size()-1; col_step=-1; col_end=0;
            break;
        case 2:        // left down
            row_start=matrix.size()-1; row_step=-1; row_end=0;
            col_start=0; col_step=1; col_end=matrix[0].size()-1;
            break;
        case 3:        // right down
            row_start=matrix.size()-1; row_step=-1; row_end=0;
            col_start=matrix[0].size()-1; col_step=-1; col_end=0;
            break;
        default:
            ERROR("corner fail : %d\n",corner);
            exit(1);
            return -1;
    }

    unsigned int rank=0;
    for(unsigned int col_idx=col_start; col_idx!=col_end+col_step; col_idx+=col_step){
        int target_row = -1;
        for(unsigned int row_idx = row_start + rank*row_step; row_idx!=row_end+row_step; row_idx+=row_step){
            if(matrix[row_idx][col_idx]){
                target_row = row_idx;
                for(unsigned int subcol_idx=col_start; subcol_idx!=col_idx; subcol_idx+=col_step){
                    if(matrix[row_idx][subcol_idx]){
                        target_row = -1;
                        break;
                    }
                }
                if(target_row != -1)
                    break;
            }
        }
        if(target_row == -1)
            continue;

        if(target_row != (int)(row_start + rank*row_step))
            for(unsigned int subcol_idx=col_idx; subcol_idx!=col_end+col_step; subcol_idx+=col_step)
                matrix[row_start + rank*row_step][subcol_idx] ^= matrix[target_row][subcol_idx];

        target_row = row_start + rank*row_step;

        for(unsigned int subrow_idx=row_start; subrow_idx!=row_end+row_step; subrow_idx+=row_step){
            if((int)subrow_idx==target_row)
                continue;
            if(matrix[subrow_idx][col_idx])
                for(unsigned int subcol_idx=col_idx; subcol_idx!=col_end+col_step; subcol_idx+=col_step)
                    matrix[subrow_idx][subcol_idx] ^= matrix[target_row][subcol_idx];
        }

        rank+=1;
    }
    return rank;
}

void findNullSpace(std::vector<std::vector<char> > inMatrix, std::vector<std::vector<char> > &outMatrix){
    // assert(inMatrix.size() <= inMatrix[0].size());
    unsigned int rank = GaussianJordanElimination(inMatrix, 0);
    unsigned int outMatrix_col = inMatrix[0].size()-rank;
    std::vector<unsigned int> leader(rank, 0), others(outMatrix_col, 0);
    INFO("rank %u outMatrix_col %u", rank ,outMatrix_col);
    for(unsigned col_idx=0, leader_idx=0, others_idx=0; col_idx<inMatrix[0].size(); col_idx++){
        if(leader_idx == rank){
            others[others_idx] = col_idx;
            others_idx++;
        }else if(inMatrix[leader_idx][col_idx]){
            leader[leader_idx] = col_idx;
            leader_idx++;
        }else{
            others[others_idx] = col_idx;
            others_idx++;
        }
    }
    outMatrix.resize(inMatrix[0].size());
    
    for(unsigned int row_idx=0, leader_idx=0, others_idx=0; row_idx<inMatrix[0].size(); row_idx++){
        outMatrix[row_idx].assign(outMatrix_col,0);
        if(row_idx == leader[leader_idx]){
            for(unsigned int col_idx=0; col_idx<outMatrix_col; col_idx++)
                outMatrix[row_idx][col_idx] = inMatrix[leader_idx][others[col_idx]];
            leader_idx++;
        }else{
            outMatrix[row_idx][others_idx] = 1;
            others_idx++;
        }
    }
}

bool transposeMatrix(std::vector<std::vector<char> > &in, std::vector<std::vector<char> > &out){
    out.resize(in[0].size());
    for(unsigned int idx=0; idx<in[0].size(); idx++){
        out[idx].resize(in.size());
        for(unsigned int jdx=0; jdx<in.size(); jdx++)
            out[idx][jdx] = in[jdx][idx];
    }
    return 0;
}
bool copyMatrix(std::vector<std::vector<char> > &in, std::vector<std::vector<char> > &out){
    out.resize(in.size());
    for(unsigned int idx=0; idx<in.size(); idx++){
        out[idx].resize(in[0].size());
        for(unsigned int jdx=0; jdx<in[0].size(); jdx++)
            out[idx][jdx] = in[idx][jdx];
    }
    return 0;
}

bool permutationMatrix(std::vector<std::vector<char> > &matrix, unsigned int size, long *random_seed){
    std::vector<unsigned int> site;
    site.resize(size);
    matrix.resize(size);
    for(unsigned int idx=0; idx<size; idx++){
        site[idx] = idx;
        matrix[idx].assign(size,0);
    }
    unsigned int loop_count = 0;
    while(loop_count<size*3){
        loop_count++;
        unsigned int A = (unsigned int)(ran0(random_seed)*size);
        unsigned int B = (unsigned int)(ran0(random_seed)*size);
        if(A!=B){
            site[A] ^= site[B];
            site[B] ^= site[A];
            site[A] ^= site[B];
        }
    }
    for(unsigned int idx=0; idx<size; idx++){
        matrix[idx].assign(size,0);
        matrix[idx][site[idx]] = 1;
    }
    return 0;
}

bool permutationMatrix(std::vector<std::vector<char> > &matrix, std::vector<unsigned int> &permutation_array){
    matrix.resize(permutation_array.size());
    for(unsigned int idx=0; idx<matrix.size(); idx++)
        matrix[idx].resize(permutation_array.size());
    for(unsigned int idx=0; idx<matrix.size(); idx++)
        matrix[permutation_array[idx]][idx] = 1;
    return 0;
}

// finite field, field poly x^3 + 2x^2 = [1, 2, 0, 0]
bool isPrime(unsigned int in){
    for(unsigned int idx=2; idx*idx<=in; idx++){
        if(isPrime(idx))
            if(in % idx == 0)
                return false;
    }
    return true;
}

bool num2FieldPoly(unsigned int in, unsigned int prime_num, std::vector<unsigned int> &out){
    std::list<unsigned int> ans;
    while(in > 0){
        ans.insert(ans.begin(),in%prime_num);
        in = in/prime_num;
    }
    out.assign(ans.begin(),ans.end());
    return 0;
}

unsigned int fieldPoly2Num(std::vector<unsigned int> &in, unsigned int prime_num){
    unsigned int ans=0;
    for(unsigned int idx=0; idx<in.size(); idx++){
        ans = ans*prime_num + in[idx];
    } 
    return ans;
}

bool creatFieldPoly(unsigned int prime_num, unsigned int dimension, std::vector<unsigned int> &out){
    std::vector<unsigned int> target_poly;
    std::vector<unsigned int> temp_polyA;
    std::vector<unsigned int> temp_polyB;

    for(unsigned int poly_num = (unsigned int)pow(prime_num,dimension)+1;
                     poly_num < (unsigned int)pow(prime_num,dimension+1)-1; poly_num++){
        bool irreducible = true;
        num2FieldPoly(poly_num,prime_num,target_poly);
        // check irreducible
        for(unsigned int sub_poly_num = 2; sub_poly_num < (unsigned int)pow(prime_num,dimension); sub_poly_num++){
            num2FieldPoly(sub_poly_num,prime_num,temp_polyA);
            fieldPolyMod(target_poly,temp_polyA,prime_num,temp_polyB);
            if(fieldPoly2Num(temp_polyB,prime_num) == 0){
                irreducible = false;
                break;
            }
        }
        // check primetive
        if(irreducible){
            for(int order=1; order < (int)pow(prime_num,dimension); order++){
                temp_polyA.assign(order,prime_num-1);
                fieldPolyMod(temp_polyA,target_poly,prime_num,temp_polyB);
                if(fieldPoly2Num(temp_polyB,prime_num) == 0){
                    if(order == (int)pow(prime_num,dimension)-1){
                        out.assign(target_poly.begin(),target_poly.end());
                        return 0;
                    }else{
                        break;
                    }
                }
            }
        }
    }
    return 0;
}

bool fieldPolyMod(std::vector<unsigned int> poly_A, std::vector<unsigned int> &poly_B, unsigned int prime_num, std::vector<unsigned int> &out){
    if(poly_B.size() == 0){
        ERROR("poly_B.size() == 0");    exit(1);
    }
    if(!isPrime(prime_num)){
        ERROR("!isPrime(prime_num)");    exit(1);
    }
    unsigned int poly_A_idx = 0;
    while(poly_A.size() - poly_A_idx >= poly_B.size()){
        if(poly_A[poly_A_idx] == 0){
            poly_A_idx++;
        }else{
            unsigned int order=0;
            for(unsigned int idx=0; idx<prime_num; idx++){
                if((idx*poly_B[0]) % prime_num == poly_A[poly_A_idx]){
                    order = idx;
                    break;
                }
            }
            for(unsigned int idx=0; idx<poly_B.size(); idx++)
                poly_A[idx+poly_A_idx] = (poly_A[idx+poly_A_idx] - (poly_B[idx] * order % prime_num) + prime_num) % prime_num;
            poly_A_idx++;
        }
    }
    while(poly_A.size()-poly_A_idx > 0){
        if(poly_A[poly_A_idx] == 0)
            poly_A_idx++;
        else
            break;
    }
    out.assign(poly_A.begin()+poly_A_idx,poly_A.end());
    return 0;
}

unsigned int fieldAdd(unsigned int A, unsigned int B){
    return A^B;
}

unsigned int fieldMulti(unsigned int A, unsigned int B, std::vector<unsigned int> &AtoD, std::vector<unsigned int> &DtoA){
    if(A==0 || B== 0)
        return 0;
    return DtoA[(AtoD[A] + AtoD[B]) % (AtoD.size()-1)];
}

unsigned int fieldDivide(unsigned int A, unsigned int B, std::vector<unsigned int> &AtoD, std::vector<unsigned int> &DtoA){
    if(B == 0){
        ERROR("B == 0");    exit(1);
    }else if(A==0){
        return 0;
    }
    return DtoA[(AtoD[A] + AtoD.size()-1 - AtoD[B]) %(AtoD.size())];
}

bool fieldPolyMulti(std::vector<unsigned int> &poly_A, std::vector<unsigned int> &poly_B, std::vector<unsigned int> &out, std::vector<unsigned int> &AtoD, std::vector<unsigned int> &DtoA){
    out.assign(poly_A.size()+poly_B.size()-1,0);
    for(unsigned int A_idx=0; A_idx<poly_A.size(); A_idx++)
        for(unsigned int B_idx=0; B_idx<poly_B.size(); B_idx++){
            out[A_idx+B_idx] = fieldAdd(out[A_idx+B_idx], fieldMulti(poly_A[A_idx],poly_B[B_idx],AtoD,DtoA));
        }
    return 0;
}

bool fieldPolyMod(std::vector<unsigned int> poly_A, std::vector<unsigned int> &poly_B, std::vector<unsigned int> &out, std::vector<unsigned int> &AtoD, std::vector<unsigned int> &DtoA){
    if(poly_B.size() == 0){
        ERROR("poly_B.size() == 0");    exit(1);
    }
    unsigned int poly_A_idx = 0;
    while(poly_A.size() - poly_A_idx >= poly_B.size()){
        if(poly_A[poly_A_idx] == 0){
            poly_A_idx++;
        }else{
            unsigned int order=fieldDivide(poly_A[poly_A_idx],poly_B[0],AtoD,DtoA);
            for(unsigned int idx=0; idx<poly_B.size(); idx++)
                poly_A[idx+poly_A_idx] = fieldAdd(poly_A[idx+poly_A_idx],fieldMulti(poly_B[idx],order,AtoD,DtoA));
            poly_A_idx++;
        }
    }
    while(poly_A.size()-poly_A_idx > 0){
        if(poly_A[poly_A_idx] == 0)
            poly_A_idx++;
        else
            break;
    }
    out.assign(poly_A.begin()+poly_A_idx,poly_A.end());
    return 0;
}

unsigned int oneCount(std::vector<char> &in){
     int ans = 0;
     for(unsigned int idx=0; idx<in.size(); idx++)
          if(in[idx] == 1)
               ans += 1;
     return ans;
}

double EuclideanDistance(std::vector<char> &codeword, std::vector<double> &received){
     assert(codeword.size() == received.size()) ;
     double ans = 0;
     for(unsigned int idx=0; idx<codeword.size(); idx++)
          ans += pow((codeword[idx] ? -1:1)-received[idx],2);
     return ans;
}
