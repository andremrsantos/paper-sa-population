#include <Rcpp.h>
#include <iostream>
#include <fstream>

#define PACK_DENSITY 4
/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
#define MASK0 3   /* 3 << 2 * 0 */
#define MASK1 12  /* 3 << 2 * 1 */
#define MASK2 48  /* 3 << 2 * 2 */
#define MASK3 192 /* 3 << 2 * 3 */

void decode_plink(
    const unsigned char *in,
    const uint n,
    Rcpp::IntegerMatrix::Column X)
{
    uint N = X.size();

    uint i, k;
    unsigned char tmp, geno;
    uint a1, a2;

    for (i = 0; i < n; ++i)
    {
        tmp = in[i];
        k = PACK_DENSITY * i;

        geno = (tmp & MASK0);
        a1 = !(geno & 1);
        a2 = !(geno >> 1);
        if (k < N)
            X[k] = ((geno == 1) ? NA_INTEGER : a1 + a2);
        k++;

        geno = (tmp & MASK1) >> 2;
        a1 = !(geno & 1);
        a2 = !(geno >> 1);
        if (k < N)
            X[k] = ((geno == 1) ? NA_INTEGER : a1 + a2);
        k++;

        geno = (tmp & MASK2) >> 4;
        a1 = !(geno & 1);
        a2 = !(geno >> 1);
        if (k < N)
            X[k] = ((geno == 1) ? NA_INTEGER : a1 + a2);
        k++;

        geno = (tmp & MASK3) >> 6;
        a1 = !(geno & 1);
        a2 = !(geno >> 1);
        if (k < N)
            X[k] = ((geno == 1) ? NA_INTEGER : a1 + a2);
    }
}

//[[Rcpp::export]]
Rcpp::IntegerMatrix read_bed(const std::string bedfile, const uint N)
{
    std::ifstream in(bedfile.c_str(), std::ios::in | std::ios::binary);
    if (!in) {
        throw std::runtime_error("io erorr");
    }
    // Compute number of variants
    in.seekg(0, std::ifstream::end);
    uint len = (unsigned int) in.tellg() - 3;
    uint np = (unsigned int)ceil((double) N / PACK_DENSITY);
    uint nsnps = len / np;

    Rcpp::IntegerMatrix X(N, nsnps);
    unsigned char *tmp = new unsigned char[np * PACK_DENSITY];

    in.seekg(3, std::ifstream::beg);
    for (uint j = 0; j < nsnps; j++)
    {
        in.read((char*) tmp, sizeof(char) * np);
        decode_plink(tmp, np, X.column(j));
    }

    delete[] tmp;
    in.close();


    return(X);
}