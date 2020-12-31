// Compile it with
// g++ -O -I/path/to/coin/include -L/path/to/coin/lib -lCoinUtils -o lp2compact lp2compact.cpp
//
// then call with
//
// lp2compact < myfile.lp | bzip2 > myfile.bz2

#include "CoinFileIO.hpp"
#include "CoinLpIO.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVectorBase.hpp"

const char* d2s(const double d)
{
  static char str[1000];
  if (d > 1e29) {
    sprintf(str, "1e30");
  } else if (d < -1e29) {
    sprintf(str, "-1e30");
  } else {
    sprintf(str, "%f", d);
    int i;
    for (i = strlen(str) - 1; i > 0; --i) {
      if (str[i] != '0')
        break;
      str[i] = 0;
    }
    if (str[i] == '.')
      str[i] = 0;
  }
  return str;
}

int main(int argc, char* argv[])
{
  CoinLpIO lp;
  lp.readLp(stdin);
  const CoinPackedMatrix& m = *lp.getMatrixByRow();
  printf("%i\t%i\t%i\t0\n",
         m.getNumCols(), m.getNumRows(), m.getNumElements());
  const int ncol = lp.getNumCols();
  const double* clb = lp.getColLower();
  const double* cub = lp.getColUpper();
  for (int i = 0; i < ncol; ++i) {
    printf("%s", d2s(clb[i]));
    printf("\t%s", d2s(cub[i]));
    //		printf("\t%i", i);
    //		printf("\t%s", lp.getColNames()[i]);
    printf("\n");
  }
  const int nrow = lp.getNumRows();
  const double* rlb = lp.getRowLower();
  const double* rub = lp.getRowUpper();
  for (int i = 0; i < nrow; ++i) {
    const CoinPackedVectorBase& v = m.getVector(i);
    const int size = v.getNumElements();
    const int* ind = v.getIndices();
    const double* val = v.getElements();
    printf("%i", size);
    printf("\t%s", d2s(rlb[i]));
    printf("\t%s", d2s(rub[i]));
    for (int j = 0; j < size; ++j) {
      printf("\t%i\t%s", ind[j], d2s(val[j]));
    }
    printf("\n");
  }
  return 0;
}
