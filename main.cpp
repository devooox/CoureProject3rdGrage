#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>


class TMatrixView {
 public:
  int* M;
  bool type;
  int i;
  int j;
  int rowLength;
  int colLength;
  int n;
  int m;
  inline int get(int ind_i, int ind_j) {
    if (type) {
      return M[(i + ind_i) + (j+ind_j)*n];
    } else {
      return M[(i + ind_i)*m + (j + ind_j)];
    }
  }
  inline void set(int ind_i, int ind_j, int value) {
    if (type) {
      M[(i + ind_i) + (j+ind_j)*n] = value;
    } else {
      M[(i + ind_i)*m + (j + ind_j)] = value;
    }
  }
};

int* viewMultiplication(TMatrixView& A, TMatrixView& B, int n, int m, int l) {
  int* C = new int[n*l];
  for (int i = 0; i < n*l; ++i) {
    C[i] = 0;
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < l; ++j) {
      for (int k = 0; k < m; ++k) {
        C[i*l+j] += A.get(i, k)*B.get(k, j);
      }
    }
  }

  return C;
}

int* matrixMultiplication(const int* A, const int* B, int n, int m, int l) {
  int* C = new int[n*l];
  for (int i = 0; i < n*l; ++i) {
    C[i] = 0;
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < l; ++j) {
      for (int k = 0; k < m; ++k) {
        C[i*l+j] += A[i*m+k]*B[k+j*m];
      }
    }
  }
  return C;
}



int* blockMatrixMultiplication(TMatrixView& A, TMatrixView& B, int n, int m, int l, int border = 0) {
  int* C;
  if (m > border && border > 1) {
    C = new int[n*l];
    int newM = m/2;
    int newN = n/2;
    int newL = l/2;
    //BLOCK SPLIT
    TMatrixView A11;
    A11.M = A.M;
    A11.i = A.i;
    A11.j = A.j;
    A11.colLength = A.colLength/2;
    A11.rowLength = A.rowLength/2;
    A11.type = false;
    A11.n = A.n;
    A11.m = A.m;

    TMatrixView A21;
    A21.M = A.M;
    A21.i = A.i + A.colLength/2;
    A21.j = A.j;
    A21.colLength = A.colLength - A11.colLength;
    A21.rowLength = A.rowLength/2;
    A21.type = false;
    A21.n = A.n;
    A21.m = A.m;

    TMatrixView A12;
    A12.M = A.M;
    A12.i = A.i;
    A12.j = A.j + A.rowLength/2;
    A12.colLength = A.colLength/2;
    A12.rowLength = A.rowLength - A11.rowLength;
    A12.type = false;
    A12.n = A.n;
    A12.m = A.m;

    TMatrixView A22;
    A22.M = A.M;
    A22.i = A.i + A.rowLength/2;
    A22.j = A.j + A.colLength/2;
    A22.colLength = A.colLength - A11.colLength;
    A22.rowLength = A.rowLength - A11.rowLength;
    A22.type = false;
    A22.n = A.n;
    A22.m = A.m;


    TMatrixView B11;
    B11.M = B.M;
    B11.i = B.i;
    B11.j = B.j;
    B11.colLength = B.colLength/2;
    B11.rowLength = B.rowLength/2;
    B11.type = true;
    B11.n = B.n;
    B11.m = B.m;

    TMatrixView B12;
    B12.M = B.M;
    B12.i = B.i;
    B12.j = B.j + B.rowLength/2;
    B12.colLength = B.colLength/2;
    B12.rowLength = B.rowLength - B11.rowLength;
    B12.type = true;
    B12.n = B.n;
    B12.m = B.m;

    TMatrixView B21;
    B21.M = B.M;
    B21.i = B.i + B.colLength/2;
    B21.j = B.j;
    B21.colLength = B.colLength - B11.colLength;
    B21.rowLength = B.rowLength/2;
    B21.type = true;
    B21.n = B.n;
    B21.m = B.m;

    TMatrixView B22;
    B22.M = B.M;
    B22.i = B.i + B.colLength/2;
    B22.j = B.j + B.rowLength/2;
    B22.colLength = B.colLength - B11.colLength;
    B22.rowLength = B.rowLength - B11.rowLength;
    B22.type = true;
    B22.n = B.n;
    B22.m = B.m;

    //BLOCK MULTIPLICATION
    int* tmpMatrix1 = blockMatrixMultiplication(A11, B11, newN, newM, newL, border);
    int* tmpMatrix2 = blockMatrixMultiplication(A12, B21, newN, m-newM, newL, border);
    for (int i = 0; i < newN*newL; ++i) {
      tmpMatrix1[i] += tmpMatrix2[i];
    }
    for (int i = 0; i < newN; ++i) {
      for (int j = 0; j < newL; ++j) {
        C[i*l + j] = tmpMatrix1[i*newL + j];
      }
    }

    tmpMatrix1 = blockMatrixMultiplication(A11, B12, newN, newM, l-newL, border);
    tmpMatrix2 = blockMatrixMultiplication(A12, B22, newN, m-newM, l-newL, border);
    for (int i = 0; i < newN*(l-newL); ++i) {
      tmpMatrix1[i] += tmpMatrix2[i];
    }
    for (int i = 0; i < newN; ++i) {
      for (int j = newL; j < l; ++j) {
        C[i*l+j] = tmpMatrix1[i*(l-newL) + j-newL];
      }
    }

    tmpMatrix1 = blockMatrixMultiplication(A21, B11, n-newN, newM, newL, border);
    tmpMatrix2 = blockMatrixMultiplication(A22, B21, n-newN, m-newM, newL, border);
    for (int i = 0; i < (n-newN)*newL; ++i) {
      tmpMatrix1[i] += tmpMatrix2[i];
    }
    for (int i = newN; i < n; ++i) {
      for (int j = 0; j < newL; ++j) {
        C[i*l+j] = tmpMatrix1[(i-newN)*newL + j];
      }
    }

    tmpMatrix1 = blockMatrixMultiplication(A21, B12, n-newN, newM, l-newL, border);
    tmpMatrix2 = blockMatrixMultiplication(A22, B22, n-newN, m-newM, l-newL, border);
    for (int i = 0; i < (n-newN)*(l-newL); ++i) {
      tmpMatrix1[i] += tmpMatrix2[i];
    }
    for (int i = newN; i < n; ++i) {
      for (int j = newL; j < l; ++j) {
        C[i*l+j] = tmpMatrix1[(i-newN)*(l-newL) + j-newL];
      }
    }

  } else {
    C = viewMultiplication(A, B, n, m, l);
  }

  return C;
}


int* blockMatrixMultiplicationPrepare(int *A, int *B, int n, int m, int l, int border = 0) {
  TMatrixView a;
  a.M = A;
  a.i = 0;
  a.j = 0;
  a.colLength = n;
  a.rowLength = m;
  a.type = 0;
  a.n = n;
  a.m = m;
  TMatrixView b;
  b.M = B;
  b.i = 0;
  b.j = 0;
  b.colLength = m;
  b.rowLength = l;
  b.type = 1;
  b.n = m;
  b.m = l;

  return blockMatrixMultiplication(a, b, n, m, l, border);
}



int* shtrassenMatrixMultiplication(TMatrixView& A, TMatrixView& B, int size, int border = 0) {

  int* C;
  if (size > border && border > 1) {
    int* D;
    int* D1;
    int* D2;
    int* H1;
    int* H2;
    int* V1;
    int* V2;
    TMatrixView tmpMatrix1;
    tmpMatrix1.type = 0;
    tmpMatrix1.rowLength = A.rowLength/2;
    tmpMatrix1.colLength = A.colLength/2;
    tmpMatrix1.i = 0;
    tmpMatrix1.j = 0;
    tmpMatrix1.n = A.rowLength/2;
    tmpMatrix1.m = A.colLength/2;
    tmpMatrix1.M = new int[tmpMatrix1.n*tmpMatrix1.m];
    TMatrixView tmpMatrix2;
    tmpMatrix2.type = 1;
    tmpMatrix2.rowLength = B.rowLength/2;
    tmpMatrix2.colLength = B.colLength/2;
    tmpMatrix2.i = 0;
    tmpMatrix2.j = 0;
    tmpMatrix2.n = B.rowLength/2;
    tmpMatrix2.m = B.colLength/2;
    tmpMatrix2.M = new int[tmpMatrix2.n*tmpMatrix2.m];
    for (int i = 0; i < A.rowLength/2; ++i) {
      for (int j = 0; j < A.colLength/2; ++j) {
        tmpMatrix1.set(i, j, A.get(i, j) + A.get(i + A.colLength/2, j + A.rowLength/2));
        tmpMatrix2.set(i, j, B.get(i, j) + B.get(i + B.colLength/2, j + B.rowLength/2));
      }
    }
    D = shtrassenMatrixMultiplication(tmpMatrix1, tmpMatrix2, size/2, border);

    for (int i = 0; i < A.rowLength/2; ++i) {
      for (int j = 0; j < A.colLength/2; ++j) {
        tmpMatrix1.set(i, j, A.get(i, j + A.rowLength/2) - A.get(i + A.colLength/2, j + A.rowLength/2));
        tmpMatrix2.set(i, j, B.get(i + A.colLength, j) + B.get(i + B.colLength/2, j + B.rowLength/2));
      }
    }
    D1 = shtrassenMatrixMultiplication(tmpMatrix1, tmpMatrix2, size/2, border);

    for (int i = 0; i < A.rowLength/2; ++i) {
      for (int j = 0; j < A.colLength/2; ++j) {
        tmpMatrix1.set(i, j, A.get(i + A.colLength/2, j) - A.get(i, j));
        tmpMatrix2.set(i, j, B.get(i + A.colLength, j) + B.get(i, j + B.rowLength/2));
      }
    }
    D2 = shtrassenMatrixMultiplication(tmpMatrix1, tmpMatrix2, size/2, border);

    for (int i = 0; i < A.rowLength/2; ++i) {
      for (int j = 0; j < A.colLength/2; ++j) {
        tmpMatrix1.set(i, j, A.get(i, j) + A.get(i, j + A.rowLength/2));
        tmpMatrix2.set(i, j, B.get(i + B.colLength/2, j + B.rowLength/2));
      }
    }
    H1 = shtrassenMatrixMultiplication(tmpMatrix1, tmpMatrix2, size/2, border);

    for (int i = 0; i < A.rowLength/2; ++i) {
      for (int j = 0; j < A.colLength/2; ++j) {
        tmpMatrix1.set(i, j, A.get(i + A.colLength/2, j) + A.get(i + A.colLength/2, j + A.rowLength/2));
        tmpMatrix2.set(i, j, B.get(i, j));
      }
    }
    H2 = shtrassenMatrixMultiplication(tmpMatrix1, tmpMatrix2, size/2, border);

    for (int i = 0; i < A.rowLength/2; ++i) {
      for (int j = 0; j < A.colLength/2; ++j) {
        tmpMatrix1.set(i, j, A.get(i + A.colLength/2, j + A.rowLength/2));
        tmpMatrix2.set(i, j, B.get(i + B.colLength/2, j) - B.get(i, j));
      }
    }
    V1 = shtrassenMatrixMultiplication(tmpMatrix1, tmpMatrix2, size/2, border);

    for (int i = 0; i < A.rowLength/2; ++i) {
      for (int j = 0; j < A.colLength/2; ++j) {
        tmpMatrix1.set(i, j, A.get(i, j));
        tmpMatrix2.set(i, j, B.get(i, j + B.rowLength/2) - B.get(i + B.colLength/2, j + B.rowLength/2));
      }
    }
    V2 = shtrassenMatrixMultiplication(tmpMatrix1, tmpMatrix2, size/2, border);

    int* ans = new int[size*size];
    for (int i = 0; i < size*size; ++i) {
      ans[i] = 0;
    }

    for (int i = 0; i < size/2; ++i) {
      for (int j = 0; j < size/2; ++j) {
        ans[i*size + j] = D[i*size/2 + j] + D1[i*size/2 + j] + V1[i*size/2 + j] - H1[i*size/2 + j];
        ans[i*size + j + size/2] = V2[i*size/2 + j] + H1[i*size/2 + j];
        ans[(i+size/2)*size + j] = V1[i*size/2 + j] + H2[i*size/2 + j];
        ans[(i+size/2)*size + j + size/2] = D[i*size/2 + j] + D2[i*size/2 + j] + V2[i*size/2 + j] - H2[i*size/2 + j];
      }
    }
    return ans;
  } else {
    C = viewMultiplication(A, B, A.n, A.m, B.m);
//    C = matrixMultiplication(A.M, B.M, A.n, A.m, B.m);
  }

  return C;
}

int* shtrassenPrepare(int* A, int* B, int n, int m, int l, int border = 0) {
  int twoPower = 2;
  while (twoPower < n || twoPower < m || twoPower < l) {
    twoPower *= 2;
  }
  int* newA = new int[twoPower*twoPower];
  int* newB = new int[twoPower*twoPower];
  for (int i = 0; i < twoPower*twoPower; ++i) {
    newA[i] = 0;
    newB[i] = 0;
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      newA[i*twoPower + j] = A[i*m + j];
    }
  }
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < l; ++j) {
      newB[i + j*twoPower] = B[i + j*m];
    }
  }

  TMatrixView a;
  a.M = newA;
  a.i = 0;
  a.j = 0;
  a.colLength = twoPower;
  a.rowLength = twoPower;
  a.type = 0;
  a.n = twoPower;
  a.m = twoPower;
  TMatrixView b;
  b.M = newB;
  b.i = 0;
  b.j = 0;
  b.colLength = twoPower;
  b.rowLength = twoPower;
  b.type = 1;
  b.n = twoPower;
  b.m = twoPower;


  int* Ans = shtrassenMatrixMultiplication(a, b, twoPower, border);

  int* C = new int[n*l];
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < l; ++j) {
      C[i*l + j] = Ans[i*twoPower + j];
    }
  }
  return C;
}

int main() {
//  freopen("output.txt", "w", stdout);
//  INITIALIZATION
  int n, m, l;
  std::ifstream fin("code3.txt");
  fin >> n >> m >> l;
  int *A = new int[n * m];
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      fin >> A[i * m + j];
    }
  }

  int *B = new int[m * l];
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < l; ++j) {
      fin >> B[i + j * m];
    }
  }
  TMatrixView a;
  a.rowLength = n;
  a.colLength = m;
  a.i = 0;
  a.j = 0;
  a.n = n;
  a.m = m;
  a.type = false;
  a.M = A;

  TMatrixView b;
  b.rowLength = m;
  b.colLength = l;
  b.i = 0;
  b.j = 0;
  b.n = m;
  b.m = l;
  b.type = true;
  b.M = B;
  int* Ans;
  int time = std::clock();

//  Ans = viewMultiplication(a, b, n, m, l);
//  Ans = matrixMultiplication(A, B, n, m, l);
//  Ans = blockMatrixMultiplicationPrepare(A, B, n, m, l, 512);
//  Ans = shtrassenPrepare(A, B, n, m, l, 512);
//  for (int i = 0; i < n; ++i) {
//    for (int j = 0; j < l; ++j) {
//      std::cout << Ans[i*l+j] << ' ';
//    }
//    std::cout << "\n";
//  }
  std::cout << clock() - time;
  return 0;
}

//
//#include <fstream>
//
//int main() {
//  std::ofstream fout("code3.txt");
//  fout << "1024 2 1024\n";
//  for (int i = 0; i < 1024; ++i) {
//    for (int j = 0; j < 2; ++j) {
//      fout << "1 ";
//    }
//    fout << '\n';
//  }
//  for (int i = 0; i < 2; ++i) {
//    for (int j = 0; j < 1024; ++j) {
//      fout << "1 ";
//    }
//    fout << '\n';
//  }
//  return 0;
//}