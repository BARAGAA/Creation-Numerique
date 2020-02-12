#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "matrix.hpp"

using namespace std;

int main(int argc, char** argv){

    cout << "Do or do not. There is no try." << endl;
        double *m1 = allocateMatrix(4,4);
        double *m2 = allocateMatrix(4,4);
        double *m3 = allocateMatrix(4,4);
        setMatrixZero (m3,4,4);

    return 0;
}