//
//  main.cpp
//  LMSolver
//
//  Created by Claude Rogers on 7/8/11.
//  Copyright 2011 California Institute of Technology. All rights reserved.
//

#include "LMSolver.h"

int main (int argc, const char * argv[])
{

    char *myf = (char *) "/Users/cjrogers/Desktop/curvefit/cprog/examples/ic50_ex.txt";
    int myMod = ic50;
    double v1 = 1.0;
    double v2 = 1.0;
    double v3 = 1.0;
    double v4 = 0.0;
    
    LMSolver solver (myf, myMod, v1, v2, v3, v4);
    solver.levenberg_marquardt();
    return 0;
}

