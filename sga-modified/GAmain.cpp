/***************************************************************************
 *   Copyright (C) 2004 by Tian-Li Yu                                      *
 *   tianliyu@cc.ee.ntu.edu.tw                                             *
 *                                                                         *
 *   You can redistribute it and/or modify it as you like                  *
 ***************************************************************************/

#include <cmath>
#include <cstdio>
#include <iostream>
#include <cstdlib>

#include "statistics.h"
#include "ga.h"
#include "chromosome.h"
#include "global.h"

using namespace std;

int main (int argc, char *argv[])
{

    if (argc != 10 && argc != 11 && argc != 12 && argc != 13) {
        printf ("GA function ell nInitial selectionPressure pc pm maxGen maxFe repeat XO_type display rand_seed\n");
        printf ("                     [XO_type]: \n");
        printf ("                  one-point XO:  0\n");
        printf ("                    uniform XO:  1\n");
        printf ("     population-wise shuffling:  2\n");
        printf ("\n");
        printf ("                    [function]: \n");
        printf ("                        ONEMAX:  0\n");
        printf ("                        MK    :  1\n");
        printf ("                        FTRAP :  2\n");
        printf ("                        CYC   :  3\n");
        printf ("                        HTRAP :  4\n");
        return -1;
    }
    int fffff = atoi (argv[1]);

    int ell = atoi (argv[2]);    // problem size
                                 // initial population size
    int nInitial = atoi (argv[3]);
                                 // selection pressure
    int selectionPressure = atoi (argv[4]);
    double pc = atof (argv[5]);  // pc
    double pm = atof (argv[6]);  // pm
    int maxGen = atoi (argv[7]); // max generation
    int maxFe = atoi (argv[8]);  // max fe
    int repeat = atoi (argv[9]); // how many time to repeat
    int XO_type = 0;
    int display = 0;
    int rand_seed = -1;
    if (argc > 10) XO_type = atoi (argv[10]);
    if (argc > 11) display = atoi (argv[11]); // display each generation or not
    if (argc > 12) rand_seed = atoi (argv[12]);  // rand seed

    if (rand_seed != -1)  // time
        myRand.seed((unsigned long)rand_seed);

    int i;

    Statistics stGenS, stGenF;
    int usedGen;

    int failNum = 0;

    for (i = 0; i < repeat; i++) {

        GA ga (ell, nInitial, selectionPressure, pc, pm, maxGen, maxFe, XO_type, fffff);

        if (display == 0)
            usedGen = ga.doIt (false);
        else
            usedGen = ga.doIt (true);

        Chromosome ch(ell);
        if (ga.stFitness.getMax() < ch.getMaxFitness()) {
            printf ("-");
            failNum++;
            stGenF.record (usedGen);
        }
        else {
            printf ("+");
            stGenS.record (usedGen);
        }

        fflush (NULL);

    }

    printf ("\nAverage Gen of Success: %f\n", stGenS.getMean());
    printf ("Average Gen of Failure: %f\n", stGenF.getMean());
    printf ("Total number of Failure: %d\n", failNum);

    return EXIT_SUCCESS;
}
