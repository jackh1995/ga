/***************************************************************************
 *   Copyright (C) 2004 by Tian-Li Yu                                      *
 *   tianliyu@cc.ee.ntu.edu.tw                                             *
 ***************************************************************************/

#include <stdio.h>
#include "global.h"
#include "chromosome.h"

#define TRAP_K 5
#define HTRAP_K 3

Chromosome::Chromosome ()
{
    length = 0;
    gene = NULL;
    evaluated = false;
}


Chromosome::Chromosome (int n_length)
{
    gene = NULL;
    init (n_length);
}


Chromosome::~Chromosome ()
{
    delete[]gene;
}


void Chromosome::init (int n_length)
{
    length = n_length;

    if (gene != NULL)
        delete[]gene;

    gene = new bool[length];
    evaluated = false;
}

int Chromosome::getVal (int index) const
{
    if (index < 0 || index > length)
        outputErrMsg ("Index overrange in Chromosome::operator[]");

    return (gene[index])? 1:0;
}


void Chromosome::setVal (int index, int val)
{
    if (index < 0 || index > length)
        outputErrMsg ("Index overrange in Chromosome::operator[]");

    gene[index] = (val==1)? true:false;
    evaluated = false;
}


double Chromosome::getFitness ()
{
    if (evaluated)
        return fitness;
    else
        return (fitness = evaluate ());
}


bool Chromosome::isEvaluated () const
{
    return evaluated;
}


double Chromosome::evaluate ()
{
    evaluated = true;


    double accum = 0.0;

    switch (FITNESS_FUNCTION) {
        case 0:
            accum = oneMax();
            break;
        case 1:
            accum = mkTrap(1, 0.8);
            break;
        case 3:
            accum = cycTrap(1, 0.8);
            break;
        case 2:
            accum = fTrap();
            break;
        case 4:
            accum = hTrap(1, 0.8);
            break;
        default:
            outputErrMsg ("fitness function does not defined");
            break;
    }

    return accum;
}


// OneMax
double Chromosome::oneMax () const {

    double result = 0;

    for (int i = 0; i < length; ++i)
        result += getVal(i);

    return result;
}


double Chromosome::trap (int unitary, double fHigh, double fLow, int trapK) const {
    if (unitary > trapK)
        return 0;

    if (unitary == trapK)
        return fHigh;
    else
        return fLow - unitary * fLow / (trapK-1);
}


double Chromosome::fTrap() const {

    double result = 0.0;

    for (int i=0; i<length/6; ++i) {
        int u=0;
        for (int j=0; j<6; ++j)
            u += getVal(i*6+j);

        if (u==0)
            result += 1.0;
        else if (u==1)
            result += 0.0;
        else if (u==2)
            result += 0.4;
        else if (u==3)
            result += 0.8;
        else if (u==4)
            result += 0.4;
        else if (u==5)
            result += 0.0;
        else // u == 6
            result += 1.0;
    }

    return result;
}

double Chromosome::cycTrap(double fHigh, double fLow) const {
    int i, j;
    int u;
    int TRAP_M = length / (TRAP_K-1);
    if (length % (TRAP_K-1) != 0)
        outputErrMsg ("TRAP_k doesn't divide length for Cyclic Setting");
    double result = 0;
    for (i = 0; i < TRAP_M; i++) {
        u = 0;
        int idx = i * TRAP_K - i;
        for (j = 0; j < TRAP_K; j++) {
            int pos = idx + j;
            if (pos == length)
                pos = 0;
            else if (pos > length)
                outputErrMsg ("CYCLIC BUG");
            //
            u += getVal(pos);
        }
        result += trap (u, fHigh, fLow, TRAP_K);
    }
    return result;
}



double Chromosome::mkTrap (double fHigh, double fLow) const {
    int i, j;
    int u;

    int TRAP_M = length / TRAP_K;

    if (length % TRAP_K != 0)
        outputErrMsg ("TRAP_K doesn't divide length");

    double result = 0;

    for (i = 0; i < TRAP_M; i++) {
        u = 0;
        for (j = 0; j < TRAP_K; j++)
            u += getVal(i * TRAP_K + j);

        result += trap (u, fHigh, fLow, TRAP_K);
    }

    return result;
}

double Chromosome::hTrap (double high, double low) const {
    int u, _u;

    int t = HTRAP_K;
    while (length>t) {
        t=t*HTRAP_K;
    }
    if (length!=t)
        outputErrMsg ("length != HTRAP_K^k, ex: 27, 81, 243, ...");

    double result = 0;
    
    t = HTRAP_K;
    while (length>=t) {
        int HTRAP_M = length / t;
        int HTRAP_T = int(t/HTRAP_K);

        for (int _i=0; _i<HTRAP_M; _i++) {
            u = 0;
            for (int _j=0; _j<HTRAP_K; _j++) {
                _u = 0;
                for (int _k=0; _k<HTRAP_T; _k++) {
                    _u += getVal(_i*t+_j*HTRAP_T+_k);
                }
                if (_u == HTRAP_T) u += 1;
                else if (_u > 0) {
                    u += (HTRAP_K+1);
                    break;
                }
            }
            if (length>t) {
                result += (int(t/HTRAP_K)*trap (u, high, high, HTRAP_K));
            }
            else {
                result += (int(t/HTRAP_K)*trap (u, high, low, HTRAP_K));
            }
        }
        t=t*HTRAP_K;
    }

    return result;
}

Chromosome & Chromosome::operator= (const Chromosome & c)
{
    int i;

    if (length != c.length) {
        length = c.length;
        delete[]gene;
        init (length);
    }

    evaluated = c.evaluated;
    fitness = c.fitness;

    for (i = 0; i < length; i++)
        gene[i] = c.gene[i];

    return *this;
}


void Chromosome::printf () const
{
    int i;
    for (i = 0; i < length; i++)
        ::printf ("%d", gene[i]);
}


int Chromosome::getLength () const
{
    return length;
}


double Chromosome::getMaxFitness () const {

    double maxF;
    int t,HTRAP_M;

    switch (FITNESS_FUNCTION) {
        case 0:
            maxF = length;
            break;
        case 1:
            maxF = length/TRAP_K;
            break;
        case 2:
            maxF = length/6;
            break;
        case 3:
            maxF =  length/(TRAP_K - 1);
            break;
        case 4:
            t = HTRAP_K;
            maxF = 0;
            while (length>=t) {
                HTRAP_M = int(length/t);
                maxF += HTRAP_M*int(t/HTRAP_K);
                t = t*HTRAP_K;
            }
            break;
        default:
            // Never converge
            maxF = INF;
    }

    return maxF - EPSILON;

}