/***************************************************************************
 *   Copyright (C) 2011 by TEIL                                        *
 *                                                                         *
 ***************************************************************************/
#ifndef _CHROMOSOME_H
#define _CHROMOSOME_H

#include <unordered_map>
#include "global.h"
#include "nk-wa.h"

using namespace std;

class Chromosome {

public:

    static enum Function {
        ONEMAX=0,
        MKTRAP=1,
        FTRAP=2,
        CYCTRAP=3,
        NK=4,
        SPINGLASS=5,
        SAT=6
    } function;


    Chromosome ();
    Chromosome (int n_ell);

    ~Chromosome ();

    bool hasSeen() const;

    // Greedy hillclibing
    bool GHC();
    
    // ! ????????????????? 這個有被implement嗎
    void steepestDescent();

    // ! 這三個差在哪裡
    void init (int _ell);
    void init0 (int _ell);
    void initR (int _ell);

    // ! 必須確認一下這在幹嘛
    bool tryFlipping (int index);

    int getVal (int index) const {
        assert (index >= 0 && index < length);

        int q = quotientLong(index);
        int r = remainderLong(index);

        if ( (gene[q] & (1lu << r)) == 0 )
            return 0;
        else
            return 1;
    }

    // Set the value of a ch at index as val
    void setVal (int index, int val) {

        assert (index >= 0 && index < length);

        if (getVal(index) == val) return;

        setValF(index, val);
        key ^= zKey[index];
    }

    // Get the key of a ch
    unsigned long getKey () const {
        return key;
    }


    // ! 這和setVal差在哪裡
    void setValF (int index, int val) {

        assert (index >= 0 && index < length);
        //outputErrMsg ("Index overrange in Chromosome::operator[]");

        int q = quotientLong(index);
        int r = remainderLong(index);

        if (val == 1)
            gene[q] |= (1lu<<r);
        else
            gene[q] &= ~(1lu<<r);

        evaluated = false;
    }

    // ! 這和try flip差在哪裡
    void flip (int index) {
        assert (index >= 0 && index < length);

        int q = quotientLong(index);
        int r = remainderLong(index);

        gene[q] ^= (1lu<<r);
        key ^= zKey[index];

        evaluated = false;
    }

    /** real evaluator */
    
    // Call the fitness function once
    double evaluate ();

    // Verify whether the ch has been evaluated or not
    bool isEvaluated () const;

    // Verify whether two chormosomes are equivalent or not
    bool operator== (const Chromosome & c) const;
    Chromosome & operator= (const Chromosome & c);

    // Get the fitness of a ch
    double getFitness ();
    
    // A collection of benchmark problems
    double trap (int u, double high, double low, int trapK) const;
    double oneMax () const;
    double mkTrap (double high, double low) const;
    double cycTrap(double fHigh, double fLow) const;
    double fTrap () const;
    double spinGlass () const;
    double nkFitness() const;
    double satFitness() const;

    // Get the length of a ch
    int getLength () const;
    
    // Set the length of a ch (this function is undefined)
    void setLength ();

    // Return the maximum possible fitness of each problem
    double getMaxFitness () const;


public:
    
    // ! What is the diff between these four members?!
    static int nfe;
    static int lsnfe;
    static int hitnfe;
    static bool hit;
    
    // ! What the fuck is this?    
    static unordered_map<unsigned long, double> cache;

protected:

    // These are the properties of a ch
    unsigned long *gene;
    int length;
    
    // ! vs length
    int lengthLong;
    double fitness;
    bool evaluated;
    
    // ! What is this
    unsigned long key;

};

#endif
