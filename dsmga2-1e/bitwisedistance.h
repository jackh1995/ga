
#ifndef _BITWISEDISTANCE_H_
#define _BITWISEDISTANCE_H_

class BitwiseDistance {

public:

    // Count the number of ones in inBinary(x) where x is some nonnegative integer.
    int countOne(unsigned long bitString) {
        return __builtin_popcountl(bitString);
    }
    
    // Get the Hamming distance of two nonnegative integers.
    int getHammingDistance(unsigned long a, unsigned long b) {
        return countOne(a^b);
    }

};


#endif
