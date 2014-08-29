/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Basic functionality for base manipulation.
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_MODEL_NUCLEOTIDES_HH
#define EAGLE_MODEL_NUCLEOTIDES_HH

#include <math.h>
#include <string>
#include <vector>
#include <cassert>


namespace eagle
{
namespace model
{

/*
 * Symbol | Base Represented | Normalized
 *     A  |   A              |    A
 *     C  |       C          |    C
 *     G  |           G      |    G
 *     T  |               T  |    T
 *     U  |               U  |    T
 *     W  |   A           T  |    A
 *     S  |       C   G      |    C
 *     M  |   A   C          |    C
 *     K  |           G   T  |    G
 *     R  |   A       G      |    G
 *     Y  |       C       T  |    C
 *     B  |       C   G   T  |    T
 *     D  |   A       G   T  |    G
 *     H  |   A   C       T  |    C
 *     V  |   A   C   G      |    A
 *     N  |   A   C   G   T  |    N
 *
 * (*) from the International Union of Pure and Applied Chemistry
 */


class IUPAC
{
public:
    IUPAC(bool neg=false) : neg_(neg)
    {
        if (encode.empty())  // uninitialized
        {
            encode.resize(6,std::vector<char>(256,'-'));
            const std::string bases   ("abcdghkmnrstuvwy");  // bases before transformation
            const std::string nbases  ("atcggcgcngcttaac");  // [0] normalized bases
            const std::string cbases  ("tvghcdmknysaabwr");  // [1] complemented bases
            const std::string ncbases ("tagccgcgncgaattg");  // [2] normalized complemented bases
                                                             // [3] binary values, where bit0='A', bit1='C', bit2='G', bit3='T'
#define A 1
#define C 2
#define G 4
#define T 8
            const char binbases[] = {A, C|G|T ,C, A|G|T, G, A|C|T, G|T, A|C, A|C|G|T, A|G, C|G, T, T, A|C|G, A|T, C|T};
#undef A
#undef C
#undef G
#undef T
                                                             // [4] normalized bcl-like values: A=0, C=1, G=2, T=3 + N=4
#define A 0
#define C 1
#define G 2
#define T 3
#define N 4
            const char nbclbases[] = {A,T,C,G,G,C,G,C,N,G,C,T,T,A,A,C};
#undef A
#undef C
#undef G
#undef T
#undef N

            assert( bases.size() == nbases.size() );
            assert( bases.size() == cbases.size() );
            assert( bases.size() == ncbases.size() );
            assert( bases.size() == sizeof(binbases) );
            assert( bases.size() == sizeof(nbclbases) );

            for (unsigned int i=0; i<bases.size(); ++i)
            {
                encode[0][tolower(bases[i])] = tolower(nbases[i]);
                encode[0][toupper(bases[i])] = toupper(nbases[i]);
                encode[1][tolower(bases[i])] = tolower(cbases[i]);
                encode[1][toupper(bases[i])] = toupper(cbases[i]);
                encode[2][tolower(bases[i])] = tolower(ncbases[i]);
                encode[2][toupper(bases[i])] = toupper(ncbases[i]);
                encode[3][tolower(bases[i])] = binbases[i];
                encode[3][toupper(bases[i])] = binbases[i];
                encode[4][tolower(bases[i])] = nbclbases[i];
                encode[4][toupper(bases[i])] = nbclbases[i];
            }

            // Initialise ASCII-complemented locations, so that complemented ASCII values
            // such as ~'A' are understood as the complement of the corresponding base
            for (unsigned int i=0; i<bases.size(); ++i)
            {
                encode[0][static_cast<unsigned char>(~ tolower(cbases[i]))] = tolower(nbases[i]);
                encode[0][static_cast<unsigned char>(~ toupper(cbases[i]))] = toupper(nbases[i]);
                encode[1][static_cast<unsigned char>(~ tolower(cbases[i]))] = tolower(cbases[i]);
                encode[1][static_cast<unsigned char>(~ toupper(cbases[i]))] = toupper(cbases[i]);
                encode[2][static_cast<unsigned char>(~ tolower(cbases[i]))] = tolower(ncbases[i]);
                encode[2][static_cast<unsigned char>(~ toupper(cbases[i]))] = toupper(ncbases[i]);
                encode[3][static_cast<unsigned char>(~ tolower(cbases[i]))] = binbases[i];
                encode[3][static_cast<unsigned char>(~ toupper(cbases[i]))] = binbases[i];
                encode[4][static_cast<unsigned char>(~ tolower(cbases[i]))] = nbclbases[i];
                encode[4][static_cast<unsigned char>(~ toupper(cbases[i]))] = nbclbases[i];
            }

            // Initialise bin->IUPAC
            encode[5][0] = '=';
            for (unsigned int i=0; i<bases.size(); ++i)
            {
                encode[5][binbases[i]] = toupper(bases[i]);
            }
            assert( encode[5][8] == 'U' );
            encode[5][8] = 'T'; // Change the 'U' back to a 'T'

            // Global "initialised" marker
            encode[0][0] = '=';
        }
    }
    char norm(const char c)                     { return encode[0][(unsigned char)c]; }
    char cmpl(const char c)                     { return encode[1][(unsigned char)c]; }
    // normalization => returns either encode[0] or encode[2] depending whether we want the complement or not
    char normalize(const char c,const bool neg) { return encode[static_cast<int>(neg << 1)][(int)c]; }
    char bin(const char c)                      { return encode[3][(unsigned char)c]; }
    char normalizedBcl(const char c)            { return encode[4][(unsigned char)c]; }
    // for completeness, when used as a functor it performs a normalization:
    char operator()(const char c)               { return normalize(c,neg_); }

    char normFromBcl(const char c)              { static const char bases[]="ACGT"; return c?bases[((unsigned char)c)%4]:'N'; }
    char binToIupac(const char c)               { return encode[5][(unsigned char)c]; }

private:
    static std::vector< std::vector<char> > encode;
    bool neg_;
};



} // namespace model
} // namespace eagle

#endif // EAGLE_MODEL_NUCLEOTIDES_HH
