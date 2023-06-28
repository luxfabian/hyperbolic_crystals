/*
 * ./src/fuchsian.cpp
 *
 * Author:  Fabian R. Lux
 * Date:    05/01/2023
 *
 * Implementation of the Fuchsian class. An instance of the class represents
 * a member of the Fuchsian group of genus 2.
 * 
 * TODO: keep track of generated words
 */
#include "fuchsian_group.h"
#include "fuchsian_ring.h"
#include <vector>

Field generator[8][2][2] = {
    {{Field(2, 2), Field(-3, 0)}, {Field(3, 0), Field(2, -2)}},     // a1
    {{Field(2, 0), Field(-3, -2)}, {Field(3, -2), Field(2, 0)}},    // a2
    {{Field(2, 0), Field(0, -1)}, {Field(0, -1), Field(2, 0)}},     // b1
    {{Field(-2, 2), Field(-6, -3)}, {Field(6, -3), Field(-2, -2)}}, // b2
    {{Field(2, -2), Field(3, 0)}, {Field(-3, 0), Field(2, 2)}},     // a1 inverse
    {{Field(2, 0), Field(3, 2)}, {Field(-3, 2), Field(2, 0)}},      // a2 inverse
    {{Field(2, 0), Field(0, 1)}, {Field(0, 1), Field(2, 0)}},       // b1 inverse
    {{Field(-2, -2), Field(6, 3)}, {Field(-6, 3), Field(-2, 2)}}    // b2 inverse
};

vector<Fuchsian> get_generators(const int &m)
{
    vector<Fuchsian> generators;
    Fuchsian G;

    for (int g = 0; g < 8; g++)
    {
	G.word.clear();
	G.word.push_back(g+1);
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                G.mat[2 * i + j] = generator[g][i][j] % m;
            }
        }
        generators.push_back(G);
    }

    return generators;
}

void Fuchsian::assign(Field F[4])
{
    for (int i = 0; i < 4; i++)
        this->mat[i] = F[i];
}

Fuchsian::Fuchsian(void)
{
    // initialize to e
    this->mat[0] = Field(1, 0);
    this->mat[1] = Field(0, 0);
    this->mat[2] = Field(0, 0);
    this->mat[3] = Field(1, 0);
}

void Fuchsian::right_action(int generator_id)
{

    Field buffer[4];

    for (int i = 0; i < 4; i++)
    {
        buffer[i] = Field(0, 0);
    }

    // Matrix multiplication
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                buffer[i * 2 + j] = buffer[i * 2 + j] + (mat[i * 2 + k] * generator[generator_id][k][j]);
            }
        }
    }

    // Update word
    this->word.push_back(generator_id);

    // Update mat
    this->assign(buffer);
}

void Fuchsian::operator=(const Fuchsian &other)
{
    for (int i = 0; i < 4; i++)
    {
        this->mat[i] = other.mat[i];
    }

   this->word = other.word;
}

Fuchsian Fuchsian::operator%(const int &m) const
{
    Fuchsian result;

    result.word = this->word;

    for (int i = 0; i < 4; i++)
    {
        result.mat[i] = (this->mat[i]) % m;
    }

    return result;
}

bool Fuchsian::operator==(const Fuchsian &other) const
{
    bool identical = true;
    for (int i = 0; i < 4; i++)
    {
        identical &= (this->mat[i] == other.mat[i]);
        if (!identical)
        {
            break;
        }
    }

    return identical;
}

//Group multiplication
Fuchsian Fuchsian::operator*(const Fuchsian &other) const
{
    Fuchsian result;

    // -- initialize to zero
    for (int i = 0; i < 4; i++)
    {
        result.mat[i] = Field(0, 0);
    }

    // -- matrix multiplication
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                result.mat[2 * i + j] = result.mat[2 * i + j] + (this->mat[2 * i + k] * other.mat[2 * k + j]);
            }
        }
    }

    // -- join the words
    result.word = this->word;
    result.word.insert(result.word.end(), other.word.begin(), other.word.end());

    return result;
}

/*
 * Generate a new hash from two known hashes lhs and rhs.
 *
 * The constants which are added represent
 * 
 *      1. 0x517cc1b727220a95 = inverse of pi 
 *      2. 0x9e3779b9 = inverse golden ratio
 * 
 * Their bit representations consist of a random sequence of
 * 0's and 1's which is ideal for scrambling up the hashes.
 * 
 * Together, the bit-shifts and the addition of the constant above
 * break the symmetry of the XOR operator ^ and lead to a reasonable
 * way of combining to hashes.
*/
size_t hash_combine(size_t lhs, size_t rhs)
{
    // https://stackoverflow.com/a/27952689/7236657
    if constexpr (sizeof(size_t) >= 8)
    {
        // 64-bit architecture
        lhs ^= rhs + 0x517cc1b727220a95 + (lhs << 6) + (lhs >> 2);
    }
    else
    {
        // 32-bit architecture
        lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
    }
    return lhs;
}


size_t FuchsianHash::operator()(const Fuchsian &fuchs) const
{
    // https://stackoverflow.com/a/2595226/7236657

    std::hash<int> hasher;
    size_t seed = 0;

    for (int i = 0; i < 4; i++)
    {
        seed = hash_combine(seed, hasher(fuchs.mat[i].a));
        seed = hash_combine(seed, hasher(fuchs.mat[i].b));
    }

    return seed;
} 
