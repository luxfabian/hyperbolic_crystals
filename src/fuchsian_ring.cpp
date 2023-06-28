/*
 * ./src/fuchsian_ring.cpp
 *
 * Author:  Fabian R. Lux
 * Date:    6/28/23
 *
 * Adjoins the square root of 3 to the ring of integers.
 * 
 */
#include "fuchsian_ring.h"
#include <iostream>

//Default constructor.
Field::Field(void)
{
    this->a = 0;
    this->b = 0;
}

//Initialize a value in the field as a_in + b_in * sqrt(3)
Field::Field(int a_in, int b_in)
{
    this->a = a_in;
    this->b = b_in;
}

//Add two members of the field
Field Field::operator+(const Field &other) const
{
    return Field(this->a + other.a, this->b + other.b);
}

//Multiply two members of the field
Field Field::operator*(const Field &other) const
{
    Field result;
    result.a = this->a * other.a + 3 * this->b * other.b;
    result.b = this->a * other.b + this->b * other.a;
    return result;
}

//Take modulus of the field value w.r.t to m
Field Field::operator%(const int &m) const
{   
    //https://stackoverflow.com/a/44197900/7236657
    return Field((m + ((this->a)%m)) % m, (m + ((this->b)%m)) % m);
}

//Copy one field value to another
void Field::operator=(const Field other)
{
    this->a = other.a;
    this->b = other.b;
}

//Check if two field values are the same
bool Field::operator==(const Field &other) const
{
    return this->a == other.a && this->b == other.b;
}

//Print field value in terminal
void Field::print(void)
{
    std::cout << a << " + (" << b << "*sqrt(3))" << std::endl;
}