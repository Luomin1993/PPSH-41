#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

// ------------------- DEF SYMBOL ----------------------------
#ifndef SYMBOL_H
#define SYMBOL_H

// Symbol 的属性
typedef struct {
    char* NAME; 
    int16_t VALUE; 
} Symbol;

// Symbol 的操作函数，接口函数
void Symbol_make(Symbol * const pThisSymbol, char*  const NAME, int16_t VALUE);
void Symbol_change(Symbol * const pThisSymbol, char*  const NAME, int16_t VALUE);
char* Symbol_getName(Symbol const * const pThisSymbol);
int16_t Symbol_getValue(Symbol const * const pThisSymbol);

#endif /* SYMBOL_H */


// ------------------- DEF MONOMIAL ----------------------------
#ifndef MONOMIAL_H
#define MONOMIAL_H

#define MAX_LEN_STR_MONO 130

// Monomial 的属性
typedef struct {
    Symbol* ARR_SYMBOLS; 
    int16_t* ARR_COEFF;
    int16_t* ARR_POWER;
    int16_t DIM_SYMBOLS; 
    char STR_MONOMIAL[MAX_LEN_STR_MONO];
} Monomial;

// Monomial 的操作函数，接口函数
void Monomial_make(Monomial * const pThisMonomial, Symbol*  const ARR_SYMBOLS, int16_t* ARR_COEFF,int16_t* ARR_POWER, int16_t DIM_SYMBOLS);
char* Monomial_getMonomial(Monomial const * const pThisMonomial);
int16_t Monomial_getMonomialValue(Monomial const * const pThisMonomial);

#endif /* MONOMIAL_H */


// ------------------- DEF POLYNOMIAL ----------------------------
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#define MAX_LEN_STR_POLY 1000

// Polynomial 的属性
typedef struct {
    Symbol* POLY_NAME; 
    Monomial* ARR_MONOMIALS; 
    int16_t* ARR_MONOMIALS_COEFF;
    int16_t DIM_MONOMIALS; 
    char STR_POLYNOMIAL[MAX_LEN_STR_POLY];
} Polynomial;

// Polynomial 的操作函数，接口函数
void Polynomial_make(Polynomial * const pThisPolynomial, Symbol*  const POLY_NAME, Monomial* ARR_MONOMIALS,int16_t* ARR_MONOMIALS_COEFF, int16_t DIM_MONOMIALS);
char* Polynomial_getPolynomial(Polynomial const * const pThisPolynomial);
int16_t Polynomial_getPolynomialValue(Polynomial const * const pThisPolynomial);

#endif /* POLYNOMIAL_H */