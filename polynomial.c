#include "polynomial.h"

// ------------------- BASIC FUNCTIONS -----------------------
char* int2char(int16_t NUMBER)
{
    char STRING[15];
    //itoa(NUMBER,STRING,10); 
    sprintf(STRING,"%d", NUMBER);
    char* pSTRING = STRING;
    return pSTRING;
}

void clean_str(char* STR_MONOMIAL)
{
    char STRING[1]={0};
    STR_MONOMIAL = STRING;
}

// ------------------- API SYMBOL ----------------------------
// 构造函数
void Symbol_make(Symbol * const pThisSymbol, char*  const NAME, int16_t VALUE)
{
    pThisSymbol->NAME = NAME;
    pThisSymbol->VALUE = VALUE;
}

void Symbol_change(Symbol * const pThisSymbol, char*  const NAME, int16_t VALUE)
{
    pThisSymbol->NAME = NAME;
    pThisSymbol->VALUE = VALUE;
}

// 获取属性值函数
char* Symbol_getName(Symbol const * const pThisSymbol)
{
    return pThisSymbol->NAME;
}
int16_t Symbol_getValue(Symbol const * const pThisSymbol)
{
    return pThisSymbol->VALUE;
}


// ------------------- API MONOMIAL ----------------------------
void Monomial_make(Monomial * const pThisMonomial, Symbol*  const ARR_SYMBOLS, int16_t* ARR_COEFF,int16_t* ARR_POWER, int16_t DIM_SYMBOLS)
{
    pThisMonomial->ARR_SYMBOLS = ARR_SYMBOLS;
    pThisMonomial->ARR_COEFF = ARR_COEFF;
    pThisMonomial->ARR_POWER = ARR_POWER;
    pThisMonomial->DIM_SYMBOLS = DIM_SYMBOLS;
    //printf("Monomial fx=%s)\n",pThisMonomial->STR_MONOMIAL); 
    //clean_str(pThisMonomial->STR_MONOMIAL);
    strcpy(pThisMonomial->STR_MONOMIAL, "");
    //printf("Monomial fx=%s)\n",pThisMonomial->STR_MONOMIAL); 
    for (int INDEX = 0; INDEX < pThisMonomial->DIM_SYMBOLS; ++INDEX)
    {
        if (pThisMonomial->ARR_COEFF[INDEX]==0){continue;}
        if (strlen(pThisMonomial->STR_MONOMIAL)>0){strcat(pThisMonomial->STR_MONOMIAL,"*");}
        strcat(pThisMonomial->STR_MONOMIAL,Symbol_getName(&ARR_SYMBOLS[INDEX]));
        if (ARR_POWER[INDEX]>1)
        {
            strcat(pThisMonomial->STR_MONOMIAL,"^");
            strcat(pThisMonomial->STR_MONOMIAL,int2char(ARR_POWER[INDEX]));
        }
    }
}

char* Monomial_getMonomial(Monomial const * const pThisMonomial)
{
    return pThisMonomial->STR_MONOMIAL;
}

int16_t Monomial_getMonomialValue(Monomial const * const pThisMonomial)
{
    int16_t VALUE = 1;
	for (int INDEX = 0; INDEX < pThisMonomial->DIM_SYMBOLS; ++INDEX)
    {
        if (pThisMonomial->ARR_COEFF[INDEX]==0){continue;}
        VALUE *= Symbol_getValue( &(pThisMonomial->ARR_SYMBOLS[INDEX]) );
    }
    return VALUE;
}


// ------------------- API POLYNOMIAL ----------------------------
void Polynomial_make(Polynomial * const pThisPolynomial, Symbol*  const POLY_NAME, Monomial* ARR_MONOMIALS,int16_t* ARR_MONOMIALS_COEFF, int16_t DIM_MONOMIALS)
{
    pThisPolynomial->POLY_NAME = POLY_NAME;
    pThisPolynomial->ARR_MONOMIALS = ARR_MONOMIALS;
    pThisPolynomial->ARR_MONOMIALS_COEFF = ARR_MONOMIALS_COEFF;
    pThisPolynomial->DIM_MONOMIALS = DIM_MONOMIALS;
    strcpy(pThisMonomial->STR_POLYNOMIAL, "");
    for (int INDEX = 0; INDEX < pThisMonomial->DIM_MONOMIALS; ++INDEX)
    {
        if(pThisPolynomial->ARR_MONOMIALS_COEFF[INDEX]==0){continue;}
        if (strlen(pThisPolynomial->STR_POLYNOMIAL)>0){strcat(pThisPolynomial->STR_POLYNOMIAL," + ");}
        strcat(pThisPolynomial->STR_POLYNOMIAL,Monomial_getMonomial(&(pThisPolynomial->ARR_MONOMIALS[INDEX])));
    }
}

char* Polynomial_getPolynomial(Polynomial const * const pThisPolynomial)
{
    return pThisPolynomial->STR_POLYNOMIAL;
}

int16_t Polynomial_getPolynomialValue(Polynomial const * const pThisPolynomial)
{
    int16_t VALUE = 0;
	for (int INDEX = 0; INDEX < pThisPolynomial->DIM_MONOMIALS; ++INDEX)
    {
        if (pThisPolynomial->ARR_MONOMIALS_COEFF[INDEX]==0){continue;}
        VALUE += Monomial_getMonomialValue( &(pThisPolynomial->ARR_MONOMIALS[INDEX]) );
    }
    return VALUE%2;
}