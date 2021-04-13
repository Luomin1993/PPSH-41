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
    strcpy(pThisPolynomial->STR_POLYNOMIAL, "");
    for (int INDEX = 0; INDEX < pThisPolynomial->DIM_MONOMIALS; ++INDEX)
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


// ------------------------ OPERATIONS API --------------------------------
void init_monomial_by(Monomial* pMonomial_res,Monomial const * const pMonomial_1st)
{
    pMonomial_res->DIM_SYMBOLS = pMonomial_1st->DIM_SYMBOLS;
    pMonomial_res->ARR_SYMBOLS = pMonomial_1st->ARR_SYMBOLS;
    pMonomial_res->ARR_POWER = (int16_t*)malloc( sizeof(int16_t)*pMonomial_res->DIM_SYMBOLS );
    memcpy(pMonomial_res->ARR_POWER,pMonomial_1st->ARR_POWER,sizeof(int16_t)*pMonomial_1st->DIM_SYMBOLS);
    pMonomial_res->ARR_COEFF = (int16_t*)malloc( sizeof(int16_t)*pMonomial_res->DIM_SYMBOLS );
    memcpy(pMonomial_res->ARR_COEFF,pMonomial_1st->ARR_COEFF,sizeof(int16_t)*pMonomial_1st->DIM_SYMBOLS);
    // int16_t TMP_ARR_POWER[pMonomial_res->DIM_SYMBOLS];
    // pMonomial_res->ARR_POWER = TMP_ARR_POWER;
    // for (int INDEX = 0; INDEX < pMonomial_res->DIM_SYMBOLS; ++INDEX){pMonomial_res->ARR_POWER[INDEX] = pMonomial_1st->ARR_POWER[INDEX];}
    // int16_t TMP_ARR_COEFF[pMonomial_res->DIM_SYMBOLS];
    // pMonomial_res->ARR_COEFF = TMP_ARR_COEFF;
    // for (int INDEX = 0; INDEX < pMonomial_res->DIM_SYMBOLS; ++INDEX){pMonomial_res->ARR_COEFF[INDEX] = pMonomial_1st->ARR_COEFF[INDEX];}
}

void print_arr_num(int16_t const * const ARR_COEFF,int16_t DIM_SYMBOLS)
{
    //int DIM_SYMBOLS=sizeof(ARR_COEFF);
    for (int INDEX = 0; INDEX < DIM_SYMBOLS; ++INDEX){printf("%d,",ARR_COEFF[INDEX]);}
    printf("\n");
}

void Monomial_print_info(Monomial const * const pThisMonomial)
{
    printf("---------- INFO of Monomial ----------\n");
}

void Monomial_optProduct(Monomial const * const pMonomial_1st,Monomial const * const pMonomial_2nd,Monomial* pMonomial_res)
{
    if (pMonomial_1st->ARR_SYMBOLS != pMonomial_2nd->ARR_SYMBOLS){return;}
    init_monomial_by(pMonomial_res,pMonomial_1st);
    print_arr_num(pMonomial_2nd->ARR_POWER,pMonomial_2nd->DIM_SYMBOLS);
    Monomial_make(pMonomial_res,pMonomial_res->ARR_SYMBOLS,pMonomial_res->ARR_COEFF,pMonomial_res->ARR_POWER,pMonomial_res->DIM_SYMBOLS);
    print_arr_num(pMonomial_res->ARR_POWER,pMonomial_res->DIM_SYMBOLS);
    print_arr_num(pMonomial_res->ARR_POWER,pMonomial_res->DIM_SYMBOLS);
    for (int INDEX = 0; INDEX < pMonomial_res->DIM_SYMBOLS; ++INDEX)
    {
        //printf("%d\n",INDEX);
        //print_arr_num(pMonomial_res->ARR_POWER,pMonomial_res->DIM_SYMBOLS);
        pMonomial_res->ARR_POWER[INDEX] = pMonomial_res->ARR_POWER[INDEX]+pMonomial_2nd->ARR_POWER[INDEX];
        //(pMonomial_res->ARR_POWER)[INDEX] = (pMonomial_res->ARR_POWER)[INDEX]+1;
    }
    print_arr_num(pMonomial_res->ARR_POWER,pMonomial_res->DIM_SYMBOLS);
    Monomial_make(pMonomial_res,pMonomial_res->ARR_SYMBOLS,pMonomial_res->ARR_COEFF,pMonomial_res->ARR_POWER,pMonomial_res->DIM_SYMBOLS);
}