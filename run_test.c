#include "polynomial.h"  /* Polynomial class interface */
#include <stdio.h>  /* for printf() */

// =========== TEST OOP ===========
// int main() 
// {
//     Symbol x1, x2; /* multiple instances of Symbol */
// 
//     Symbol_make(&x1, "x1", 1);
//     Symbol_make(&x2, "x2", 0);
// 
//     printf("Symbol x1(NAME=%s,VALUE=%d)\n", Symbol_getName(&x1), Symbol_getValue(&x1));
//     printf("Symbol x2(NAME=%s,VALUE=%d)\n", Symbol_getName(&x2), Symbol_getValue(&x2));
// 
//     Symbol_change(&x1, "x1", 0);
//     Symbol_change(&x2, "x2", 1);
// 
//     printf("Symbol x1(NAME=%s,VALUE=%d)\n", Symbol_getName(&x1), Symbol_getValue(&x1));
//     printf("Symbol x2(NAME=%s,VALUE=%d)\n", Symbol_getName(&x2), Symbol_getValue(&x2));
// 
//     return 0;
// }

// =========== TEST MONOMIAL ===========
// int main() 
// {
//     Symbol x1, x2, x3;
// 
//     Symbol_make(&x1, "x1", 1);
//     Symbol_make(&x2, "x2", 0);
//     Symbol_make(&x3, "x3", 1);
// 
//     Monomial fx;
// 
//     Symbol ARR_SYMBOLS[3] = {x1,x2,x3};
//     int16_t ARR_COEFF[3] = {1,0,1};
//     int16_t ARR_POWER[3] = {1,0,2};
//     Monomial_make(&fx, ARR_SYMBOLS, ARR_COEFF,ARR_POWER, 3);
// 
//     printf("Monomial fx =%s, VALUE=%d)\n", Monomial_getMonimial(&fx), Monomial_getMonimialValue(&fx)); 
//     printf("Monomial fx LEN=%d)\n", strlen(Monomial_getMonimial(&fx))); 
// 
//     return 0;
// }


// =========== TEST POLYNOMIAL ===========
int main() 
{
    Symbol x1, x2, x3;

    Symbol_make(&x1, "x1", 1);
    Symbol_make(&x2, "x2", 0);
    Symbol_make(&x3, "x3", 1);

    Monomial m1;

    Symbol ARR_SYMBOLS[3] = {x1,x2,x3};
    int16_t ARR_COEFF_1[3] = {1,0,1};
    int16_t ARR_POWER_1[3] = {1,0,2};
    Monomial_make(&m1, ARR_SYMBOLS, ARR_COEFF_1,ARR_POWER_1, 3);

    Monomial m2;

    Symbol ARR_SYMBOLS[3] = {x1,x2,x3};
    int16_t ARR_COEFF_2[3] = {0,1,1};
    int16_t ARR_POWER_2[3] = {0,1,1};
    Monomial_make(&m2, ARR_SYMBOLS, ARR_COEFF_2,ARR_POWER_2, 3);

    Monomial m3;

    Symbol ARR_SYMBOLS[3] = {x1,x2,x3};
    int16_t ARR_COEFF_3[3] = {1,1,0};
    int16_t ARR_POWER_3[3] = {1,1,0};
    Monomial_make(&m3, ARR_SYMBOLS, ARR_COEFF_3,ARR_POWER_3, 3);

    printf("Monomial fx =%s, VALUE=%d)\n", Monomial_getMonimial(&fx), Monomial_getMonimialValue(&fx)); 
    printf("Monomial fx LEN=%d)\n", strlen(Monomial_getMonimial(&fx))); 

    return 0;
}