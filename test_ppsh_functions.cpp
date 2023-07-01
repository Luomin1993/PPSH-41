#include "ppsh_functions.hpp"
// g++ test_ppsh_functions.cpp -lginac -lcln -L. -lppsh -o exe

// =======================  Testing Functions  =======================
void test__pseudo_divide()
{
    GiNaC::ex POLYNOMIAL = 4*a*a*a*a*a + 2*a*a*a*a + a*a*a + a + 1;
    GiNaC::ex POLYNOMIAL_DIV = a*a*a + a*a + 1;
    cout<<pseudo_divide(POLYNOMIAL, POLYNOMIAL_DIV, 3)<<endl;
}

void test__operator_add()
{
    GiNaC::ex POLYNOMIAL_1 = 2*a*a*a*a*a + 2*a*a*a*a + a*a*a + a + 1;
    GiNaC::ex POLYNOMIAL_2 = 2*a*a*a*a + a*a*a + 2*a*a + 1;
    GiNaC::ex IRRED_POLY = a*a*a + 2*a*a + 1;
    struct extension_field EXT_FIELD;
    EXT_FIELD.PRIME_NUM = 3;EXT_FIELD.EXTEN_NUM = 3;EXT_FIELD.IRRED_POLY = IRRED_POLY;
    extension_field_element EXT_FIELD_ELE(&EXT_FIELD,POLYNOMIAL_1);
    extension_field_element EXT_FIELD_ELE_OPS(&EXT_FIELD,POLYNOMIAL_2);
    extension_field_element EXT_FIELD_ELE_RES = EXT_FIELD_ELE + EXT_FIELD_ELE_OPS;
    cout<<EXT_FIELD_ELE_RES.POLYNOMIAL<<endl;
}

void test__operator_times()
{
    GiNaC::ex POLYNOMIAL_1 = 2*a*a + a + 1;
    GiNaC::ex POLYNOMIAL_2 = a*a + 2;
    GiNaC::ex IRRED_POLY = a*a*a + 2*a*a + 1;
    struct extension_field EXT_FIELD;
    EXT_FIELD.PRIME_NUM = 3;EXT_FIELD.EXTEN_NUM = 3;EXT_FIELD.IRRED_POLY = IRRED_POLY;
    extension_field_element EXT_FIELD_ELE(&EXT_FIELD,POLYNOMIAL_1);
    extension_field_element EXT_FIELD_ELE_OPS(&EXT_FIELD,POLYNOMIAL_2);
    extension_field_element EXT_FIELD_ELE_RES = EXT_FIELD_ELE * EXT_FIELD_ELE_OPS;
    cout<<EXT_FIELD_ELE_RES.POLYNOMIAL<<endl;
}

void test__order_thm()
{
    GiNaC::ex POLYNOMIAL_1 = 2*a*a + a + 1;
    GiNaC::ex IRRED_POLY = a*a*a*a + a*a*a + a*a + 1;
    struct extension_field EXT_FIELD;
    EXT_FIELD.PRIME_NUM = 3;EXT_FIELD.EXTEN_NUM = 4;EXT_FIELD.IRRED_POLY = IRRED_POLY;
    extension_field_element EXT_FIELD_ELE(&EXT_FIELD,POLYNOMIAL_1);
    extension_field_element EXT_FIELD_ELE_RES = EXT_FIELD_ELE;
    for (int INDEX_i = 1; INDEX_i <= 9; ++INDEX_i)
    {
        if (INDEX_i>1){EXT_FIELD_ELE_RES = EXT_FIELD_ELE_RES*EXT_FIELD_ELE;}
        cout<<"g^"<<INDEX_i<<" = "<<EXT_FIELD_ELE_RES.POLYNOMIAL<<endl;
    }
}

void test__polynomial_ring()
{
    GiNaC::ex IRRED_POLY = a*a*a + 2*a*a + 1;
    struct extension_field EXT_FIELD;
    EXT_FIELD.PRIME_NUM = 3;EXT_FIELD.EXTEN_NUM = 3;EXT_FIELD.IRRED_POLY = IRRED_POLY;
    extension_field_class EXT_FIELD_CLASS(&EXT_FIELD,IRRED_POLY);
    //cout<<EXT_FIELD_CLASS.multi_computing(EXT_FIELD_CLASS.ELEMENTS_SET[3],EXT_FIELD_CLASS.ELEMENTS_SET[4]).POLYNOMIAL<<endl;
    //cout<<EXT_FIELD_CLASS.add_computing(EXT_FIELD_CLASS.ELEMENTS_SET[3],EXT_FIELD_CLASS.ELEMENTS_SET[4]).POLYNOMIAL<<endl;
    //cout<<EXT_FIELD_CLASS.sub_computing(EXT_FIELD_CLASS.ELEMENTS_SET[3],EXT_FIELD_CLASS.ELEMENTS_SET[4]).POLYNOMIAL<<endl;
    //cout<<EXT_FIELD_CLASS.div_computing(EXT_FIELD_CLASS.ELEMENTS_SET[3],EXT_FIELD_CLASS.ELEMENTS_SET[4]).POLYNOMIAL<<endl;
    EXT_FIELD_CLASS.print_multi_tab();
}

int main(int argc, char const *argv[])
{
    //test__pseudo_divide();
    //test__operator_add();
    //test__operator_times();
    test__order_thm();
    //test__polynomial_ring();
    return 0;
}