#include "ppsh_functions.hpp"
#include <string>
#include <cstring>
#include <ctime>

// g++ ppsh_functions.cpp -lginac -lcln -fPIC -shared -o libppsh.so

extension_field_element::extension_field_element(extension_field* pt_EXT_FIELD,const GiNaC::ex& POLYNOMIAL)
{
    this->pt_EXT_FIELD = pt_EXT_FIELD;
    this->POLYNOMIAL = POLYNOMIAL;
}

bool extension_field_element::operator<(const extension_field_element& EXT_FIELD_ELE_OPS) const
{
    return this->ELEMENT_INDEX < EXT_FIELD_ELE_OPS.ELEMENT_INDEX;
}

extension_field_element extension_field_element::operator+(const extension_field_element& EXT_FIELD_ELE_OPS) const
{
    GiNaC::ex POLYNOMIAL_RES=0;
    extension_field_element EXT_FIELD_ELE_RES(this->pt_EXT_FIELD,POLYNOMIAL_RES);
    if(this->pt_EXT_FIELD != EXT_FIELD_ELE_OPS.pt_EXT_FIELD){return EXT_FIELD_ELE_RES;}
    EXT_FIELD_ELE_RES.POLYNOMIAL = pseudo_divide(this->POLYNOMIAL+EXT_FIELD_ELE_OPS.POLYNOMIAL,this->pt_EXT_FIELD->IRRED_POLY,this->pt_EXT_FIELD->PRIME_NUM);
    return EXT_FIELD_ELE_RES;
}

extension_field_element extension_field_element::operator*(const extension_field_element& EXT_FIELD_ELE_OPS) const
{
    GiNaC::ex POLYNOMIAL_RES=0;
    extension_field_element EXT_FIELD_ELE_RES(this->pt_EXT_FIELD,POLYNOMIAL_RES);
    if(this->pt_EXT_FIELD != EXT_FIELD_ELE_OPS.pt_EXT_FIELD){return EXT_FIELD_ELE_RES;}
    EXT_FIELD_ELE_RES.POLYNOMIAL = pseudo_divide((this->POLYNOMIAL*EXT_FIELD_ELE_OPS.POLYNOMIAL).expand(),this->pt_EXT_FIELD->IRRED_POLY,this->pt_EXT_FIELD->PRIME_NUM);
    return EXT_FIELD_ELE_RES;
}

extension_field_class::extension_field_class(extension_field* pt_EXT_FIELD,const GiNaC::ex& IRRED_POLY)
{
    this->pt_EXT_FIELD = pt_EXT_FIELD;
    this->IRRED_POLY   =   IRRED_POLY;
    this->ORDER_NUM    = pow(this->pt_EXT_FIELD->PRIME_NUM,this->pt_EXT_FIELD->EXTEN_NUM);
    // --------- Constructing ELEMENTS_SET; ---------
    std::vector<int> COEFF_VEC;std::vector<GiNaC::ex> MONOMIALS;
    for (int INDEX_j = 0; INDEX_j < this->pt_EXT_FIELD->EXTEN_NUM; ++INDEX_j){COEFF_VEC.push_back(0);}
    for (int INDEX_j = 0; INDEX_j < this->pt_EXT_FIELD->EXTEN_NUM; ++INDEX_j)
    {
        GiNaC::ex MONOMIAL = 1;
        if (INDEX_j==0){MONOMIALS.push_back(MONOMIAL);continue;}
        for (int INDEX_k = 1; INDEX_k <= INDEX_j; ++INDEX_k){MONOMIAL*=a;}
        MONOMIALS.push_back(MONOMIAL);
    }
    for (int INDEX_i = 1; INDEX_i < this->ORDER_NUM; ++INDEX_i)
    {
        GiNaC::ex POLYNOMIAL = 0;
        func_decimal_conv(INDEX_i,this->pt_EXT_FIELD->PRIME_NUM,this->pt_EXT_FIELD->EXTEN_NUM,COEFF_VEC);
        for (int INDEX_j = 0; INDEX_j < this->pt_EXT_FIELD->EXTEN_NUM; ++INDEX_j)
        {
            POLYNOMIAL += COEFF_VEC[INDEX_j]*MONOMIALS[INDEX_j];
        }
        extension_field_element EXT_FIELD_ELE(this->pt_EXT_FIELD,POLYNOMIAL);
        EXT_FIELD_ELE.ELEMENT_INDEX = INDEX_i-1;
        this->ELEMENTS_SET.push_back(EXT_FIELD_ELE);
        this->TAB_POLY2INDEX[POLYNOMIAL] = EXT_FIELD_ELE.ELEMENT_INDEX;
        this->TAB_ELEMENT2INDEX[EXT_FIELD_ELE] = EXT_FIELD_ELE.ELEMENT_INDEX;
    }
    // --------- Constructing TAB_MULTI and TAB_DIV; ---------
    std::vector<std::vector<int> > TAB_NEW = init_all_tab(this->ORDER_NUM-1);
    this->TAB_MULTI = TAB_NEW;this->TAB_ADD = TAB_NEW;this->TAB_DIV = TAB_NEW;this->TAB_SUB = TAB_NEW;
    for (int INDEX_i = 0; INDEX_i < this->ORDER_NUM-1; ++INDEX_i)
    {
        for (int INDEX_j = 0; INDEX_j <= INDEX_i; ++INDEX_j)
        {
            GiNaC::ex POLYNOMIAL = (this->ELEMENTS_SET[INDEX_i]*this->ELEMENTS_SET[INDEX_j]).POLYNOMIAL;
            this->TAB_MULTI[INDEX_i][INDEX_j] = GiNaC::ex_to<GiNaC::numeric>(this->TAB_POLY2INDEX[POLYNOMIAL]).to_int();
            this->TAB_MULTI[INDEX_j][INDEX_i] = GiNaC::ex_to<GiNaC::numeric>(this->TAB_POLY2INDEX[POLYNOMIAL]).to_int();
            this->TAB_DIV[ TAB_MULTI[INDEX_j][INDEX_i] ][INDEX_i] = INDEX_j;
            this->TAB_DIV[ TAB_MULTI[INDEX_j][INDEX_i] ][INDEX_j] = INDEX_i;
        }
    }
    // --------- Constructing TAB_ADD; ---------
    for (int INDEX_i = 0; INDEX_i < this->ORDER_NUM-1; ++INDEX_i)
    {
        for (int INDEX_j = 0; INDEX_j <= INDEX_i; ++INDEX_j)
        {
            GiNaC::ex POLYNOMIAL = (this->ELEMENTS_SET[INDEX_i] + this->ELEMENTS_SET[INDEX_j]).POLYNOMIAL;
            this->TAB_ADD[INDEX_i][INDEX_j] = GiNaC::ex_to<GiNaC::numeric>(this->TAB_POLY2INDEX[POLYNOMIAL]).to_int();
            this->TAB_ADD[INDEX_j][INDEX_i] = GiNaC::ex_to<GiNaC::numeric>(this->TAB_POLY2INDEX[POLYNOMIAL]).to_int();
            this->TAB_SUB[ TAB_ADD[INDEX_j][INDEX_i] ][INDEX_i] = INDEX_j;
            this->TAB_SUB[ TAB_ADD[INDEX_j][INDEX_i] ][INDEX_j] = INDEX_i;
        }
    }
    return;
}

extension_field_element extension_field_class::multi_computing(const extension_field_element& ELEMENT_1,const extension_field_element& ELEMENT_2)
{
    return this->ELEMENTS_SET[ this->TAB_MULTI[ this->TAB_ELEMENT2INDEX[ELEMENT_1] ][ this->TAB_ELEMENT2INDEX[ELEMENT_2] ] ];
}

extension_field_element extension_field_class::add_computing(const extension_field_element& ELEMENT_1,const extension_field_element& ELEMENT_2)
{
    return this->ELEMENTS_SET[ this->TAB_ADD[ this->TAB_ELEMENT2INDEX[ELEMENT_1] ][ this->TAB_ELEMENT2INDEX[ELEMENT_2] ] ];
}

extension_field_element extension_field_class::div_computing(const extension_field_element& ELEMENT_1,const extension_field_element& ELEMENT_2)
{
    return this->ELEMENTS_SET[ this->TAB_DIV[ this->TAB_ELEMENT2INDEX[ELEMENT_1] ][ this->TAB_ELEMENT2INDEX[ELEMENT_2] ] ];
}

extension_field_element extension_field_class::sub_computing(const extension_field_element& ELEMENT_1,const extension_field_element& ELEMENT_2)
{
    if (this->TAB_ELEMENT2INDEX[ELEMENT_1] == this->TAB_ELEMENT2INDEX[ELEMENT_2])
    {
        GiNaC::ex POLYNOMIAL_ZERO = 0;
        extension_field_element EXT_FIELD_ELE(this->pt_EXT_FIELD,POLYNOMIAL_ZERO);
        return EXT_FIELD_ELE;
    }
    return this->ELEMENTS_SET[ this->TAB_SUB[ this->TAB_ELEMENT2INDEX[ELEMENT_1] ][ this->TAB_ELEMENT2INDEX[ELEMENT_2] ] ];
}

void extension_field_class::print_multi_tab()
{
    cout<<" ====================== Multiplication Table of GF("<<this->pt_EXT_FIELD->PRIME_NUM<<"^"<<this->pt_EXT_FIELD->EXTEN_NUM<<") ======================"<<endl;
    cout<<"*    | ";
    for (int INDEX_i = 0; INDEX_i < this->ORDER_NUM-1; ++INDEX_i){cout<<"g_"<<print_align(INDEX_i)<<"|";}
    cout<<endl;
    for (int INDEX_i = 0; INDEX_i < this->ORDER_NUM-1; ++INDEX_i)
    {
        cout<<"g_"<<print_align(INDEX_i)<<"| ";
        for (int INDEX_j = 0; INDEX_j < this->ORDER_NUM-1; ++INDEX_j)
        {
            cout<<"g_"<< print_align(this->TAB_MULTI[INDEX_i][INDEX_j]) <<" ";
        }
        cout<<endl;
    }
    return;
}