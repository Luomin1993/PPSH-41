#include "ppsh_functions.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <cstring>
#include <fstream>
#include <map>
#include <ctime>
#include <math.h>
// g++ ppsh_extension_field_computing.cpp -o app_exfc -lginac -lcln
// g++ ppsh_extension_field_computing.cpp -o app_exfc -I$HOME/usr/include -L$HOME/usr/lib -lginac -lcln

using namespace std;

// =======================  Symbols Define  =======================
GiNaC::symbol a("a"),x0("x0"), x1("x1"), x2("x2"), x3("x3"), x4("x4"), x5("x5"), x6("x6"), x7("x7"), x8("x8"), x9("x9"), x10("x10"), x11("x11"), x12("x12"), x13("x13"), x14("x14"), x15("x15"), x16("x16"), x17("x17"), x18("x18"), x19("x19"), x20("x20"), x21("x21"), x22("x22"), x23("x23"), x24("x24"), x25("x25"), x26("x26"), x27("x27"), x28("x28"), x29("x29"), x30("x30"), x31("x31"), x32("x32"), x33("x33"), x34("x34"), x35("x35"), x36("x36"), x37("x37"), x38("x38"), x39("x39"), x40("x40"), x41("x41"), x42("x42"), x43("x43"), x44("x44"), x45("x45"), x46("x46"), x47("x47"), x48("x48"), x49("x49"), x50("x50"), x51("x51"), x52("x52"), x53("x53"), x54("x54"), x55("x55"), x56("x56"), x57("x57"), x58("x58"), x59("x59"), x60("x60"), x61("x61"), x62("x62"), x63("x63"), x64("x64"), x65("x65"), x66("x66"), x67("x67"), x68("x68"), x69("x69"), x70("x70"), x71("x71"), x72("x72"), x73("x73"), x74("x74"), x75("x75"), x76("x76"), x77("x77"), x78("x78"), x79("x79"), x80("x80"), x81("x81"), x82("x82"), x83("x83"), x84("x84"), x85("x85"), x86("x86"), x87("x87"), x88("x88"), x89("x89"), x90("x90"), x91("x91"), x92("x92"), x93("x93"), x94("x94"), x95("x95"), x96("x96"), x97("x97"), x98("x98"), x99("x99"), x100("x100");
std::vector<GiNaC::symbol> ITEMS_VARS {x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63,x64,x65,x66,x67,x68,x69,x70,x71,x72,x73,x74,x75,x76,x77,x78,x79,x80,x81,x82,x83,x84,x85,x86,x87,x88,x89,x90,x91,x92,x93,x94,x95,x96,x97,x98,x99,x100};

// =======================  Basic Functions  =======================
void func_decimal_conv( int INDEX_i, int PRIME_NUM, int EXTEN_NUM, std::vector<int> &COEFF_VEC)
{
    int TRANS_NUM = INDEX_i;
    for (int INDEX_j = 0; INDEX_j < EXTEN_NUM; ++INDEX_j)
    {
        COEFF_VEC[INDEX_j] = TRANS_NUM%PRIME_NUM;
        TRANS_NUM = TRANS_NUM/PRIME_NUM;
    }
    return;
}

struct extension_field
{
    int PRIME_NUM;
    int EXTEN_NUM;
    GiNaC::ex IRRED_POLY;
};

extension_field* construct_extension_field(int const& PRIME_NUM,int const& EXTEN_NUM,GiNaC::ex const& IRRED_POLY)
{
    struct extension_field* pt_EXT_FIELD;
    pt_EXT_FIELD->PRIME_NUM = PRIME_NUM;pt_EXT_FIELD->EXTEN_NUM = EXTEN_NUM;
    pt_EXT_FIELD->IRRED_POLY = IRRED_POLY;
    return pt_EXT_FIELD;
}

GiNaC::ex coeff_mod(const GiNaC::ex& POLYNOMIAL,const int& PRIME_NUM)
{
    GiNaC::ex POLYNOMIAL_RES = 0;GiNaC::ex COEFF_LEAD;int COEFF_LEAD_NUM;
    std::vector<GiNaC::ex> MONOMIALS = ppsh_monomials_of(POLYNOMIAL);
    //for (int INDEX_i = 0; INDEX_i < MONOMIALS.size(); ++INDEX_i)
    for (int INDEX_i = 0; INDEX_i <= POLYNOMIAL.degree(a); ++INDEX_i)
    {
        COEFF_LEAD = POLYNOMIAL.coeff(a,INDEX_i);
        COEFF_LEAD_NUM = GiNaC::ex_to<GiNaC::numeric>(COEFF_LEAD).to_int() % PRIME_NUM;
        if (COEFF_LEAD_NUM<0){COEFF_LEAD_NUM = PRIME_NUM + COEFF_LEAD_NUM;}
        POLYNOMIAL_RES += COEFF_LEAD_NUM*GiNaC::pow(a,INDEX_i);
    }
    return POLYNOMIAL_RES;
}

GiNaC::ex pseudo_divide(const GiNaC::ex& POLYNOMIAL, const GiNaC::ex& POLYNOMIAL_DIV, const int& PRIME_NUM)
{
    GiNaC::ex POLYNOMIAL_RES = POLYNOMIAL;
    //for (int i=POLYNOMIAL_RES.ldegree(a); i<=POLYNOMIAL_RES.degree(a); ++i) {cout << "The a^" << i << "-coefficient is "<< POLYNOMIAL_RES.coeff(a,i) << endl;}
    int DIS_DEG;GiNaC::ex COEFF_LEAD;
    while(POLYNOMIAL_RES.degree(a)>=POLYNOMIAL_DIV.degree(a))
    {
        DIS_DEG = POLYNOMIAL_RES.degree(a) - POLYNOMIAL_DIV.degree(a);
        COEFF_LEAD = POLYNOMIAL_RES.coeff(a,POLYNOMIAL_RES.degree(a));
        POLYNOMIAL_RES = coeff_mod(POLYNOMIAL_RES - (COEFF_LEAD*POLYNOMIAL_DIV*GiNaC::pow(a, DIS_DEG)).expand(), PRIME_NUM);
        //cout<<"POLYNOMIAL_RES = "<<POLYNOMIAL_RES<<endl;
    }
    POLYNOMIAL_RES = coeff_mod(POLYNOMIAL_RES, PRIME_NUM);
    return POLYNOMIAL_RES;
}

class extension_field_element
{
public:
    extension_field* pt_EXT_FIELD;
    int ELEMENT_INDEX;
    string SYMBOL;
    GiNaC::ex POLYNOMIAL;
    extension_field_element(extension_field* pt_EXT_FIELD,const GiNaC::ex& POLYNOMIAL);
    //~extension_field_element();
    extension_field_element operator+(const extension_field_element& EXT_FIELD_ELE_OPS) const;
    extension_field_element operator*(const extension_field_element& EXT_FIELD_ELE_OPS) const;
    bool operator<(const extension_field_element& EXT_FIELD_ELE_OPS) const;
};

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

/*extension_field_element extension_field_element::operator*=(const extension_field_element& EXT_FIELD_ELE_OPS) const
{
    extension_field_element EXT_FIELD_ELE_RES(this->pt_EXT_FIELD,POLYNOMIAL_RES);
    return EXT_FIELD_ELE_RES*EXT_FIELD_ELE_OPS;
}*/

struct index_point
{
    int ELEMENT_INDEX_1;
    int ELEMENT_INDEX_2;
    index_point(int ELEMENT_INDEX_1,int ELEMENT_INDEX_2){this->ELEMENT_INDEX_1 = ELEMENT_INDEX_1;this->ELEMENT_INDEX_2 = ELEMENT_INDEX_2;}
};

index_point make_index_point(int ELEMENT_INDEX_1,int ELEMENT_INDEX_2)
{
    index_point NEW_POINT(ELEMENT_INDEX_1,ELEMENT_INDEX_2);return NEW_POINT;
}

std::vector<std::vector<int> > init_all_tab(int ORDER_NUM)
{
    std::vector<std::vector<int> > TAB_NEW;
    TAB_NEW.resize(ORDER_NUM);
    for (int INDEX_i = 0; INDEX_i < ORDER_NUM; ++INDEX_i){TAB_NEW[INDEX_i].resize(ORDER_NUM);}
    return TAB_NEW;
}

string print_align(int INDEX_i)
{
    string SYMBOL = std::to_string(INDEX_i);
    int INDEX_NUM = INDEX_i;int BIT_NUM = 0;
    if (INDEX_NUM==0){INDEX_NUM=1;}
    while(INDEX_NUM!=0){INDEX_NUM = INDEX_NUM/10;BIT_NUM++;}
    for (int INDEX_j = 0; INDEX_j < 3 - BIT_NUM; ++INDEX_j){SYMBOL += " ";}
    return SYMBOL;
}

class extension_field_class
{
public:
    extension_field* pt_EXT_FIELD;
    GiNaC::ex IRRED_POLY;
    int ORDER_NUM;
    std::vector<extension_field_element> ELEMENTS_SET;
    GiNaC::exmap TAB_POLY2INDEX;
    std::map<extension_field_element, int> TAB_ELEMENT2INDEX;
    std::vector<std::vector<int> > TAB_MULTI;
    std::vector<std::vector<int> > TAB_DIV;
    std::vector<std::vector<int> > TAB_ADD;
    std::vector<std::vector<int> > TAB_SUB;
    extension_field_class(extension_field* pt_EXT_FIELD,const GiNaC::ex& IRRED_POLY);
    extension_field_element multi_computing(const extension_field_element& ELEMENT_1,const extension_field_element& ELEMENT_2);
    extension_field_element add_computing(const extension_field_element& ELEMENT_1,const extension_field_element& ELEMENT_2);
    extension_field_element div_computing(const extension_field_element& ELEMENT_1,const extension_field_element& ELEMENT_2);
    extension_field_element sub_computing(const extension_field_element& ELEMENT_1,const extension_field_element& ELEMENT_2);
    void print_multi_tab();
};

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
        //cout<<POLYNOMIAL<<endl;
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
    /*for (int INDEX_i = 0; INDEX_i < this->ORDER_NUM-1; ++INDEX_i)
    {
        for (int INDEX_j = 0; INDEX_j <= INDEX_i; ++INDEX_j)
        {
            cout<<"g_"<<INDEX_i<<" * g_"<<INDEX_j<<" = g_"<< this->TAB_MULTI[INDEX_i][INDEX_j] <<endl;
        }
    }*/
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
    /*for (int INDEX_i = 0; INDEX_i < this->ORDER_NUM-1; ++INDEX_i)
    {
        for (int INDEX_j = 0; INDEX_j <= INDEX_i; ++INDEX_j)
        {
            cout<<"g_"<<INDEX_i<<" + g_"<<INDEX_j<<" = g_"<< this->TAB_ADD[INDEX_i][INDEX_j] <<endl;
            cout<<"f1 "<<ELEMENTS_SET[INDEX_i].POLYNOMIAL<<" + f2 "<<ELEMENTS_SET[INDEX_j].POLYNOMIAL<<" = f "<< this->ELEMENTS_SET[this->TAB_ADD[INDEX_i][INDEX_j]].POLYNOMIAL <<endl;
            cout<<"g_"<<INDEX_i<<" - g_"<<INDEX_j<<" = g_"<< this->TAB_SUB[INDEX_i][INDEX_j] <<endl;
        }
    }*/
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

// ================================ Polynomial Ring on Extension Field ======================================
/*class extension_field_polynomial
{
public:
    std::vector<GiNaC::ex> MONOMIALS;
    std::vector<extension_field_element> COEFF_VEC;
    extension_field_polynomial(extension_field_class* pt_EXT_FIELD_CLASS,std::vector<GiNaC::symbol> EXT_SYMBOLS);
    print();
};

extension_field_polynomial::extension_field_polynomial(extension_field_class* pt_EXT_FIELD_CLASS,std::vector<GiNaC::symbol> EXT_SYMBOLS)
{
    int ORDER_NUM = pow(pt_EXT_FIELD_CLASS->pt_EXT_FIELD->PRIME_NUM,pt_EXT_FIELD_CLASS->pt_EXT_FIELD->EXTEN_NUM);
    for (int INDEX_i = 0; INDEX_i < EXT_SYMBOLS.size(); ++INDEX_i)
    {
        // ... code ...
    }
}

void extension_field_polynomial::print()
{
    ;
}

class extension_field_polynomial_ring
{
public:
    extension_field_polynomial_ring(extension_field_class* pt_EXT_FIELD_CLASS,std::vector<GiNaC::symbol> EXT_SYMBOLS);
};

extension_field_polynomial_ring::extension_field_polynomial_ring(extension_field_class* pt_EXT_FIELD_CLASS,std::vector<GiNaC::symbol> EXT_SYMBOLS)
{
    pass;
}*/

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
    //test__order_thm();
    test__polynomial_ring();
    return 0;
}


