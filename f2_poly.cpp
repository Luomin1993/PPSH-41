#include <ginac/ginac.h>
#include <vector>
#include <map>
// g++ f2_poly.cpp -o hello -lginac -lcln

using namespace std;

bool ppsh_has(GiNaC::ex const& POLYNOMIAL,GiNaC::ex const& MONINOMIAL)
{
    GiNaC::ex REPR = GiNaC::wild();
    GiNaC::exmap LIST;
    if(GiNaC::has(POLYNOMIAL,REPR*MONINOMIAL) && !GiNaC::match(POLYNOMIAL,MONINOMIAL,LIST) && GiNaC::match(POLYNOMIAL,MONINOMIAL+REPR,LIST)){return true;}
    return false;
}

bool ppsh_has_mul(GiNaC::ex const& POLYNOMIAL,GiNaC::ex const& MONINOMIAL)
{
    GiNaC::ex REPR = GiNaC::wild();
    GiNaC::exset FOUND_SET;
    if(GiNaC::find(POLYNOMIAL,REPR*MONINOMIAL,FOUND_SET)){return true;}
    return false;
}

// ================== PPSH Core Functions ======================
std::vector<int> ppsh_eq2co(GiNaC::ex const& POLYNOMIAL,std::vector<GiNaC::ex> const& ITEMS_3rd)
{
    std::vector<int> COEFF;
    for (int INDEX_ITEM = 0; INDEX_ITEM < ITEMS_3rd.size(); ++INDEX_ITEM)
    {
        bool HAS_ITEM = false;
        for (int INDEX_MONO = 0; INDEX_MONO < POLYNOMIAL.nops(); ++INDEX_MONO)
        {
            if (  GiNaC::is_a<GiNaC::numeric>(POLYNOMIAL.op(INDEX_MONO)/ITEMS_3rd[INDEX_ITEM]) )
            {
                COEFF.push_back( GiNaC::ex_to<GiNaC::numeric>(POLYNOMIAL.op(INDEX_MONO)/ITEMS_3rd[INDEX_ITEM]).to_int() );
                HAS_ITEM = true;
                break;
            }
        }
        if(!HAS_ITEM){COEFF.push_back(0);}
    }
    return COEFF;
}

GiNaC::ex ppsh_co2eq(std::vector<int> const& COEFF,std::vector<GiNaC::ex> const& ITEMS_3rd)
{
    GiNaC::ex POLYNOMIAL;
    if (COEFF.size() != ITEMS_3rd.size()){return POLYNOMIAL;}
    for (int INDEX_ITEM = 0; INDEX_ITEM < ITEMS_3rd.size(); ++INDEX_ITEM){POLYNOMIAL += ITEMS_3rd[INDEX_ITEM] * COEFF[INDEX_ITEM];}
    return POLYNOMIAL;
}

std::vector<GiNaC::ex> ppsh_monomials_of(GiNaC::ex const& POLYNOMIAL)
{
    std::vector<GiNaC::ex> MONINOMIALS;
    for (int INDEX_MONO = 0; INDEX_MONO < POLYNOMIAL.nops(); ++INDEX_MONO){MONINOMIALS.push_back(POLYNOMIAL.op(INDEX_MONO));}
    return MONINOMIALS;
}

GiNaC::ex ppsh_subs_by_dict_mono(GiNaC::ex const& MONINOMIAL,GiNaC::exmap TAB_MONOMIALS_PAIR)
{
    //if (GiNaC::is_a<GiNaC::numeric>(MONINOMIAL)){return MONINOMIAL;}
    cout<<MONINOMIAL<<endl;
    bool OVER_SUBS = true;
    GiNaC::ex REPR = GiNaC::wild();
    GiNaC::ex POLYNOMIAL_SUBS = MONINOMIAL;
    GiNaC::ex POLYNOMIAL_AFT_SUBS = MONINOMIAL;
    for (GiNaC::exmap::iterator MONINOMIALS_PAIR = TAB_MONOMIALS_PAIR.begin();MONINOMIALS_PAIR != TAB_MONOMIALS_PAIR.end();MONINOMIALS_PAIR++)
    {
        POLYNOMIAL_AFT_SUBS = POLYNOMIAL_SUBS;
        if (POLYNOMIAL_AFT_SUBS==MONINOMIALS_PAIR->first){POLYNOMIAL_AFT_SUBS = MONINOMIALS_PAIR->second;}
        //if (POLYNOMIAL_AFT_SUBS!=MONINOMIALS_PAIR->first){cout<<POLYNOMIAL_AFT_SUBS<<endl;cout<<MONINOMIALS_PAIR->first<<endl;}
        // cout<<POLYNOMIAL_AFT_SUBS<<" == "<<MONINOMIALS_PAIR->first<<endl;
        POLYNOMIAL_AFT_SUBS = GiNaC::subs(POLYNOMIAL_AFT_SUBS,MONINOMIALS_PAIR->first+REPR == MONINOMIALS_PAIR->second+REPR).expand();
        POLYNOMIAL_AFT_SUBS = GiNaC::subs(POLYNOMIAL_AFT_SUBS,MONINOMIALS_PAIR->first*REPR == MONINOMIALS_PAIR->second*REPR).expand();
        if(POLYNOMIAL_AFT_SUBS!=POLYNOMIAL_SUBS){OVER_SUBS=false;}
        if(!OVER_SUBS)
        {
            cout<<" --- "<<POLYNOMIAL_AFT_SUBS<<endl;
            for (int INDEX_MONO = 0; INDEX_MONO < POLYNOMIAL_AFT_SUBS.nops(); ++INDEX_MONO)
            {
                POLYNOMIAL_AFT_SUBS = POLYNOMIAL_AFT_SUBS - POLYNOMIAL_AFT_SUBS.op(INDEX_MONO) 
                                      + ppsh_subs_by_dict_mono(POLYNOMIAL_AFT_SUBS.op(INDEX_MONO),TAB_MONOMIALS_PAIR);
            }
        }
        POLYNOMIAL_SUBS = POLYNOMIAL_AFT_SUBS;
    }
    if(OVER_SUBS){cout<<" return "<<POLYNOMIAL_SUBS<<endl;return POLYNOMIAL_SUBS;}
}

void ppsh_subs_by_dict(GiNaC::ex & POLYNOMIAL,GiNaC::exmap TAB_MONOMIALS_PAIR)
{
    for (int INDEX_MONO = 0; INDEX_MONO < POLYNOMIAL.nops(); ++INDEX_MONO)
    {
        POLYNOMIAL = POLYNOMIAL - POLYNOMIAL.op(INDEX_MONO) + ppsh_subs_by_dict_mono(POLYNOMIAL.op(INDEX_MONO),TAB_MONOMIALS_PAIR);
    }
    return;
}

std::vector<int> ppsh_coeff_pow_of(GiNaC::ex const& MONINOMIAL,std::vector<GiNaC::symbol> const& ITEMS_VARS)
{
    std::vector<int> COEFF_POW;
    //if(MONINOMIAL.nops()>2){return COEFF_POW;}
    cout<<MONINOMIAL<<" ";
    for (int INDEX_ITEM = 0; INDEX_ITEM < ITEMS_VARS.size(); ++INDEX_ITEM)
    {
        COEFF_POW.push_back(MONINOMIAL.degree(ITEMS_VARS[INDEX_ITEM]));
        cout<<MONINOMIAL.degree(ITEMS_VARS[INDEX_ITEM])<<" ";
    }
    cout<<endl;
    return COEFF_POW;
}


GiNaC::ex ppsh_simplify_mono(GiNaC::ex MONINOMIAL,std::vector<GiNaC::symbol> const& ITEMS_VARS)
{
    GiNaC::ex MONINOMIAL_SMF = MONINOMIAL;
    std::vector<int> COEFF_POW = ppsh_coeff_pow_of(MONINOMIAL,ITEMS_VARS);
    for (int INDEX_ITEM = 0; INDEX_ITEM < COEFF_POW.size(); ++INDEX_ITEM)
    {
        if(COEFF_POW[INDEX_ITEM]<2){continue;}
        MONINOMIAL = ITEMS_VARS[INDEX_ITEM]*(MONINOMIAL/GiNaC::pow(ITEMS_VARS[INDEX_ITEM],COEFF_POW[INDEX_ITEM]));
        //cout<<" ------ "<<GiNaC::pow(ITEMS_VARS[INDEX_ITEM],COEFF_POW[INDEX_ITEM])<<endl;
    }
    //return MONINOMIAL_SMF;
    cout<<MONINOMIAL_SMF<<" -- smf --> "<<MONINOMIAL<<endl;
    return MONINOMIAL;
}

int ppsh_coeff_of_mono(GiNaC::ex const& MONINOMIAL,std::vector<GiNaC::symbol> const& ITEMS_VARS,std::vector<GiNaC::ex> const& ITEMS_3rd)
{
    int COEFF;
    GiNaC::ex MONINOMIAL_TMP = ppsh_simplify_mono(MONINOMIAL,ITEMS_VARS);
    for (int INDEX_ITEM = 0; INDEX_ITEM < ITEMS_3rd.size(); ++INDEX_ITEM)
    {
        if (  GiNaC::is_a<GiNaC::numeric>(MONINOMIAL_TMP/ITEMS_3rd[INDEX_ITEM]) )
        {
            return GiNaC::ex_to<GiNaC::numeric>(MONINOMIAL_TMP/ITEMS_3rd[INDEX_ITEM]).to_int();
        }
    }        
    return 0;
}

void ppsh_simplify(GiNaC::ex & POLYNOMIAL,std::vector<GiNaC::symbol> const& ITEMS_VARS,std::vector<GiNaC::ex> const& ITEMS_3rd)
{
    GiNaC::ex MONINOMIAL_TMP;
    GiNaC::ex POLYNOMIAL_SMF = POLYNOMIAL;
    cout<<POLYNOMIAL<<endl;
    // descend order;
    for (int INDEX_MONO = 0; INDEX_MONO < POLYNOMIAL.nops(); ++INDEX_MONO)
    {
        MONINOMIAL_TMP = POLYNOMIAL.op(INDEX_MONO);
        POLYNOMIAL_SMF = POLYNOMIAL_SMF - MONINOMIAL_TMP 
                                + ppsh_simplify_mono( MONINOMIAL_TMP,ITEMS_VARS );
        cout<<" f3 = "<<POLYNOMIAL_SMF<<endl;
    }
    // descend coeff;
    std::vector<int> COEFF = ppsh_eq2co(POLYNOMIAL_SMF,ITEMS_3rd);
    for (int INDEX_ITEM = 0; INDEX_ITEM < COEFF.size(); ++INDEX_ITEM)
    {
        COEFF[INDEX_ITEM] = COEFF[INDEX_ITEM]%2;
    }
    POLYNOMIAL = ppsh_co2eq(COEFF,ITEMS_3rd);
    return;
}

int ppsh_subs_by_solv(GiNaC::ex const& POLYNOMIAL,GiNaC::exmap TAB_SOLV_8I2O)
{
    GiNaC::ex POLYNOMIAL_SMF = POLYNOMIAL;
    /*for (GiNaC::exmap::iterator SOLV_PAIR = TAB_SOLV_8I2O.begin();SOLV_PAIR != TAB_SOLV_8I2O.end();SOLV_PAIR++)
    {
        POLYNOMIAL_SMF = POLYNOMIAL_SMF.subs(SOLV_PAIR->first == SOLV_PAIR->second);
    }*/
    POLYNOMIAL_SMF = POLYNOMIAL_SMF.subs(TAB_SOLV_8I2O);
    return ( GiNaC::ex_to<GiNaC::numeric>(POLYNOMIAL_SMF).to_int() )%2;
}


// ==================> PPSH Core Functions END

int main(int argc, char const *argv[])
{
    GiNaC::symbol x1("x1"), x2("x2"), x3("x3"), x4("x4"), x5("x5"), x6("x6"), x7("x7"), x8("x8"), x9("x9");
    GiNaC::ex f1 = x1*x2*x3 + x1*x3 + x1 + x2 + x3;
    GiNaC::ex f2 = 4*GiNaC::pow(x1,3)*x2 + 5*x1*GiNaC::pow(x2,2) + 3*x2- GiNaC::pow(x1+x2,2) + 2*GiNaC::pow(x2+2,2) - 8;
    f2 = f2.expand();
    cout<< f1 <<endl;
    cout<< f2 <<endl;
    GiNaC::ex sss = GiNaC::wild();
    GiNaC::ex ddd = GiNaC::wild();
    GiNaC::exmap LIST;
    GiNaC::exset FOUND_SET;
    
    //cout<< GiNaC::has(f1,sss*x1*x3) <<endl;
    //cout<< GiNaC::has(f1,sss*x1*x2*x3) <<endl;
    //cout<< GiNaC::has(f1,sss*x3) <<endl;
    
    cout<< "=====================" <<endl;
    //cout<< GiNaC::match(f1,sss+ddd*x1*x3,LIST) <<endl;
    //cout<< LIST <<endl;
    //cout<< GiNaC::match(f1,sss*x1*x2*x3,LIST) <<endl;
    //cout<< GiNaC::match(f1,sss*x3,LIST) <<endl;
    //cout<< GiNaC::match(f1,sss*x1*x3,LIST) <<endl;
    //cout<< LIST <<endl;
    //cout<< GiNaC::match(f1,x1,LIST) <<endl;
    //cout<< LIST <<endl;
    //cout<< GiNaC::has(f1,x1*x2) <<endl;
    //cout<< f1.coeff(x1*x3,1) <<endl;
    //cout<< GiNaC::subs(GiNaC::subs(f1,x1==0),x2==1) <<endl;
    
    //cout<< ppsh_has(f1,x1*x3) <<endl;
    //cout<< ppsh_has(f1,x1*x2) <<endl;
    //cout<< ppsh_has(f1,x1) <<endl;
    //cout<< ppsh_has(f1,x2*x3) <<endl;
    //cout<< ppsh_has_mul(f2,x2) <<endl;
    //cout<< ppsh_has_mul(f2,x2*GiNaC::pow(x1,3)) <<endl;
    //cout<< ppsh_has_mul(f2,x2*x1) <<endl;

    // cout<< f2.coeff(x1,2) <<endl;
    // cout<< f2.coeff(x1*x2,1) <<endl;
    // cout<< GiNaC::find(f2,sss*x1*x2,FOUND_SET) <<endl;
    // cout<< FOUND_SET <<endl;

    //for(size_t i=0; i<f2.nops(); i++){cout<< f2.op(i) <<endl;}
    
    //std::vector<GiNaC::ex> MONINOMIALS;MONINOMIALS.resize(f2.nops());
    //std::copy(f2.begin(), f2.end(), MONINOMIALS.begin());
    //for(size_t i=0; i<MONINOMIALS.size(); i++){cout<< MONINOMIALS[i] <<endl;}
    
    // -------------- TEST ppsh_eq2co(...) ---------------
    // std::vector<GiNaC::ex> ITEMS_3rd = {x1*x2*x3, x1*x3, x1*x2, x2*x3, x1, x2, x3};
    // std::vector<int> COEFF = ppsh_eq2co(f1,ITEMS_3rd);
    // for(size_t i=0; i<COEFF.size(); i++){cout<< COEFF[i] <<endl;}
    
    // -------------- TEST ppsh_co2eq(...) ---------------
    // std::vector<GiNaC::ex> ITEMS_3rd = {x1*x2*x3, x1*x3, x1*x2, x2*x3, x1, x2, x3};
    // std::vector<int> COEFF = {1,0,0,0,1,1,1};
    // GiNaC::ex f3 = ppsh_co2eq(COEFF,ITEMS_3rd);
    // cout<< f3 <<endl;

    // -------------- TEST ppsh_monomials_of(...) ---------------
    // std::vector<GiNaC::ex> MONINOMIALS = ppsh_monomials_of(f1);
    // for(size_t i=0; i<MONINOMIALS.size(); i++){cout<< MONINOMIALS[i] <<endl;}
    
    // -------------- TEST ppsh_subs_by_dict(...) ---------------
    //std::map<GiNaC::ex,GiNaC::ex> TAB_MONOMIALS_PAIR{{x1*x2,x1+x3},{x1*x3,x2+x3}};
    //  GiNaC::exmap TAB_MONOMIALS_PAIR;TAB_MONOMIALS_PAIR[x1*x2]=x1+x3;TAB_MONOMIALS_PAIR[x1*x3]=x2+x3;
    //  ppsh_subs_by_dict(f1,TAB_MONOMIALS_PAIR);
    //  cout<< f1 <<endl;
    //cout<< ppsh_subs_by_dict_mono(x1*x3,TAB_MONOMIALS_PAIR) <<endl;
    //cout<< f1.subs(sss*x1*x2==sss*(x1+x3)) <<endl;
    //cout<< f1.subs(x1*x3+sss==sss+x2+x3) <<endl;

    // -------------- TEST ppsh_simplify(...) ---------------
    // std::vector<GiNaC::symbol> ITEMS_VARS = {x1,x2,x3};
    // std::vector<GiNaC::ex> ITEMS_3rd = {x1*x2*x3, x1*x3, x1*x2, x2*x3, x1, x2, x3};
    // GiNaC::ex f3 = 2*x2+2*x3+x1+GiNaC::pow(x2,2)+x2*x3;
    //cout<< GiNaC::pow(x1*x2,2).degree(x1) <<endl;
    //ppsh_simplify(f3,ITEMS_VARS,ITEMS_3rd);
    // cout<< ppsh_simplify_mono(GiNaC::pow(x2,2),ITEMS_VARS) <<endl;
    // cout<< ppsh_simplify_mono(x2*x3,ITEMS_VARS) <<endl;
    // cout<< ppsh_simplify_mono(x1,ITEMS_VARS) <<endl;
    // cout<< ppsh_simplify_mono(2*GiNaC::pow(x3,3),ITEMS_VARS) <<endl;
    // cout<< ppsh_simplify_mono(3*GiNaC::pow(x1,-1),ITEMS_VARS) <<endl;
    //cout<< f3 <<endl;

    // -------------- TEST ppsh_subs_by_solv(...) ---------------
    GiNaC::exmap TAB_SOLV_8I2O{{x1,0}, {x6,0}, {x2,0}, {x7,0}, {x3,1}, {x8,1}, {x4,1}, {x5,0}};
    GiNaC::ex f3 = x2*x3 + x1*x2*x3 + x6*x7 + x8*x4 + x3 + 1;
    GiNaC::ex f4 = x2*x5 + x1*x2*x3*x4 + x6*x7 + x8*x4 + x7 + 1;
    cout<<"f3= "<< ppsh_subs_by_solv(f3,TAB_SOLV_8I2O) <<endl;
    cout<<"f4= "<< ppsh_subs_by_solv(f4,TAB_SOLV_8I2O) <<endl;
    return 0;
}