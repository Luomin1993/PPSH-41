#include <ginac/ginac.h>
#include <vector>
#include <map>
#include <math.h>
#include <time.h>
#include<fstream> 
#include<iostream>
// g++ f2_poly.cpp -o hello -lginac -lcln

#define MAX_LEN_LINE 3000

using namespace std;

void print_coeff(std::vector<int> const& COEFF_POW)
{
    cout<<"[ ";
    for (int INDEX_ITEM = 0; INDEX_ITEM < COEFF_POW.size(); ++INDEX_ITEM)
    {
        cout<<COEFF_POW[INDEX_ITEM]<<", ";
    }
    cout<<" ]"<<endl;
}


int ppsh_size_of(GiNaC::ex const& POLYNOMIAL)
{
    if(GiNaC::is_a<GiNaC::add>(POLYNOMIAL))
    {
        return POLYNOMIAL.nops();
    }
    return 1;
}

GiNaC::ex ppsh_index_monomial_of(GiNaC::ex const& POLYNOMIAL,int const& INDEX_MONO)
{
    GiNaC::ex MONINOMIAL;
    if(GiNaC::is_a<GiNaC::add>(POLYNOMIAL))
    {
        MONINOMIAL = POLYNOMIAL.op(INDEX_MONO);
        return MONINOMIAL;
    }
    MONINOMIAL = POLYNOMIAL;return MONINOMIAL;
}

void print_items(GiNaC::ex const& POLYNOMIAL)
{
    cout<<"[ ";
    for (int INDEX_ITEM = 0; INDEX_ITEM < ppsh_size_of(POLYNOMIAL); ++INDEX_ITEM)
    {
        cout<<ppsh_index_monomial_of(POLYNOMIAL,INDEX_ITEM)<<", ";
    }
    cout<<" ]"<<endl;
}

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
    // If it's a monomial;
    for (int INDEX_ITEM = 0; INDEX_ITEM < ITEMS_3rd.size(); ++INDEX_ITEM)
    {
        bool HAS_ITEM = false;
        for (int INDEX_MONO = 0; INDEX_MONO < ppsh_size_of(POLYNOMIAL); ++INDEX_MONO)
        {
            if (  GiNaC::is_a<GiNaC::numeric>(ppsh_index_monomial_of(POLYNOMIAL,INDEX_MONO)/ITEMS_3rd[INDEX_ITEM]) )
            {
                COEFF.push_back( GiNaC::ex_to<GiNaC::numeric>(ppsh_index_monomial_of(POLYNOMIAL,INDEX_MONO)/ITEMS_3rd[INDEX_ITEM]).to_int() );
                HAS_ITEM = true;
                //cout<< ppsh_index_monomial_of(POLYNOMIAL,INDEX_MONO) <<" should be "<<ITEMS_3rd[INDEX_ITEM] <<endl;
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
    if(POLYNOMIAL.nops()==0 && (POLYNOMIAL+11331).nops() >1){MONINOMIALS.push_back((POLYNOMIAL+11331).op(0));}
    return MONINOMIALS;
}

std::vector<int> ppsh_coeff_pow_of(GiNaC::ex const& MONINOMIAL,std::vector<GiNaC::symbol> const& ITEMS_VARS)
{
    std::vector<int> COEFF_POW;
    //if(MONINOMIAL.nops()>1){return COEFF_POW;}
    for (int INDEX_ITEM = 0; INDEX_ITEM < ITEMS_VARS.size(); ++INDEX_ITEM)
    {
        COEFF_POW.push_back(MONINOMIAL.degree(ITEMS_VARS[INDEX_ITEM]));
    }
    return COEFF_POW;
}

int ppsh_coeff_of_mono(GiNaC::ex const& MONINOMIAL)
{
    if (  GiNaC::is_a<GiNaC::numeric>(MONINOMIAL) ){return GiNaC::ex_to<GiNaC::numeric>(MONINOMIAL).to_int();}
    if(GiNaC::is_a<GiNaC::symbol>( MONINOMIAL )){return 1;}
    if(GiNaC::is_a<GiNaC::symbol>( MONINOMIAL.op( MONINOMIAL.nops()-1 ) )){return 1;}
    return GiNaC::ex_to<GiNaC::numeric>( MONINOMIAL.op( MONINOMIAL.nops()-1 ) ).to_int();
}

GiNaC::ex ppsh_simplify_mono(GiNaC::ex MONINOMIAL,std::vector<GiNaC::symbol> const& ITEMS_VARS)
{
    GiNaC::ex MONINOMIAL_SMF = MONINOMIAL;
    std::vector<int> COEFF_POW = ppsh_coeff_pow_of(MONINOMIAL,ITEMS_VARS);
    for (int INDEX_ITEM = 0; INDEX_ITEM < COEFF_POW.size(); ++INDEX_ITEM)
    {
        if(COEFF_POW[INDEX_ITEM]<2){continue;}
        MONINOMIAL = ITEMS_VARS[INDEX_ITEM]*(MONINOMIAL/GiNaC::pow(ITEMS_VARS[INDEX_ITEM],COEFF_POW[INDEX_ITEM]));
    }
    //return MONINOMIAL_SMF;
    //cout<<MONINOMIAL_SMF<<" -- smf --> "<<MONINOMIAL<<endl;
    return MONINOMIAL;
}

void ppsh_simplify(GiNaC::ex & POLYNOMIAL,std::vector<GiNaC::symbol> const& ITEMS_VARS)
{
    if(POLYNOMIAL==0 or POLYNOMIAL==1){return;}
    GiNaC::ex MONINOMIAL_TMP; int COEFF_TMP;
    GiNaC::ex POLYNOMIAL_SMF = POLYNOMIAL;
    //cout<<POLYNOMIAL<<endl;
    // descend order;
    for (int INDEX_MONO = 0; INDEX_MONO < ppsh_size_of(POLYNOMIAL); ++INDEX_MONO)
    {
        MONINOMIAL_TMP = ppsh_index_monomial_of(POLYNOMIAL,INDEX_MONO);
        POLYNOMIAL_SMF = POLYNOMIAL_SMF - MONINOMIAL_TMP 
                                + ppsh_simplify_mono( MONINOMIAL_TMP,ITEMS_VARS );
        //cout<<" f3 = "<<POLYNOMIAL_SMF<<endl;
    }
    // descend coeff;
    // std::vector<int> COEFF = ppsh_eq2co(POLYNOMIAL_SMF,ITEMS_3rd);
    // for (int INDEX_ITEM = 0; INDEX_ITEM < COEFF.size(); ++INDEX_ITEM)
    // {
    //     COEFF[INDEX_ITEM] = COEFF[INDEX_ITEM]%2;
    // }
    // POLYNOMIAL = ppsh_co2eq(COEFF,ITEMS_3rd);
    POLYNOMIAL = POLYNOMIAL_SMF;
    for (int INDEX_MONO = 0; INDEX_MONO < ppsh_size_of(POLYNOMIAL_SMF); ++INDEX_MONO)
    {
        //cout<<" f3 = "<<POLYNOMIAL<<endl;
        MONINOMIAL_TMP = ppsh_index_monomial_of(POLYNOMIAL_SMF,INDEX_MONO);
        COEFF_TMP = ppsh_coeff_of_mono(MONINOMIAL_TMP);
        if(COEFF_TMP==1){continue;}
        POLYNOMIAL = POLYNOMIAL - MONINOMIAL_TMP + (MONINOMIAL_TMP/COEFF_TMP)*(COEFF_TMP%2);
    }
    return;
}

GiNaC::ex ppsh_subs_by_solv(GiNaC::ex const& POLYNOMIAL, GiNaC::exmap TAB_SOLV_8I2O, std::vector<GiNaC::symbol> const& ITEMS_VARS)
{
    GiNaC::ex POLYNOMIAL_SMF = POLYNOMIAL;
    /*for (GiNaC::exmap::iterator SOLV_PAIR = TAB_SOLV_8I2O.begin();SOLV_PAIR != TAB_SOLV_8I2O.end();SOLV_PAIR++)
    {
        POLYNOMIAL_SMF = POLYNOMIAL_SMF.subs(SOLV_PAIR->first == SOLV_PAIR->second);
    }*/
    POLYNOMIAL_SMF = POLYNOMIAL_SMF.subs(TAB_SOLV_8I2O);
    if(POLYNOMIAL_SMF==0){return 0;}
    ppsh_simplify(POLYNOMIAL_SMF,ITEMS_VARS);
    return POLYNOMIAL_SMF;
    //return ( GiNaC::ex_to<GiNaC::numeric>(POLYNOMIAL_SMF).to_int() )%2;
}

GiNaC::ex ppsh_subs_by_dict_mono(GiNaC::ex const& MONINOMIAL,GiNaC::exmap TAB_MONOMIALS_PAIR,std::vector<GiNaC::symbol> const& ITEMS_VARS)
{
    if (GiNaC::is_a<GiNaC::numeric>(MONINOMIAL)){return (GiNaC::ex_to<GiNaC::numeric>( MONINOMIAL ).to_int())%2;}
    if (1==TAB_MONOMIALS_PAIR.count(MONINOMIAL)){return TAB_MONOMIALS_PAIR[MONINOMIAL];}
    //cout<<MONINOMIAL<<endl;
    bool OVER_SUBS = true;
    GiNaC::ex REPR = GiNaC::wild();
    GiNaC::ex POLYNOMIAL_SUBS = MONINOMIAL;
    GiNaC::ex POLYNOMIAL_AFT_SUBS = MONINOMIAL;
    GiNaC::ex POLYNOMIAL_TMP = MONINOMIAL;
    for (GiNaC::exmap::iterator MONINOMIALS_PAIR = TAB_MONOMIALS_PAIR.begin();MONINOMIALS_PAIR != TAB_MONOMIALS_PAIR.end();MONINOMIALS_PAIR++)
    {
        POLYNOMIAL_AFT_SUBS = POLYNOMIAL_SUBS;
        //if (POLYNOMIAL_AFT_SUBS==MONINOMIALS_PAIR->first){POLYNOMIAL_AFT_SUBS = MONINOMIALS_PAIR->second;}
        //if (POLYNOMIAL_AFT_SUBS!=MONINOMIALS_PAIR->first){cout<<POLYNOMIAL_AFT_SUBS<<endl;cout<<MONINOMIALS_PAIR->first<<endl;}
        //cout<<POLYNOMIAL_AFT_SUBS<<" ==> "<<MONINOMIALS_PAIR->first<<endl;
        ppsh_simplify(POLYNOMIAL_AFT_SUBS,ITEMS_VARS);
        //cout<<POLYNOMIAL_AFT_SUBS<<" ==> "<<MONINOMIALS_PAIR->first<<endl;
        if (GiNaC::is_a<GiNaC::numeric>(POLYNOMIAL_AFT_SUBS)){return (GiNaC::ex_to<GiNaC::numeric>( POLYNOMIAL_AFT_SUBS ).to_int())%2;}
        POLYNOMIAL_AFT_SUBS = GiNaC::subs(POLYNOMIAL_AFT_SUBS,MONINOMIALS_PAIR->first == MONINOMIALS_PAIR->second);
        //cout<<POLYNOMIAL_AFT_SUBS<<" ==> "<<MONINOMIALS_PAIR->first<<endl;
        POLYNOMIAL_AFT_SUBS = GiNaC::subs(POLYNOMIAL_AFT_SUBS,MONINOMIALS_PAIR->first+REPR == MONINOMIALS_PAIR->second+REPR);
        //cout<<POLYNOMIAL_AFT_SUBS<<" ==> "<<MONINOMIALS_PAIR->first<<endl;
        POLYNOMIAL_AFT_SUBS = GiNaC::subs(POLYNOMIAL_AFT_SUBS,MONINOMIALS_PAIR->first*REPR == (MONINOMIALS_PAIR->second)*REPR).expand();
        //cout<<POLYNOMIAL_AFT_SUBS<<" ==> "<<MONINOMIALS_PAIR->first<<endl;
        if(POLYNOMIAL_AFT_SUBS!=POLYNOMIAL_SUBS){OVER_SUBS=false;} // 发生了替换;
        //if(POLYNOMIAL_AFT_SUBS!=POLYNOMIAL_SUBS){cout<<MONINOMIALS_PAIR->first<<" subs to ==> "<<POLYNOMIAL_AFT_SUBS<<endl;}
        POLYNOMIAL_SUBS = POLYNOMIAL_AFT_SUBS;
    }
    //cout<<" ********* ******** "<<endl;
    //cout<<POLYNOMIAL_AFT_SUBS<<" ===> "<<POLYNOMIAL_SUBS<<endl;
    ppsh_simplify(POLYNOMIAL_SUBS,ITEMS_VARS);
    if (GiNaC::is_a<GiNaC::numeric>(POLYNOMIAL_SUBS)){return (GiNaC::ex_to<GiNaC::numeric>( POLYNOMIAL_SUBS ).to_int())%2;}
    //cout<<POLYNOMIAL_AFT_SUBS<<" ===> "<<POLYNOMIAL_SUBS<<endl;
    if(POLYNOMIAL_AFT_SUBS!=POLYNOMIAL_SUBS){OVER_SUBS=false;} // 发生了简化;
    if(OVER_SUBS){return POLYNOMIAL_SUBS;}
    if(!OVER_SUBS) // 需要继续递归一次;
    {
        POLYNOMIAL_TMP = POLYNOMIAL_SUBS;
        for (int INDEX_MONO = 0; INDEX_MONO < ppsh_size_of(POLYNOMIAL_TMP); ++INDEX_MONO)
        {
            POLYNOMIAL_SUBS = POLYNOMIAL_SUBS - ppsh_index_monomial_of(POLYNOMIAL_TMP,INDEX_MONO) 
                                  + ppsh_subs_by_dict_mono(ppsh_index_monomial_of(POLYNOMIAL_TMP,INDEX_MONO),TAB_MONOMIALS_PAIR,ITEMS_VARS);
        }
        return POLYNOMIAL_SUBS;
    }
}

void ppsh_subs_by_dict(GiNaC::ex & POLYNOMIAL,GiNaC::exmap TAB_MONOMIALS_PAIR,std::vector<GiNaC::symbol> const& ITEMS_VARS)
{
    // for (int INDEX_MONO = 0; INDEX_MONO < POLYNOMIAL.nops(); ++INDEX_MONO)
    // {
    //     POLYNOMIAL = POLYNOMIAL - POLYNOMIAL.op(INDEX_MONO) + ppsh_subs_by_dict_mono(POLYNOMIAL.op(INDEX_MONO),TAB_MONOMIALS_PAIR,ITEMS_VARS);
    // }
    bool OVER_SIMPLYFY = false; GiNaC::ex POLYNOMIAL_TMP;
    GiNaC::ex POLYNOMIAL_SUBS;
    while(!OVER_SIMPLYFY)
    {
        POLYNOMIAL_TMP = POLYNOMIAL;
        for (int INDEX_MONO = 0; INDEX_MONO < ppsh_size_of(POLYNOMIAL_TMP); ++INDEX_MONO)
        {
            POLYNOMIAL = POLYNOMIAL - ppsh_index_monomial_of(POLYNOMIAL_TMP,INDEX_MONO) 
                                    + ppsh_subs_by_dict_mono(ppsh_index_monomial_of(POLYNOMIAL_TMP,INDEX_MONO),TAB_MONOMIALS_PAIR,ITEMS_VARS);
        }
        POLYNOMIAL_SUBS = POLYNOMIAL;
        ppsh_simplify(POLYNOMIAL_SUBS,ITEMS_VARS);
        if(POLYNOMIAL_SUBS==POLYNOMIAL){OVER_SIMPLYFY=true;} // over
    }
    return;
}

// -------------- Gaussian Elimination ---------------
int first_not_0(const std::vector<int>& COEFF)
{
    for (int i = 0; i < COEFF.size(); ++i)
    {
        if (COEFF[i])
        {
            return i;
        }
    }
    return COEFF.size();
}

bool bad_gauss(const std::vector<std::vector<int> >& COEFF_MAT)
{
    for (int i = 0; i < COEFF_MAT.size()-1; ++i)
    {
        if (first_not_0(COEFF_MAT[i])>=first_not_0(COEFF_MAT[i+1]))
        {
            if (first_not_0(COEFF_MAT[i])==first_not_0(COEFF_MAT[i+1]) && first_not_0(COEFF_MAT[i])==COEFF_MAT[i].size()){continue;}
            cout<<i<<endl;
            cout<<first_not_0(COEFF_MAT[i])<<" , "<<first_not_0(COEFF_MAT[i+1])<<endl;
            return true;
        }
    }
    return false;
}

void gaussian_elimination(std::vector<std::vector<int> >& COEFF_MAT)
{
    std::vector<int> TMP_COEFF;
    int NUM_COL = COEFF_MAT[0].size();int NUM_ROW = COEFF_MAT.size();
    int ROW_INDEX = -1;int COL_INDEX = 0;
    bool STILL_THIS_ROW = false;
    // ------------------------ step1:消元 -------------------------------
    while (ROW_INDEX<NUM_ROW-2)
    {   
        if(!STILL_THIS_ROW){ROW_INDEX++;COL_INDEX=ROW_INDEX;}
        if(STILL_THIS_ROW){COL_INDEX++;}
        if(COL_INDEX>=NUM_COL){STILL_THIS_ROW=false;continue;}
        
        // --- swap one non-zero line to top;
        if (COEFF_MAT[ROW_INDEX][COL_INDEX]==0)
        {
            for (int NOT0_ROW_INDEX = ROW_INDEX; NOT0_ROW_INDEX < NUM_ROW; ++NOT0_ROW_INDEX)
            {
                if (COEFF_MAT[NOT0_ROW_INDEX][COL_INDEX]==0){continue;}
                TMP_COEFF = COEFF_MAT[ROW_INDEX];COEFF_MAT[ROW_INDEX] = COEFF_MAT[NOT0_ROW_INDEX];COEFF_MAT[NOT0_ROW_INDEX] = TMP_COEFF; // 交换这两行;
                break;
            }
        }

        if(COEFF_MAT[ROW_INDEX][COL_INDEX]==0){STILL_THIS_ROW=true;continue;} // which means 本列已经全为零;    
        
        // --- elimination;
        for (int NOT0_ROW_INDEX = ROW_INDEX+1; NOT0_ROW_INDEX < NUM_ROW; ++NOT0_ROW_INDEX)
        {
            if(COEFF_MAT[NOT0_ROW_INDEX][COL_INDEX])
            {
                // 对本行进行消元;
                for (int COL_INDEX_ITER = COL_INDEX; COL_INDEX_ITER < NUM_COL; ++COL_INDEX_ITER)
                {   
                    COEFF_MAT[NOT0_ROW_INDEX][COL_INDEX_ITER] = COEFF_MAT[NOT0_ROW_INDEX][COL_INDEX_ITER] ^ COEFF_MAT[ROW_INDEX][COL_INDEX_ITER]; // 0-0 or 1-0/0-1 or 1-1;
                }    
            }        
        }
        STILL_THIS_ROW=false;
    }

    // ---------------------- step2:消元Again -------------------------------   
    cout <<"Once Finished"<<endl;
    int FIRST_NOT_0_ROW = 0;int FIRST_NOT_0_THAT = 0;
    int WHILE_TIMES = 0;
    while(bad_gauss(COEFF_MAT))
    {    
        //if(WHILE_TIMES>10){break;}
        WHILE_TIMES++;
        cout <<"bad!!!!"<<endl;
        for (int ROW_INDEX = 0; ROW_INDEX < NUM_ROW; ++ROW_INDEX)
        {
            FIRST_NOT_0_ROW = first_not_0(COEFF_MAT[ROW_INDEX]);
            for (int NOT0_ROW_INDEX = ROW_INDEX+1; NOT0_ROW_INDEX < NUM_ROW; ++NOT0_ROW_INDEX)
            {
                FIRST_NOT_0_THAT = first_not_0(COEFF_MAT[NOT0_ROW_INDEX]);
                if (FIRST_NOT_0_THAT<FIRST_NOT_0_ROW)
                {
                    TMP_COEFF = COEFF_MAT[ROW_INDEX];COEFF_MAT[ROW_INDEX] = COEFF_MAT[NOT0_ROW_INDEX];COEFF_MAT[NOT0_ROW_INDEX] = TMP_COEFF; // 交换这两行;
                }
            }
        }
        for (int ROW_INDEX = 0; ROW_INDEX < NUM_ROW; ++ROW_INDEX)
        {
            FIRST_NOT_0_ROW = first_not_0(COEFF_MAT[ROW_INDEX]);
            for (int NOT0_ROW_INDEX = ROW_INDEX+1; NOT0_ROW_INDEX < NUM_ROW; ++NOT0_ROW_INDEX)
            {
                FIRST_NOT_0_THAT = first_not_0(COEFF_MAT[NOT0_ROW_INDEX]);
                if (FIRST_NOT_0_THAT==FIRST_NOT_0_ROW)
                {
                    // 对本行进行消元;
                    for (int COL_INDEX_ITER = 0; COL_INDEX_ITER < NUM_COL; ++COL_INDEX_ITER)
                    {   
                        if(COEFF_MAT[NOT0_ROW_INDEX][COL_INDEX_ITER]==0 && COEFF_MAT[ROW_INDEX][COL_INDEX_ITER]==1)
                        {
                            COEFF_MAT[NOT0_ROW_INDEX][COL_INDEX_ITER]=1;
                            continue;
                        }
                        COEFF_MAT[NOT0_ROW_INDEX][COL_INDEX_ITER] = COEFF_MAT[NOT0_ROW_INDEX][COL_INDEX_ITER] - COEFF_MAT[ROW_INDEX][COL_INDEX_ITER]; // 0-0 or 1-0 or 1-1;
                    }    
                }
            }
        }
    }    
}

void ppsh_write_coeff(const std::vector<std::vector<int> >& COEFF_MAT)
{
    ofstream FILEOUT("coeff_ge.mat");
    for (int i = 0; i < COEFF_MAT.size(); ++i)
    {
        std::vector<int> COEFF;
        for (int j = 0; j < COEFF_MAT[0].size(); ++j)
        {
            FILEOUT << COEFF_MAT[i][j];
            FILEOUT << ' ';
        }
        FILEOUT << '\n';
    }
    FILEOUT.close();
}

// =======================  Reading Functions  =======================
std::vector<string> read_formulas(const string& FILE_NAME)
{
    char THIS_LINE[MAX_LEN_LINE];
    string FORMULA_STRING;
    string TMP_VAR;
    std::ifstream FILEIN(FILE_NAME);
    std::vector<string> FORMULA_STRING_SET;
    for (;;)
    {
        FILEIN.getline( THIS_LINE, sizeof(THIS_LINE));
        if ( FILEIN.eof()){break;}
        std::vector<int> COEFF;
        for (int j = 0; j < MAX_LEN_LINE; ++j)
        {
            if ( THIS_LINE[j]=='\0'){break;}
            TMP_VAR = THIS_LINE[j];
            FORMULA_STRING += TMP_VAR;
        }
        FORMULA_STRING_SET.push_back(FORMULA_STRING);
        FORMULA_STRING.clear();
    }
    FILEIN.close();
    return FORMULA_STRING_SET;
}

vector<string> split(const string& PROCESS_STRING, const string& DELIM) {  
    vector<string> res;  
    if("" == PROCESS_STRING) return res;  
    char * strs = new char[PROCESS_STRING.length() + 1] ;
    strcpy(strs, PROCESS_STRING.c_str());   
   
    char * d = new char[DELIM.length() + 1];  
    strcpy(d, DELIM.c_str());  
   
    char *p = strtok(strs, d);  
    while(p) {  
      string s = p; //分割得到的字符串转换为string类型  
      res.push_back(s); //存入结果数组  
      p = strtok(NULL, d);  
    }
    return res;  
}

GiNaC::ex make_poly_from_string(const string& PROCESS_STRING,std::map<string,GiNaC::symbol> MAP_STR_SYMBOL)
{
    std::vector<string> STR_MONOMIALS = split(PROCESS_STRING, " + ");
    GiNaC::ex POLYNOMIAL = 0;
    for (int INDEX_i = 0; INDEX_i < STR_MONOMIALS.size(); ++INDEX_i)
    {
        std::vector<string> STR_SYMBOLS = split(STR_MONOMIALS[INDEX_i], "*");
        if (STR_SYMBOLS.size()==1)
        {
            if(strstr(STR_SYMBOLS[0].c_str(), "^"))
            {
                std::vector<string> VAR_SQUARE = split(STR_SYMBOLS[0], "^");
                POLYNOMIAL+=MAP_STR_SYMBOL[VAR_SQUARE[0] ]*MAP_STR_SYMBOL[VAR_SQUARE[0] ];
                continue;
            }
            if (STR_SYMBOLS[0]=="1")
            {
                POLYNOMIAL+=1;continue;
            }
            POLYNOMIAL+=MAP_STR_SYMBOL[STR_SYMBOLS[0] ];
            continue;
        }
        GiNaC::ex MONOMIAL = 1;
        for (int INDEX_j = 0; INDEX_j < STR_SYMBOLS.size(); ++INDEX_j)
        {
            MONOMIAL *= MAP_STR_SYMBOL[STR_SYMBOLS[INDEX_j] ];
        }
        POLYNOMIAL += MONOMIAL;
    }
    return POLYNOMIAL;
}

std::map<string,GiNaC::symbol> construct_map(const std::vector<GiNaC::symbol>& ITEMS_X)
{
    std::map<string,GiNaC::symbol> MAP_STR_SYMBOL;
    for (int INDEX_i = 0; INDEX_i < ITEMS_X.size(); ++INDEX_i)
    {
        MAP_STR_SYMBOL[ITEMS_X[INDEX_i].get_name()] = ITEMS_X[INDEX_i];
    }
    return MAP_STR_SYMBOL;
}

std::vector<GiNaC::ex> ppsh_read_formulas_from_file(const string& FILE_NAME, std::map<string,GiNaC::symbol> MAP_STR_SYMBOL)
{
    std::vector<GiNaC::ex> EQUATIONS_192;
    std::vector<string> FORMULA_STRING_SET = read_formulas(FILE_NAME);
    for (int INDEX_i = 0; INDEX_i < FORMULA_STRING_SET.size(); ++INDEX_i)
    {
        EQUATIONS_192.push_back(make_poly_from_string(FORMULA_STRING_SET[INDEX_i],MAP_STR_SYMBOL));
    }
    return EQUATIONS_192;
}

GiNaC::exmap ppsh_read_solv_from_file(const string& FILE_NAME, std::vector<GiNaC::symbol> ITEMS_VARS)
{
    GiNaC::exmap TAB_SOLV;
    std::vector<string> FORMULA_STRING_SET = read_formulas(FILE_NAME);
    for (int INDEX_i = 0; INDEX_i < FORMULA_STRING_SET.size(); ++INDEX_i)
    {
        TAB_SOLV[ITEMS_VARS[INDEX_i]] = atoi(FORMULA_STRING_SET[INDEX_i].c_str());
    }
    return TAB_SOLV;
}

bool ppsh_is_sub_poly(GiNaC::ex POLYNOMIAL_SUBS,GiNaC::ex POLYNOMIAL)
{
    std::vector<GiNaC::ex> MONINOMIALS = ppsh_monomials_of(POLYNOMIAL_SUBS);
    for (int INDEX_i = 0; INDEX_i < MONINOMIALS.size(); ++INDEX_i)
    {
        //if (!GiNaC::has(POLYNOMIAL,MONINOMIALS[INDEX_i])){return false;}
        if ( (POLYNOMIAL - MONINOMIALS[INDEX_i] ).nops() > POLYNOMIAL.nops() ){return false;}
    }
    return true;
}

// ==================> PPSH Core Functions END

/* --------------- TEST  ---------------------------------
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
--------------- TEST  ---------------------------------*/
