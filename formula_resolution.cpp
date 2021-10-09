#include <vector>
#include <string>
#include <iostream>
#include <cstring>
#include <fstream>
#include <ginac/ginac.h>
#include <map>

#define MAX_LEN_LINE 200

using namespace std;


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

int main()
{
  //string s = "z1*c1 + z2*c2 + z3*c3 + z4*c4 + z5*c5 + z6*c6 + z7 + z8";
  //vector<string> v = split(s, " + "); //可按多个字符来分隔;
  //for(vector<string>::size_type i = 0; i != v.size(); ++i)
  //{cout << v[i] << " ";}
  //cout << endl;
  GiNaC::symbol x1("x1"), x2("x2"), x3("x3"), x4("x4"), x5("x5"), x6("x6"), x7("x7"), x8("x8"), x9("x9");
  std::map<string,GiNaC::symbol> MAP_STR_SYMBOL = {{"x1",x1 }, {"x2",x2 }, {"x3",x3 }, {"x4",x4 }, {"x5",x5 }, {"x6",x6 }, {"x7",x7 }, {"x8",x8 }, {"x9",x9 }};
  string STR_POLYNOMIAL = "x1*x2*x3 + x1*x3 + x1 + x2 + x3";
  cout<< make_poly_from_string(STR_POLYNOMIAL,MAP_STR_SYMBOL) <<endl;
  /*std::vector<string> FORMULA_STRING_SET = read_formulas("poly.dat");
  for (int i = 0; i < FORMULA_STRING_SET.size(); ++i)
  {
    cout << FORMULA_STRING_SET[i] <<endl;
  }*/
}
