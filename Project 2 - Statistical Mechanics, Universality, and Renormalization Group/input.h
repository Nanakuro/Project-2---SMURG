#ifndef input_h
#define input_h
#include <vector>
#include <string>
#include <map>
#include <cassert>
#include <fstream>
using namespace std;

class InputTree
{
 public:
  string name;
  vector<InputTree*> sections;
  map<string,string> input;
  InputTree* parent;

};


class InputClass
{
 public:
  int toInteger(string myString)
  {
    return atoi(myString.c_str());
  }
  double toDouble(string myString)
  {
    return atof(myString.c_str());
  }
  bool toBool(string myString)
  {
    if (myString=="True")
      return true;
    else if (myString=="False")
      return false;
    cerr<<"You have a boolean which is neither true nor false"<<endl;
    exit(1);
      
  }
  bool IsVariable(string myVar)
    {
      return tree->input.count(myVar)==1;
    }
  string GetVariable(string myVar)
  {
    if (!IsVariable(myVar)){
      cerr<<"Missing variable "<<myVar<<endl;
    }
    assert(IsVariable(myVar));
    return tree->input[myVar];
  }

  bool OpenSection(string myVar,int num)
  {
    int foundSection=0;
    for (int i=0;i<tree->sections.size();i++){
      if (tree->sections[i]->name==myVar && foundSection==num){
	tree=tree->sections[i];
	return true;
      }
      else if (tree->sections[i]->name==myVar){
	foundSection++;
      }
    }
    return false;
  }

  
  bool OpenSection(string myVar)
  {
    for (int i=0;i<tree->sections.size();i++){
      if (tree->sections[i]->name==myVar){
	tree=tree->sections[i];
	return true;
      }
    }
    return false;
  }

  void CloseSection()
  {
    cerr<<"Closing the section"<<endl;
    cerr<<"The tree is currently "<<tree->name<<" "<<tree->parent->name<<endl;
    tree=tree->parent;
    cerr<<"Now the tree is  "<<tree->name<<endl;
    return;
  }
  InputTree *tree;
  void RemoveTrailingWhiteSpace(string &s)
  {
    int lastPos=s.find_last_not_of(" \t\n\v\f\r");
    if (lastPos!=string::npos){
      s=s.substr(0,lastPos+1);
    }
  }

  void RemoveLeadingWhiteSpace(string &s)
  {
    int firstPos=s.find_first_not_of(" \t\n\v\f\r");
    if (firstPos!=string::npos)
      s=s.substr(firstPos,s.size());
  }
  void RemoveWhiteSpace(string &s){
    RemoveTrailingWhiteSpace(s);
    RemoveLeadingWhiteSpace(s);
  }
  void Read(ifstream &infile)
  { 
    cerr<<"Reading input"<<endl;
    tree=new InputTree();
    tree->name="root";
    InputTree *currentTree;
    currentTree=tree;

    while (!infile.eof()){
      string line;
      infile>>line;
      RemoveWhiteSpace(line);
      if (line[line.size()-1]==':'){
	cerr<<"Reading section "<<endl;
	line=line.substr(0,line.size()-1);
	InputTree *t=new InputTree();
	t->name=line;
	cerr<<"Pushed the tree "<<line<<endl;
	t->parent=currentTree;
	currentTree->sections.push_back(t);
	currentTree=t;
	cerr<<"Done reading section"<<endl;
      }
      else if (line=="END"){
	cerr<<"ending section "<<endl;
	currentTree=currentTree->parent;
	cerr<<"Done ending section"<<endl;
      }
      else if (line.size()==0){
      }
      else {
	int eqLoc=line.find('=');
	assert(eqLoc!=string::npos);
	string token=line.substr(0,eqLoc);
	RemoveWhiteSpace(token);
	string outVals=line.substr(eqLoc+1);
	RemoveWhiteSpace(outVals);
	cerr<<"Putting in variable "<<token<<endl;
	currentTree->input[token]=outVals;
      }
      

    }
    cerr<<"done Reading input"<<endl;
  }
  
  


};


#endif
