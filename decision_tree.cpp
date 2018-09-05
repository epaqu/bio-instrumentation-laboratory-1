#include <iostream>
#include <map>
#include <list>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
struct Data
{
 list<unsigned int> sampleList; // sample list
 list<string> clsGeneList; // candidate gene list for classification
};

struct Class
{
 unsigned int clsResult; // classification (only inputted when the node is leaf node - 0: not determined, 1: normal, 2:disease)
 string clsGene; // gene of best classification (previously defined as int (wrong))
 double cutoffValue; // cutoff value of the clsGene for classification
};

struct Node
{
 Data data; // microarray data (variable name was changed from "d")
 Class cls; // classification information
 Node* high; // son node (expression => cls_value);
 Node* low; // son node (expression < cls_value);
};

class DecisionTree
{
private:
 void drawGraph(Node* n, string pre); // private function for recursive drawing for drawGraph(Node* n) function.

public:
 Node dc; // decision tree
 map<unsigned int, bool> sampleCls; // Sample property: normal(true) vs. cancer(false)
 map<pair<unsigned int, string>, double> exp; // microarray expression data

 /* Prelab6 Work Below*/
 void readDataFromFile(string filename); // read & save microarray data into 'exp' variable of this class.
 list<string> readDEGFromFile(string filename); // read & return differentially expressed gene list.
 double findCutoffValue(string gene, list<unsigned int> sampleList); // find a best cutoff value of classification for an input gene and input samples.
 double getInfoGain(string gene, double cutoffValue, list<unsigned int> sampleList); //calculate information(entropy) for an input gene and sample list.

 /* Functions You should implement in Lab6 */
 void makeDecisionTree(Node* n); // generate Decision tree (recursively executed!)
 bool testSample(unsigned int sampleNo, Node* n); // Test the normal/disease-known sample to the generated Decision tree. -> follow the branching criteria to get final classification and compare with original sample class.
 double getAccuracy(list<unsigned int> sampleList, Node* n); // For given sample list, calculate the average and return the accuracy.
 double nFoldCV(unsigned int n, Data d); // validate with N-fold cross-validation and returns the accuracy

 /* Source-provided Functions */
 void drawGraph(Node* n);
};

void DecisionTree::readDataFromFile(string filename)
{
 //extract and save sample CIs
 ifstream reader(filename);
 bool disease = false;
 int count = 1;
 string line, temp;
 reader >> temp; //dumping the first word (b/c it's unnecessary)
 reader >> temp;
 reader >> line;
 sampleCls.insert(pair<int,bool>(count,false));
 while (temp == line)
 {
  count++;
  if (temp == "non-tumor")
   disease = true;
  sampleCls.insert(pair<int,bool>(count,disease)); //save the disease name
  reader >> line;
 }
 temp = line;
 disease = false; //resetting boolean value back to default so we can test it again
 while (temp == line)
 {
  count++;
  if (temp == "non-tumor")
   disease = true;
  sampleCls.insert(pair<int,bool>(count,disease)); //'non-tumor'
  reader >> line;
 }
 //extract and save exp variables
 while (reader.good())
 {
  temp = line; // storing the gene name in temp;
  reader >> line;
  for (int i = 0; i < count; i++)
  {
   exp.insert(pair<pair<unsigned int, string>, double>(pair<unsigned int, string>(i+1, temp), atof(line.c_str())));
   reader >> line;
  }
 }
};
list<string> DecisionTree::readDEGFromFile(string filename)
{
 list<string> result;
 ifstream reader(filename);
 string line;
 while (reader >> line)
 {
  result.push_back(string(line)); //adding a gene name to the list.
  reader >> line; //skipping the p-value
 }
 return result;
};
double DecisionTree::findCutoffValue(string gene, list<unsigned int> sampleList)
{
 list<unsigned int>::iterator reader;
 double totNorm=0, totDis=0, numNorm=0, numDis=0;
 for(reader= sampleList.begin(); reader!=sampleList.end(); reader++)
 {
  if (sampleCls.find(*reader)->second == true) //if normal
  {
   numNorm+=1; //counting the number of normal samples
   totNorm += exp.find(pair<unsigned int, string>(*reader, gene))->second; //summing up all expression values of normal samples
  }
  else //if cancer
  {
   numDis+=1; //counting the number of cancer samples
   totDis+=exp.find(pair<unsigned int, string>(*reader, gene))->second; //summing up all expression values of cancer samples
  }
 }
 return ((totNorm/numNorm)+(totDis/numDis))/2;
};
double calculator(double normal, double cancer, double lowNorm, double highNorm, double lowDis, double highDis)
{
 //long and tedious calculation of InfoGain
 double sum = normal+cancer;
 double lowSum = lowNorm+lowDis;
 double highSum = highNorm+highDis;
 double entropy = ((normal/sum)*(log(normal/sum)/log(2.0))) - ((cancer/sum)*(log(cancer/sum)/log(2.0)));
 double lowEntropy = ((lowNorm/lowSum)*(log(lowNorm/lowSum)/log(2.0))) - ((lowDis/lowSum)*(log(lowDis/lowSum)/log(2.0)));
 double highEntropy = ((highNorm/highSum)*(log(highNorm/highSum)/log(2.0))) - ((highDis/highSum)*(log(highDis/highSum)/log(2.0)));
 double result = entropy - (lowNorm+lowDis)/(sum)*lowEntropy - (highNorm+highDis)/(sum)*highEntropy;
 return result;
}
double DecisionTree::getInfoGain(string gene, double cutoffValue, list<unsigned int> sampleList)
{
 list<unsigned int>::iterator reader;
 double numNorm=0, lowerNorm=0,larNorm=0, numDis=0,lowerDis=0,larDis =0;
 for (reader= sampleList.begin(); reader!=sampleList.end(); reader++)
 {
  if (sampleCls.find(*reader)->second  == true) //if normal
   numNorm++; //counting the # of normal samples
  else
   numDis++; //counting # of cancer samples
  if (exp.find(pair<unsigned int, string>(*reader, gene))->second < cutoffValue) //if sample's exp value is lower than cutoff
  {
   if(sampleCls.find(*reader)->second == true) //if normal
    lowerNorm++; //counting the # of normal samples with lower exp values
   else
    lowerDis++; //counting the # of cancer samples with lower exp values
  }
  else //if sample's exp value is greater than or equal to cutoff
  {
   if(sampleCls.find(*reader)->second == true) //if normal
    larNorm++; //counting the # of normal samples with greater (or equal) exp values
   else
    larDis++; //counting the # of cancer samples with greater (or equal) exp values
  }  
 }
 double infogain = calculator (numNorm, numDis, lowerNorm, larNorm, lowerDis, larDis);
 return infogain;
};
void DecisionTree::drawGraph(Node* n){
 drawGraph(n, "");
};
void DecisionTree::drawGraph(Node* n, string pre){
 stringstream ss;
 ss<<(n->cls.cutoffValue);
 
 if(n->cls.clsResult!=0){
  cout << "─[" << (n->cls.clsResult==1 ? "Normal" : "Cancer") << "]" << endl;
 }else{
  cout << "┬" + n->cls.clsGene + "(" + ss.str() + ")" << endl;
  cout << pre << "├>";
  drawGraph(n->high, pre+"│ ");
  cout << pre << "└<";
  drawGraph(n->low, pre+"   ");
 }
};

void DecisionTree::makeDecisionTree(Node* n)
{
 list<unsigned int>::iterator reader;
 list<unsigned int>::iterator reader2;
 list<unsigned int> listsample = n->data.sampleList;
 bool addNode = false;
 for(reader = listsample.begin(); reader!=listsample.end(); reader++)
 {
  if (sampleCls.find(*reader)->second == true)
  {    
   for (reader2 = listsample.begin(); reader2!=listsample.end(); reader2++)
   {
    if (sampleCls.find(*reader2)->second == false)
    {
     n->cls.clsResult = 0;
     addNode = true; //we need new nodes;
    }
   }
   if (addNode == false)
   {
    n->cls.clsResult = 1;
    return;
   }
  }
  else//if first sample is tumor
  {
   for (reader2 = listsample.begin(); reader2!=listsample.end(); reader2++)
   {
    if (sampleCls.find(*reader2)->second == true)
    {
     n->cls.clsResult = 0;
     addNode = true; // we need new nodes;
    }
   }
   if (addNode == false)
   {
    n->cls.clsResult = 2;
    return;
   }
  }
  if (addNode == true) // we need new nodes
  {
   string bestGene;
   double temporary = 0;
   double highestScore = 0;
   list<string>::iterator reader3;
   list<string> geneList = n->data.clsGeneList;
   for (reader3 = geneList.begin(); reader3 != geneList.end(); reader3++)
   {
    temporary = getInfoGain(*reader3, findCutoffValue(*reader3, listsample), listsample);
    if (highestScore < temporary)
    {
     highestScore = temporary;
     bestGene = *reader3;
    }
   }
   n->cls.clsGene = bestGene;
   n->cls.cutoffValue = findCutoffValue(bestGene, listsample);
   Node lowN, highN;
   list<unsigned int>:: iterator itSampleCl;
   for(itSampleCl = listsample.begin(); itSampleCl != listsample.end(); itSampleCl++)
   {
    //sample list is passed on
    if (exp.find(pair<unsigned int, string>(*itSampleCl, n->cls.clsGene))->second > n->cls.cutoffValue)
     highN.data.sampleList.push_back(*itSampleCl);
    else
     lowN.data.sampleList.push_back(*itSampleCl);
   }
   geneList.remove(n->cls.clsGene);//geneList is passed on
   highN.data.clsGeneList = geneList;
   lowN.data.clsGeneList = geneList;
   n->high = &highN; //parental relationship
   n->low = &lowN;
   addNode = false;
   makeDecisionTree(n->high);
   makeDecisionTree(n->low);
   return;
  }
 }
};
bool DecisionTree::testSample(unsigned int sampleNo, Node* n)
{
 // 1. Check node's classification value(clsResult) is not 0, then return true for 1(normal) or return false for 2(tumor).
 if (n->cls.clsResult == 1)
 {
  //decision tree says it's normal
  if (sampleCls.find(sampleNo)->second == true)
   return true;
  else
   return false;
 }
 else if (n->cls.clsResult == 2)
 {
  //decision tree says it's normal
  if (sampleCls.find(sampleNo)->second == true)
   return false;
  else
   return true;
 }
 else
 {
  if (exp.find(pair<unsigned int, string>(sampleNo, n->cls.clsGene))->second > n->cls.cutoffValue)
   testSample(sampleNo, n->high);
  else
   testSample(sampleNo, n->low);
 }
};
double DecisionTree::getAccuracy(list<unsigned int> sampleList, Node* n)
{
 double count, precision = 0;
 list<unsigned int>::iterator reader;
 for (reader = sampleList.begin(); reader != sampleList.end(); reader++)
 {
  if (testSample(*reader, n) == true)
   precision++;
  count++;
 }
 return precision / count;
};
double DecisionTree::nFoldCV(unsigned int n, Data d)
{
 int average = 0;
 list<unsigned int> sampleNum = d.sampleList;
 int size = sampleNum.size();
 int numElements = size/n;
 vector<int> normal;
 vector<int> disease;
 for (int i = 0; i < size; i++)
 {
  if (sampleCls.find(i+1)->second == true)
  {
   normal.push_back(sampleNum.front());
   sampleNum.pop_front;
  }
  else
  {
   disease.push_back(sampleNum.front());
   sampleNum.pop_front;
  }
 }
 random_shuffle(normal.begin(), normal.end());
 random_shuffle(disease.begin(), disease.end());

 list<unsigned int> subset;
 map<unsigned int, list<unsigned int>> subsets;
 int normalsize = normal.size();
 int cancersize = disease.size();
 int mapsize;
 for (int i = 1; i < n+1; i++)
 {
  for (int i = 0; i < normalsize/n; i++)
  {
   subset.push_back(normal.back);
   normal.pop_back();
  }
  for (int i = 0; i < cancersize/n; i++)
  {
   subset.push_back(disease.back);
   disease.pop_back();
  }
  subsets.insert(pair<unsigned int,list<unsigned int>>(i, subset));
  while(!subset.empty())
   subset.pop_back;
 }
 if (!normal.empty)
 {
  for (int i = 1; i < normal.size()+1; i++)
  {
   subsets.find(i)->second.push_back(normal.back);
   normal.pop_back();
  }
 }
 if (!disease.empty)
 {
  mapsize = subsets.size();
  for (int i = 1; i < disease.size()+1; i++)
  {
   subsets.find(mapsize)->second.push_back(disease.back);
   disease.pop_back();
   mapsize--;
  }
 }

 mapsize = subsets.size();
 list <unsigned int> total;
 for (int i = 1; i < mapsize+1; i++)
 {
  for (int j = 0; j < subsets.find(i)->second.size(); j++)
  {
   total.push_back(subsets.find(i)->second);
   subsets.find(i)->second.pop_back();
  }
 }

   n-1 subset node
  makeDecisionTree(node);
  average += getAccuracy(sampleList of the one subset, root of decision tree);
 }
 return average / n;
};
void main()
{
 // 1. Read data from file (prelab6)
 DecisionTree dt;
 dt.readDataFromFile("EXAMPLE_splitted_dataset_example(astocytoma).txt");
 list<string> deg = dt.readDEGFromFile("1558201_s_at.txt");


 // 2. Prepare data part of root node for a parameter of function that generates decision tree
 Node n;
 map<unsigned int, bool>::iterator it3;
 for(it3 = dt.sampleCls.begin(); it3 != dt.sampleCls.end(); it3++)
 {
  n.data.sampleList.push_back(it3->first);
 }
 n.data.clsGeneList=deg;


 // 3. Generate Decision Tree & Draw the result
 dt.makeDecisionTree(&n);
 dt.drawGraph(&n);


 // 4. Prepare whole data for cross validation
 Data d;
 d.clsGeneList = n.data.clsGeneList;
 d.sampleList = n.data.sampleList;

 // 5. Print 10-fold cross validation result
 cout << dt.nFoldCV(10, d);
 system("PAUSE");
}

