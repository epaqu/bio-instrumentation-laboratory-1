#include <iostream>
#include <fstream>
#include <string> // getline()
#include <math.h>
#include <sstream> // istringstream


using namespace std;



// Declaration of grobal variables //
int nTotalLineIndex = 0;
int nTheNumberofPMID = 0;

// Lower case function //
void toLowerCaseSTD(std::string&s)
{
for (std::string::iterator i = s.begin(); i != s.end(); ++i)
*i = tolower(*i);
}

// Read file //
ifstream readData(string strFileName)
{
	ifstream f_data(strFileName);
	if (f_data.is_open()){ // f_filenaem.is_open() : check if a file stream was successful opening a file
		cout <<"Reading data success!!"<< endl;
		return f_data;
	}
	else{
		cout<<"Reading data is failed!!"<<endl;
		exit(0);
	}
}
string readLine(ifstream &DataName)
{
	string strLine_local;
	getline(DataName,strLine_local); // getline(f_fileaname, stirng var) : extract characters from f_filename data by line
	toLowerCaseSTD(strLine_local);
	return strLine_local;
}

/*My code from here*/
int countWords(string s)   //NEW FUNCTION made for counting total tokens in a line of the literature.
{
    int word_count(0);
    stringstream ss(s);
    string word;
    while(ss >> word) ++word_count;
    return word_count;
}
/*.............to here*/

void chekingProcess(int nTheNumberofPMID)
{
	if (nTheNumberofPMID % 1 == 0) cout << "The number of processed PMID : "<< nTheNumberofPMID<<endl;
}


string PMIDExtractor(string strLine, string strPreviousPMID)
{
	int nPMIDChecker = strLine.find("pmid");
	string strLocalPMID;
	if (nPMIDChecker == 0){
		strLocalPMID = strLine;
		nTheNumberofPMID++;
		chekingProcess(nTheNumberofPMID);
		cout<<strLocalPMID<<endl;
	} 
	else strLocalPMID = strPreviousPMID;
	return strLocalPMID;
}
string manipulateGeneName(string strGeneToken)
{
	strGeneToken.insert(strGeneToken.begin(),1,' ');
	strGeneToken.insert(strGeneToken.end(),1,' ');
	return strGeneToken;
}

int geneToeknMatching(istringstream &issGeneName, string strGeneToken, int nStandardGeneCheker, string strPubMedLine, ofstream &f_Results,string strPMID)
{
	string strStandardGeneID;

	int tokencount = 0;    //necessary to reset tokencount to zero for each line because when each line of PubMed is read, the tokencount starts from zero/////////////////////         
	
	
	while (getline(issGeneName,strGeneToken,'\t')){

		// Searching standard gene //
		nStandardGeneCheker++;
		if (nStandardGeneCheker == 1) 
		{
			strStandardGeneID = strGeneToken;
			
		}		
		strGeneToken = manipulateGeneName(strGeneToken); // manipulating gene name for exact matching. ex) 'pa010' ==> ' pa010 '  //

		// Matching genes in PubMed sentence //
		if (int(strPubMedLine.find(strGeneToken)) != -1 && (strGeneToken.length() > 5)){ // string.find("word") : return position of first matched word, if there exists matching event to word user definded in string, otherwise return -1

			int nCompareState = 0;
			string strSentenceToken;
			
			/*my codes from here*/                                                                                         //this part of the code counts words in a line until GeneToken, including the gene name
			string shortline=strPubMedLine;                                                                                //declare a new string to store shorten version of PubMedLine 
			shortline = shortline.erase(strPubMedLine.find(strGeneToken)+strGeneToken.length(), strPubMedLine.length());   //PubMedLine is erased from the end of GeneToken to the end of the line. In other words, string shortline is the PubMedLine upto GeneToken.
			for (int i =0; i< countWords(shortline); i++)                                                                  //Now, for i>0 and i<number of words in shortline, we count the words.
			{
				tokencount++;
			}
			/*........... to here*/

			// Sentence tokenization //
			istringstream issPubMedLine(strPubMedLine);		
			
			if (f_Results.is_open()){
				
				f_Results << strStandardGeneID<<"\t\t"<<strPMID<<"\t"<<nTotalLineIndex<<'\t'<<'\t'<<countWords(strPubMedLine)<<'\t'<<'\t'<<'\t'<<'\t'<<tokencount<<'\n'; //DATA VALUES for attributes are added//////////////////
				
			}
			else cout<<"Unable to open result file!!"<<endl;
		}
	}
	return 0;
}

int geneMatching(ifstream &f_GeneDic,string strPubMedLine,ofstream &f_Results,string strPMID)
{
	string strGeneDictionaryLine;

	while (f_GeneDic.good()){
			
		string strStandardGeneID;
		string strGeneToken;
		int nStandardGeneCheker = 0;
		
		// Reading gene dictionary by line //
		strGeneDictionaryLine = readLine(f_GeneDic);

		// Dictionary tokenization //
		istringstream issGeneName(strGeneDictionaryLine);
		geneToeknMatching(issGeneName,strGeneToken,nStandardGeneCheker,strPubMedLine,f_Results,strPMID);

		if (f_GeneDic.eof() == 1){
			f_GeneDic.clear();
			f_GeneDic.seekg(0,ios::beg);
			break;
		}
	}
	return 0;
}

int readAbstract(string strLiteraturefile,string strGeneDicfile,string strOutputfile)
{
	// Read dataset //
	string strFileName_PubMedData = strLiteraturefile;
	string strFileName_GeneDictionary = strGeneDicfile;

	ifstream f_PubMedData = readData(strFileName_PubMedData);
	ifstream f_GeneDic = readData(strFileName_GeneDictionary);

	// Open result data //
	ofstream f_Results (strOutputfile); // ofstream : stream class to write on files
	f_Results<<"GeneName\t\tPMID\tSentenceIndex\tThe#OfTotalTokensInSentence\tLocationOfTokenInSentence\n" ;       //ADDED ATTRIBUTES////////////////////
	
	string strPubMedLine;
	string strPMID;

	while (f_PubMedData.good()){ // f_filenaem.good() : check whether state of stream has error or not

		// Reading PubMed data by line //
		strPubMedLine = readLine(f_PubMedData);
		nTotalLineIndex++;

		

		// PMID Extraction //
		strPMID = PMIDExtractor(strPubMedLine,strPMID);
		geneMatching(f_GeneDic,strPubMedLine,f_Results,strPMID);

	}
		
	f_PubMedData.close();
	f_GeneDic.close();
	f_Results.close();
	return 0;
}

/*
string* Gene(string Table2)
{
	string geneName;
	string* result = new string[256];
	int count = 0;
	bool already = false;
	bool first = true;
	ifstream table2(Table2);
	if (!table2.is_open())
		cout << "Failed reading data 2" << endl;
	else
	{
		if (first)
		{
			table2 >> geneName;
			first = false;
		}
		for (int i = 0; i < 5; i++)
			table2 >> geneName;
		for (int i = 0; i < 256; i++)
		{
			if (geneName == result[i])
				already = true;
		}
		if (!already)
		{
			result[count] = geneName;
			count++;
		}
		already = false;
	}
	return result;
}
*/
/*int* GeneAlone(string* gene, string Table1, string Table2)
{	
	bool found = false;
	bool first = true;
	bool no = false;
	string cmp;
	int pmid, pmidDis;
	int index = 0;
	int* pmidGene = new int[1000];
	int* result = new int[256];
	for (int i = 0; i < 256; i++)
		result[i] = 0;
	ifstream table2(Table2);
	if (!table1.is_open())
		cout << "Failed reading data 1" << endl;
	if (!table2.is_open())
		cout << "Failed reading data 2" << endl;
	else if (table1.is_open())
	{
		if (first)
		{
			table2 >> cmp;
			table2 >> cmp;
			first = false;
		}
		while (table2.good())
		{
			for (int i = 0; i < 4; i++)
				table2 >> cmp;
			table2 >> pmid;
			if (cmp == gene[index])
			{
				for (int i = 0; i < 1000; i++)
				{
					if (pmidGene[i] == pmid)
						found = true;
				}
				if (!found)
				{
					pmidGene[index] = pmid;
					index++;
				}
				found = false;
			}
		}
		ifstream table1(Table1);
		int j = 0;
		bool forfirst = true;
		while (table1.good())
		{
			if (forfirst)
			{
				table1 >> pmidDis;
				table1 >> pmidDis;
				forfirst = false;
			}
			for (int i = 0; i < 5; i++)
				table1 >> pmidDis;
			for (int i = 0; i < 1000; i++)
			{
				if (pmidGene[i] == pmidDis)
					pmidGene[i] = 0;
			}
			if (!no)
			{
				result[j] = result[j]+1;
			}
			no = false;
		}
		for (int i = 0; i < 1000; i++)
		{
			if (pmidGene[i] != 0)
				count++;
		}
		result[ind] = count;
	return result;
}
*/
/*
string** chiScore(string Table1, string Table2)
{
	ifstream table1(Table1);
	if (!table1.is_open())
		cout << "Failed reading data 1" << endl;
	int numBoth=0, numDis=0, numGene=0;
	if (table1.is_open())
	{
		for (int i = 0; i < 5; i++)
		{
			table1 >> word;
		}
		while (table1.good())
		{
			table1 >> DisAtt1;
			table1 >> DisAtt2;
			table1 >> DisAtt3;
			table1 >> DisAtt4;
			table1 >> DisAtt5;
			ifstream table2(Table2);			
			if (!table2.is_open())
				cout << "Failed reading data 2" << endl;
			else
			{
				for (int i = 0; i < 5; i++)
				{
					table2 >> word;
				}
				while (table2.good())
				{
					table2 >> GenAtt1;
					table2 >> GenAtt2;
					table2 >> GenAtt3;
					table2 >> GenAtt4;
					table2 >> GenAtt5;
					if (DisAtt2 == GenAtt2)
					{
						for (int i = 0; i < 256; i++)
						{
							if ((GenAtt2 == passed[i]))
							{
								found = true;
							}
						}
						if (!found)
						{
							numBoth++;
							passed[count] = GenAtt2;
							count++;
						}
						disOnly = false;
					}
				}
				if (disOnly)
					numDis++;
			}
		}
	}
}
*/

string** bestScore(string** table)
{
	string **result;
	result = new string* [20];
	for(int i = 0; i < 20;i++)
	{
		result[i] = new string[2];
	}
	double largest = 0;
	int index = 0;
	string name = "";
	for (int j = 0; j < 20; j ++)
	{
		largest = 0;
		for (int i = 0; i < 256; i++)
		{
			if (atof(table[i][1].c_str()) > largest)
			{
				largest = atof(table[i][1].c_str());
				index = i;
				name = table[i][0];
			}
		}
		table[index][1] = "0";
		result[j][0] = name;
		stringstream ss;
		ss << largest;
		result[j][1] = ss.str();
	}
	ofstream f_Results ("bestScore.txt"); // ofstream : stream class to write on files
	f_Results<<"GeneName\n";
	for (int i = 0; i < 20; i++)
	{
		f_Results << result[i][0]  << endl;
	}
	cout << "Finished writing." << endl;
	return result;
}

string** score(string Table1, string Table2)
{
	double repetition = 0;
	double newScore = 0;
	int Switch = 0;//0 = need new entry
	int count = 0;
	string word, DisAtt1, DisAtt2;
	double DisAtt3, DisAtt4, DisAtt5;
	string GenAtt1, GenAtt2;
	double GenAtt3, GenAtt4, GenAtt5;
	double titleBonus = 1;
	string **result;
	result = new string* [256];
	for(int index=0; index < 256;index++)
	{
		result[index] = new string[3];
	}
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 3; j++)
			result[i][j] = "";
	}
	for (int i = 0; i < 256; i++)
	{
		result[i][2] = "0";
		result[i][1] = "0";
		result[i][0] = "";
	}
	double score = 0;
	ifstream table1(Table1);
	if (!table1.is_open())
		cout << "Failed reading data 1" << endl;
	if (table1.is_open())
	{
		for (int i = 0; i < 5; i++)
		{
			table1 >> word;
		}
		while (table1.good())
		{
			table1 >> DisAtt1;
			table1 >> DisAtt2;
			table1 >> DisAtt3;
			table1 >> DisAtt4;
			table1 >> DisAtt5;
			ifstream table2(Table2);			
			if (!table2.is_open())
				cout << "Failed reading data 2" << endl;
			else
			{
				for (int i = 0; i < 5; i++)
				{
					table2 >> word;
				}
				while (table2.good())
				{
					table2 >> GenAtt1;
					table2 >> GenAtt2;
					table2 >> GenAtt3;
					table2 >> GenAtt4;
					table2 >> GenAtt5;
					if ((DisAtt2 == GenAtt2) && (DisAtt3 == GenAtt3))
					{
						if (DisAtt3 == 2)
							titleBonus = 10;
						//scoring system//
						score = (titleBonus + pow(abs(GenAtt4/(DisAtt5 - GenAtt5)), 2));
						titleBonus = 1;
						//done//
						for (int i = 0; i < 256; i++)
						{
							if ((GenAtt1 == result[i][0]) && (Switch == 0))
							{
								stringstream ss;
								score += atof(result[i][1].c_str());
								ss << score;
								result[i][1] = ss.str();
								Switch = 1;
								repetition = (atof(result[i][2].c_str())+1);
								stringstream sc;
								sc << repetition;
								result[i][2] = sc.str();
							}
						}
						if (Switch == 0)//need new entry
						{
							stringstream ss;
							result[count][0] = GenAtt1;
							ss << score;
							result[count][1] = ss.str();
							result[count][2] = 1;
							repetition = (atof(result[count][2].c_str())+1);
							stringstream sc;
							sc << repetition;
							result[count][2] = sc.str();
							count++;
						}
						Switch = 0;
						score = 0;
					}
				}
			}
		}
		for (int i = 0; i < 256; i++)
		{
			newScore = atof(result[i][1].c_str()) * atof(result[i][2].c_str());
			stringstream ss;
			ss << newScore;
			result[i][1] = ss.str();
		}
	}
	return result;
}
int main ()
{
	string DiseaseTable = "Pancreatic_cancer_Disease_Tagging_table.txt";
	string GeneTable = "Pancreatic_cancer_Gene_Tagging_table.txt";
	bestScore(score(DiseaseTable, GeneTable));
}

