/*Lcannon 
make table:[probe,gene,fold change,TSS,CGI,bound,up,down,up5,down5]

inputFileFormat:

    controls tab-delimited
    experiments/\t delimited
    expression-change(float)
    hu/mm
    probelist name
    rma file path name
    windowbed file path name
*/

#include <iostream>
#include <vector>
#include <iterator>
#include <math.h>
#include <map>
using namespace std;


bool fileExists(string& file )
{
    ifstream f(file.c_str());
    return f.is_open();
}

void openFile(string file, ifstream& f )
{
    
 bool exists = fileExists(file);
 if (!exists) {
     cout<< "File " <<file <<"not found." <<"\n";
     return;   
    }
 f.open(file.c_str());   
        
};


int countNames(int x_incriment, string inLine)
{

    istringstream iss(inLine);
    int i = 0;
	  string word;

	   while(iss >> word)
	     {  
	       i = i + 1;
	   }
  return i;
    
};


//////

int main(int argc, char *argv[])
{
    
    ofstream outfile(argv[7], out);
    ifstream inFile;
    inFile.open( argv[1] );
    ofstream outvals(argv[2], out);
    ofstream fcp_table(argv[6], out);
    string line;
    
    int   x = -1;  

    string cgifile;
    string probefile;
    string windowbed;
    string tssfile;

    float up_bound = atof(argv[3]);
    float low_bound = atof(argv[4]);
    float iterations = atof(argv[5]);
    float ex_change;  
   
    getline ( inFile, line);
    x = 0;
    int numControls = countNames ( x, line);
  
    string controls_array[numControls];


//split into array  : control names
    istringstream iss(line);
    int i = 0;//intex token
    string word;
    
    
    while(iss >> word)
     {
        controls_array[i] = word;
        i++;
      }

//split into array : experiment names  
    getline ( inFile, line);
    x = 0;
    int numExperiments = countNames ( x, line);
    string experiments_array[numExperiments];
    istringstream iss1(line);
    i = 0;//intex token
    string word1;

    
    while(iss1 >> word1)
     {
        experiments_array[i] = word1;
        i++;   
     } 

    getline ( inFile, line);
    istringstream iss3(line);
    i = 0;//intex token
    string word_0;
    string rma_file_name;
    
    while(iss3 >> word_0)
     {
        rma_file_name = word_0;
     }
        
//windowbed name:
    getline ( inFile, line);
    istringstream iss_windowbed(line);
    i = 0;//intex token
    string word_windowbed;
    string windowbed_file;
    
    while(iss_windowbed >> word_windowbed)
     {
        windowbed_file = word_windowbed;
        windowbed = word_windowbed;
     }
            
//get cgi file name:
    getline ( inFile, line);
    istringstream iss_cgi(line);
    string word01;
    
    while( iss_cgi >> word01 )
    {
        cgifile = word01;                
    }              
    
//get probe->gene map file name:
    getline ( inFile, line);
    istringstream iss_inprobe(line);
    string word03;
    
    while( iss_inprobe >> word03 )
    {
        probefile = word03;                
    }                                            

//tssfile:
    getline ( inFile, line);
    istringstream iss_tssin(line);
    string word04;
    
    while( iss_tssin >> word04 )
    {
        tssfile = word04;                
    }                                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      


//open affy:
  ifstream rmaFile;
  openFile(rma_file_name, rmaFile);
   
    getline ( rmaFile, line);
    istringstream iss_rma(line);
    i = 0;//intex token
    int k = 0;
    int l = 0;
    string word4;
    int cols_location_of_controls[numControls];
    int cols_location_of_experiments[numExperiments]; 
    vector<int> exp;
    vector<int> con;
    vector<int>::iterator it;

    while(iss_rma >> word4)
    {    i++;
      for( int j = 0; j < numControls; j++)
      {
        if(word4==controls_array[j])
        {      
          cols_location_of_controls[k] = i;
          k++;
        }
      }
        for( int m = 0; m < numExperiments; m++)
          {
            if(word4== experiments_array[m])
            {
              cols_location_of_experiments[l] = i;
              l++;
            }  
          }            
    } 



 //get number of MA probes:
  int countLines = 0;
      
      while(getline ( rmaFile, line ))
      {
            countLines++;
      }
     
      rmaFile.clear();
      rmaFile.seekg(0, rmaFile.beg);
      string probes_ma[countLines + 1];
      


 //calculate differential expression:
    string probe;
    int control_signals[numControls];
    int experiment_signals[numExperiments];
    string newline;
    int madness = 1;
    int row = 0;
    int colCount;
    string p_name;
    float signal;
 
    multimap<string,float>probeExpression;
    multimap<float, string>orderedExpression;
    multimap<string,float>probeAveControl;
    multimap<string,float>probeAveTreatment;

    string probe_key;
    float calculatedSig;
    vector<string> probeVector;
    string probefirst;
    string probelast;
    
    int skipFirstLine = 0;
    float fold_change;
    multimap<float,string>fold;//removed?
while ( getline ( rmaFile, line ) )
{   
  if( skipFirstLine == 0 )
  {
    skipFirstLine++;
   }
  else
  {
    colCount = 0;
    istringstream rma_stream(line);
        rma_stream >> p_name;
        
        if(row==1)
          {
            probefirst = p_name;
          }

        int addControlSig = 0;
        int addExperimentSig = 0;
        float ave_control = 0;
        float ave_experiment = 0;

        while (rma_stream>>signal){

              colCount++;
              for (int n = 0; n <numExperiments; n++)
              {
                if(colCount == cols_location_of_experiments[n])
                {
                  experiment_signals[addExperimentSig] = signal; 
                  ave_experiment = ave_experiment + signal;
                  addExperimentSig++;  
                }
              }
              for (int n = 0; n <numControls; n++)
              {
                if(colCount == cols_location_of_controls[n])
                {
                  control_signals[addControlSig] = signal;
                  ave_control = ave_control + signal;
                  addControlSig++;  
                     
                } 
              } 

            probes_ma[row] = p_name;
            ave_control = ave_control/(addControlSig);
            ave_experiment = ave_experiment/(addExperimentSig);
            calculatedSig = ave_experiment - ave_control;
            fold_change = ave_experiment/ave_control;
            probeExpression.insert(pair<string,float>(p_name,calculatedSig));
            probeAveControl.insert(pair<string,float>(p_name,ave_control));
            probeAveTreatment.insert(pair<string,float>(p_name,ave_experiment));
            ostringstream myString;
            string numbered;

           if(calculatedSig !=-100)
           {
              myString << row;
              myString <<"_"<<calculatedSig<<numbered;
            }

            orderedExpression.insert(pair<float,string>(calculatedSig,p_name));
            fold.insert(pair<float,string>(fold_change,p_name));
            probeVector.push_back(p_name);
            row=row+1;
          
    }        
  } //else        
}//while rma
probelast= p_name;

////////////////TEST a map: 
/* 
 for ( std::multimap< string, string, std::less< int > >::const_iterator iter =uniqueExpression.begin(); iter != uniqueExpression.end(); ++iter )
     {
 //   cout << iter->first << '\t' << iter->second << '\n';
   
  
}

//*set animal file names:*/
string an = experimental_design.animal;

//cgi:
string pline;
ifstream cgi_stream;
openFile(cgifile, cgi_stream);

multimap<string, int> cgiMap;
int cgipluscount=0;
int cgiminuscount=0;
while(getline ( cgi_stream, pline))
{               
  istringstream iss_cginext(pline);
  string mykey;
  iss_cginext >> mykey;
  int myval;
  iss_cginext >> myval;

  if(myval==0)
  {
     cgiminuscount++;
  }
  if(myval==1)
  {
      cgipluscount++;
  }
  
  cgiMap.insert(pair<string,int>(mykey,myval));
}


//probe map:
ifstream probe_stream;
openFile(probefile, probe_stream);

string nexttoken;
multimap<string, string> probeGene;

while( getline( probe_stream, pline) )
{
  string mykey;
  string myval;
  istringstream iss_probe(pline);
  iss_probe >> mykey;
  iss_probe >>myval;
  probeGene.insert(pair<string,string>(mykey,myval));
}


//extract bound info from windowbed file:
ifstream window_stream;
openFile(windowbed, window_stream);
multimap<string,int>boundList;
while( getline( window_stream, pline))
{
     
  string nextcol;
  string genecol;
  string bound_val;
  istringstream iss_window(pline);
  iss_window >> nextcol;
  iss_window >>nextcol;
  iss_window >> nextcol;
  iss_window >>genecol;
  iss_window >> nextcol;
  iss_window >>nextcol;
  iss_window >>bound_val;
         
  string genename;
  istringstream genesplit(genecol);
  getline(genesplit, genename, '|');
             
    if ( bound_val == "0" )
    {
      int  int_bound_val = 0;   
      boundList.insert(pair<string,int>(genename.c_str(),int_bound_val));    
    }
      else{
            int int_bound_val = 1;
            boundList.insert(pair<string,int>(genename.c_str(),int_bound_val));
          }
            genesplit.clear();
}


//TSS map:
std::ifstream tss_stream;
openFile(tssfile, tss_stream);
std::multimap<string, int> tssMap;

while( std::getline( tss_stream, nexttoken) )
{
  string mykey;
  int myval = 1;
  istringstream iss_tss(nexttoken);
  iss_tss >> mykey;
  tssMap.insert(pair<string,int>(mykey,myval));
} 
 

//get up/down expression_change map:
int probeListSize = probeExpression.size(); 
multimap<string,int>down;
multimap<string,int>up;
multimap<float,string>downonly;
multimap<float,string>uponly;
int upcount = 0;
int downcount = 0; 


//get value at percent_up/down boundary:
string nextword;
i = 0;
 for ( std::multimap< float, string, std::less< int > >::const_iterator iter =orderedExpression.begin(); iter != orderedExpression.end(); ++iter )
  {
        string mykey = iter->second;
        float pos_or_neg = iter->first;
        
        int one = 1;
        int zero = 0;
  
        if( pos_or_neg > 0 )
        {
            up.insert(pair<string,int>(mykey,one));
            uponly.insert(pair<float,string>(pos_or_neg,mykey));
            down.insert(pair<string,int>(mykey,zero));
            upcount++;
        }
        if( pos_or_neg <= 0 )
        {
            up.insert(pair<string,int>(mykey,zero));
            down.insert(pair<string,int>(mykey,one));
            downonly.insert(pair<float,string>(pos_or_neg,mykey));
            downcount++;
        }
  i++;          
} 
 
multimap<string,string>geneProbe;
multimap<string,string>orderedGeneProbe;

///get unique gene->probe
map<string,string>uniqueMAgenes;
x=0;

for(std::multimap< string, float, std::less< int > >::const_iterator iter =probeExpression.begin(); 
    iter != probeExpression.end(); ++iter )
{
            
          string mykey = iter->first;
          multimap<string,string>::iterator it = probeGene.lower_bound(mykey);
          multimap<string,string>::iterator it2 = probeGene.upper_bound(mykey);
          string gene_found;
        
          while(it!=it2)
          {
            gene_found = it->second;
            geneProbe.insert(pair<string,string>(mykey,gene_found));
            uniqueMAgenes.insert(pair<string,string>(gene_found,mykey));
            orderedGeneProbe.insert(pair<string,string>(gene_found,mykey));
            it++;
          }  
              
 
}
    
map<string,int>countdowngenes;
map<string,int>countupgenes;    
for(multimap< float, string, std::less< int > >::const_iterator iter =downonly.begin(); iter != downonly.end(); ++iter )
{
            multimap<string,string>::iterator it = probeGene.lower_bound(iter->second);
            multimap<string,string>::iterator it2 = probeGene.upper_bound(iter->second);
            int one = 1;

            while(it!=it2)
            {
              countdowngenes.insert(pair<string,int>(it->second,one));  
              it++;
            }
}
        
        
for(multimap< float, string, std::less< int > >::const_iterator iter =uponly.begin(); iter != uponly.end(); ++iter )
{
            multimap<string,string>::iterator it = probeGene.lower_bound(iter->second);
            multimap<string,string>::iterator it2 = probeGene.upper_bound(iter->second);
            int one = 1;
            
            while(it!=it2)
            {
              countupgenes.insert(pair<string,int>(it->second,one));
              it++;
            }
  }        
      
ofstream x1("probeplus");
ofstream x2("probeminus");

int counti = 0;

multimap<string,string>probeCGIplus;
multimap<string,string>probeCGIminus; 
multimap<float,string>downgenes;
multimap<float,string>upgenes;
multimap<float,string>fc_upplus;
multimap<float,string>fc_upminus;
multimap<float,string>fc_downplus;
multimap<float,string>fc_downminus;



//////maketable//////
outfile<<"probe\tgene\tfoldchange: "<< up_bound<<"-"<<"low_bound"<<"\t"<<"TSS\tCGI\tbound\tup\tdown\tAveTreatment\tAveControl"<<endl;

        for(std::multimap< string, float, std::less< int > >::const_iterator iter =probeExpression.begin(); iter != probeExpression.end(); ++iter )
        {
            
            string myprobe = iter->first;
            float fold_found = iter->second;
            multimap<string,string>::iterator it = geneProbe.lower_bound(myprobe);
            multimap<string,string>::iterator iter2 = geneProbe.upper_bound(myprobe);
            multimap<string,int>::iterator it2 = down.find(myprobe);
            multimap<string,int>::iterator it3 = up.find(myprobe);
            multimap<string,float>::iterator iterator1 = probeAveControl.find(myprobe);
            multimap<string,float>::iterator iterator2 = probeAveTreatment.find(myprobe);
            
            float myAveControl = iterator1 -> second;
            float myAveTreatment = iterator2 -> second;
        
            string gene_found="";
            int down_found;
            int cgi_found;
            int tss_found;
            int up_found;
            int bound_found;
            
            if(it2 != down.end())
            {
              down_found = it2->second;
            }
            
            if(it3 != up.end())
            {
              up_found = it3->second;
            }
            
            while(it != iter2)
            {
              gene_found = it->second;
       
              multimap<string,int>::iterator it6 = cgiMap.lower_bound(gene_found);
              multimap<string,int>::iterator it7 = tssMap.lower_bound(gene_found);
              multimap<string,int>::iterator it9 = boundList.lower_bound(gene_found);
              multimap<string,int>::iterator iter6 = cgiMap.upper_bound(gene_found);
              multimap<string,int>::iterator iter7 = tssMap.upper_bound(gene_found);
              multimap<string,int>::iterator iter9 = boundList.upper_bound(gene_found);
            
              if(it6 != iter6)
              {
                cgi_found  =it6->second;
              }
              if(it7 != iter7)
              {
                tss_found  =it7->second;
                counti++;
                  if(cgi_found ==1)
                  {
                    probeCGIplus.insert(pair<string,string>(myprobe,gene_found));
                        //if up_found..}
                  }
                  if(cgi_found == 0)
                  {
                    probeCGIminus.insert(pair<string,string>(myprobe,gene_found));
                  }  
              }
              else 
              {
                tss_found = 0;   
              }           
            }//while it it2
            
            if(it9 != iter9)
            {
              bound_found  =it9->second;
            }  

            it++;
                
         
            if(gene_found =="")
            {
                outfile<<myprobe<<"\t"<<gene_found<<"\t"<<fold_found<<"\t"<<"N/A"<<"\t"<<"N/A"<<"\t"<<bound_found<<
                "\t"<<up_found<<"\t"<<down_found<<"\t"<<myAveTreatment<<"\t"<<myAveControl<<endl;
   
            }
        
            else if(tss_found ==0)
            {
                outfile<<myprobe<<"\t"<<gene_found<<"\t"<<fold_found<<"\t"<<"N/A"<<"\t"<<"N/A"<<"\t"<<bound_found<<
                "\t"<<up_found<<"\t"<<down_found<<"\t"<<myAveTreatment<<"\t"<<myAveControl<<endl;
   
            }
        
            else
            {
            
                outfile<<myprobe<<"\t"<<gene_found<<"\t"<<fold_found<<"\t"<<tss_found<<"\t"<<cgi_found<<"\t"<<bound_found<<
                "\t"<<up_found<<"\t"<<down_found<<"\t"<<myAveTreatment<<"\t"<<myAveControl<<endl;

            }
          }
    
    

fcp_table<<"%\tdown+ pval\tdown- pval\tFC\tdown+\t""\tdown-\t""\tu+ pval\tu- pval\tFC\tup+\t""\tup-"<<endl;
                     

//diff exvals
int times = up_bound/iterations;

outvals<<rma_file_name<<endl<<windowbed<<endl<<"rows: down+, down-, up+, up-"<<endl
    <<"cols: CGI/bound\t"<<"FC/CGI\t"<<"probeCGI-boundCGI\t"<<"bound/FC_CGI\t"<<
    "pval\texptd_enrchmt"<<"\tEnrich_or_Depleted"<<endl<<"foldchange\t"<<endl;

for( int printout = 0; printout < times; printout++ )
{
    ex_change = low_bound + (iterations*printout);
    if( ex_change <= up_bound)
    {
        fcp_table<<ex_change;
        multimap<string,int>down_percent;
        multimap<string,int>up_percent;
        multimap<string,float>down_percent_value;
        multimap<string,float>up_percent_value;

        int pos_count = (uponly.size())*ex_change;
        int neg_count = (downonly.size())*ex_change;
        int num_up = uponly.size()-((pos_count)/2)-1;
        int num_down = (neg_count)/2;
    }
}
        i=0; 
        float min;
        float lower_bound;

for ( multimap< float, string, std::less< int > >::const_iterator iter3 =downonly.begin(); iter3 != downonly.end(); ++iter3 )
{  
        float value = iter3->first;
        string key = iter3->second;
        if ( i == 0 )
        {
          min = value;
            
        }
        if(i<num_down)
        {

          down_percent_value.insert(pair<string,float>(key,value));

        }
         i++;   
}      

i=0; 
float max;
float upper_bound;
string up_bound_probe;

for ( multimap< float, string, std::less< int > >::const_iterator iter7 =uponly.begin(); iter7 != uponly.end(); ++iter7 )
{  
        float value = iter7->first;
        string key = iter7->second;
        if ( i == 0 )
        {
                up_bound_probe=key;
                min = value;
            
        }
        if(i>num_up)
        {

            up_percent_value.insert(pair<string,float>(key,value));

        }
         i++;   
            
}      


/////get intersection num///// A,B,C,D are abbreviations for the intersecting sets

//D
int downplus;
int downminus;
int upplus;
int upminus;

//C
int maplus;
int maminus;

//B
int boundplus;
int boundminus;

//A
int downplusbound;
int downminusbound;
int upplusbound;
int upminusbound;

map <string,float>gene_fold_down_minus;
map <string,float>gene_fold_down_plus;
map <string,float>gene_fold_up_minus;
map <string,float>gene_fold_up_plus;

map<string,string>cgi_plus_down;
map<string,string>cgi_minus_down;
map<string,string>cgi_plus_up;
map<string,string>cgi_minus_up;


//D down
for(std::multimap< string, float, std::less< int > >::const_iterator iter4 =down_percent_value.begin(); 
    iter4 != down_percent_value.end(); ++iter4 )
{
            string myprobe = iter4->first;
            float expression_val = iter4->second;
            multimap<string,string>::iterator it = probeCGIplus.lower_bound(myprobe); 
            multimap<string,string>::iterator it2 = probeCGIminus.lower_bound(myprobe);
            multimap<string,string>::iterator ittt = probeCGIplus.upper_bound(myprobe);
            multimap<string,string>::iterator itt2 = probeCGIminus.upper_bound(myprobe);
             
            if(it != ittt)
            {
              string match_probe = it->first;
              string mygene = it->second;
              gene_fold_down_plus.insert(pair<string,float>(mygene,expression_val));
              cgi_plus_down.insert(pair<string, string> (mygene, match_probe));
            }
               
            if(it2 != itt2)
            {
              string match_probe = it2->first;
              string mygene = it2->second;
              cgi_minus_down.insert(pair<string, string> (mygene, match_probe));
              gene_fold_down_minus.insert(pair<string,float>(mygene,expression_val));
            }
}

//D up
for(std::multimap< string, float, std::less< int > >::const_iterator iter =up_percent_value.begin();
     iter != up_percent_value.end(); ++iter )
{
            string myprobe = iter->first;
            float expression_val = iter->second;
            multimap<string,string>::iterator it = probeCGIplus.lower_bound(myprobe);
            multimap<string,string>::iterator it2 = probeCGIminus.lower_bound(myprobe);
            multimap<string,string>::iterator ittt = probeCGIplus.upper_bound(myprobe);
            multimap<string,string>::iterator itt2 = probeCGIminus.upper_bound(myprobe);
             
            if(it != ittt)
            {
              string match_probe = it->first;
              string my_gene = it->second;
              cgi_plus_up.insert(pair<string, string> (my_gene, match_probe));
              gene_fold_up_plus.insert(pair<string,float>(my_gene,expression_val));
            }
                
            if(it2 != itt2)
            {
              string match_probe = it2->first;
              string my_gene = it2->second;
              cgi_minus_up.insert(pair<string, string> (my_gene, match_probe));
              gene_fold_up_minus.insert(pair<string,float>(my_gene,expression_val));
            }
}

vector<string>boundplusdown;
vector<string>boundminusdown;
vector<string>boundplusup;
vector<string>boundminusup;

downplus = cgi_plus_down.size();
downminus = cgi_minus_down.size();
upplus = cgi_plus_up.size();
upminus = cgi_minus_up.size();     
   
int bound_count=0;
int count_boundplusdown = 0;   
int count_boundplusup = 0;
int count_boundminusdown = 0;
int count_boundminusup = 0;
                        
//A:       
for(multimap< string, int, std::less< int > >::const_iterator iter =boundList.begin(); iter != boundList.end(); ++iter )
{                                                                       
            string mygene = iter->first;   
            int bound_result = iter->second; 
           
            if(bound_result == 1)
            {
              bound_count++;
              map<string,float>::iterator it = gene_fold_down_plus.find(mygene);
              map<string,float>::iterator it2 = gene_fold_down_minus.find(mygene);
              map<string,float>::iterator it3 = gene_fold_up_plus.find(mygene);
              map<string,float>::iterator it4 = gene_fold_up_minus.find(mygene);
            
              if(it != gene_fold_down_plus.end())
              {                                           
                  count_boundplusdown++;  
                  boundplusdown.push_back(mygene); 
              }    

              if(it2 != gene_fold_down_minus.end())
              {                     
                  count_boundminusdown++;
                  boundminusdown.push_back(mygene);
              }      

              if(it3 != gene_fold_up_plus.end())
              {     
                  count_boundplusup++;
                  boundplusup.push_back(mygene);  
              }  

              if(it4 != gene_fold_up_minus.end())
              {                                           
                  count_boundminusup++;   
                  boundminusup.push_back(mygene);
              } 
            }                                                                                                                                                                                                                                                              
}  
 

//B: chIP cgi+/-
map<string,int>cgiplus_bound;
map<string,int>cgiminus_bound;
int cgiplusbound = 0;
int cgiminusbound = 0;

for(multimap< string, int, std::less< int > >::const_iterator iter =boundList.begin(); iter != boundList.end(); ++iter )
{                                                                       
            string gene_found = iter->first;
            map<string,int>::iterator it2 = tssMap.find(gene_found);
            int bound = iter->second;

            if( bound == 1)
            {
              if(it2 != tssMap.end())
              {  
                map<string,int>::iterator it = cgiMap.find(gene_found);
                  if(it != cgiMap.end())
                  {
                    int value = it->second;
                      if(value==0)
                      {
                        cgiminus_bound.insert(pair<string,int>(iter->first,iter->second));
                        cgiminusbound++;
                      }
                      else
                      {
                        cgiplusbound++;
                        cgiplus_bound.insert(pair<string,int>(iter->first,iter->second));
                      }
                                                                  
                  }        
              }
            }
}  
      
map<string,string>cgiMINUSprobe;
map<string,string>cgiPLUSprobe;

//make unique gene sets for cgi+/- :      
for(multimap< string, string, std::less< int > >::const_iterator iter =probeCGIminus.begin(); iter != probeCGIminus.end(); ++iter )
{                                                                       
    cgiMINUSprobe.insert(pair<string,string>(iter->second,iter->first));
}

for(multimap< string, string, std::less< int > >::const_iterator iter =probeCGIplus.begin(); iter != probeCGIplus.end(); ++iter )
{                                                                       
    cgiPLUSprobe.insert(pair<string,string>(iter->second,iter->first));
}

float min_fold;
float min_fold_boundary;
float max_fold_boundary;
float max_fold;
int get_lower = (num_down)-1;
int get_upper = (orderedExpression.size()-up_percent_value.size())-1;

i=0;

for(multimap<float,string, std::less< int > >::const_iterator iter5 =orderedExpression.begin(); 
      iter5 != orderedExpression.end(); ++iter5 )
{            
            float myval = iter5->first;  
            string key = iter5->second;                                                                                                                    
            
            if(i==0)
            {
              min_fold = myval;
            }

            if(i==(get_lower))
            {
              min_fold_boundary = myval;
              
            }
                
            if(i==get_upper)
            {
              max_fold_boundary = myval;
            }
                
            if(i==fold.size()-1)
            {
              max_fold = myval;
            }
            i++;
}
                              
//from other file i.e. not given as input data                
int cgi_minus_total = 5985;
int cgi_plus_total = 12218;

int sample_space_minus_down = (countdowngenes.size() -cgi_minus_total); 
int sample_space_plus_down = (countdowngenes.size() -cgi_plus_total);
int sample_space_minus_up = (countupgenes.size() -cgi_minus_total); 
int sample_space_plus_up = (countupgenes.size() -cgi_plus_total);

int A_boundplusdown = count_boundplusdown;
int A_boundplusup= count_boundplusup;
int A_boundminusdown = count_boundminusdown;
int A_boundminusup = count_boundminusup;

 
int B_minusbound = cgiminus_bound.size();
int B_plusbound = cgiplus_bound.size();



int C_ma_minus = cgiMINUSprobe.size();
int C_ma_plus = cgiPLUSprobe.size();


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
int D_downplus = gene_fold_down_plus.size();
int D_upplus = gene_fold_up_plus.size();
int D_downminus = gene_fold_down_minus.size();
int D_upminus = gene_fold_up_minus.size();

double B_plus = B_plusbound;
double B_minus = B_minusbound;

double d_plus = D_downplus;
double u_plus = D_upplus;
double d_minus = D_downminus;
double u_minus = D_upminus;



double plusTotalSpace_up = sample_space_plus_up;
double minusTotalSpace_up = sample_space_minus_up;
double plusTotalSpace_down = sample_space_plus_down;
double minusTotalSpace_down = sample_space_minus_down;

float ratioPlus_down = B_plusbound / plusTotalSpace_down;
float ratioMinus_down = B_minusbound /minusTotalSpace_down;

float ratioPlus_up = B_plusbound / plusTotalSpace_up;
float ratioMinus_up = B_minusbound /minusTotalSpace_up;



double expected_downplus = d_plus * ratioPlus_down;
double expected_downminus = d_minus * ratioMinus_down;
double expected_upplus = u_plus * ratioPlus_up;
double expected_upminus = u_minus * ratioMinus_up;


string downplusED;
string downminusED;
string upplusED;
string upminusED;

if (A_boundplusdown > expected_downplus){
    downplusED = "E";
    } //enrichment
if (A_boundplusdown < expected_downplus){
    downplusED = "D";
    } //downregulation 
if (A_boundplusdown == expected_downplus){
    downplusED = "N";
    }//no change

if (A_boundplusup > expected_upplus){
    upplusED = "E";
    }
if (A_boundplusup < expected_upplus){
    upplusED = "D";
    }
if (A_boundplusup == expected_upplus){
    upplusED = "N";
    }

if (A_boundminusdown > expected_downminus){
    downminusED = "E";
    }
if (A_boundminusdown < expected_downminus){
    downminusED = "D";
    }
if (A_boundminusdown == expected_downminus){
    downminusED = "N";
    }
                            
if (A_boundminusup > expected_upminus){
    upminusED = "E";
    }
if (A_boundminusup < expected_upminus){
    upminusED = "D";
    }
if (A_boundminusup == expected_upminus){
    upminusED = "N";
    }

int n1 = abs(B_plusbound+C_ma_plus);
int n2 = abs(B_minusbound+C_ma_minus);


boost::math::hypergeometric_distribution<double> hg_dist(B_plusbound,D_downplus,sample_space_plus_down);
   
boost::math::hypergeometric_distribution<double> hg1_dist(B_plusbound, D_upplus,sample_space_plus_up);

boost::math::hypergeometric_distribution<double> hg2_dist(B_minusbound, D_downminus,sample_space_minus_down);

 boost::math::hypergeometric_distribution<double> hg3_dist(B_minusbound,D_upminus, sample_space_minus_up);
        


outvals<<endl<<"----"<<endl<<ex_change<<endl<<"----"<<endl; 

outvals
    <<B_plusbound<<"\t"<<D_downplus<< "\t"<<sample_space_plus_down<<"\t"<<A_boundplusdown<<"\t"
    << (boost::math::pdf<double>(hg_dist, A_boundplusdown))<<"\t"<<setprecision(4)<<(expected_downplus)<<"\t"<<downplusED<<"\tdown+"<<std::endl
   
     <<B_minusbound<<"\t"<<D_downminus<< "\t"<<sample_space_minus_down<<"\t"<<A_boundminusdown<<"\t"
    <<(boost::math::pdf<double>(hg2_dist, A_boundminusdown))<<"\t"<<setprecision(4)<<(expected_downminus)<<"\t"<<downminusED<<"\tdown-"<< std::endl
   
    <<B_plusbound<<"\t"<<D_upplus<< "\t"<<sample_space_plus_up<<"\t"<<A_boundplusup<<"\t"
    <<(boost::math::pdf<double>(hg1_dist, A_boundplusup))<<"\t"<<setprecision(4)<<(expected_upplus)<<"\t"<<upplusED<<"\tup+"<< std::endl
   
    <<B_minusbound<<"\t"<<D_upminus<< "\t"<<sample_space_minus_up<<"\t"<<A_boundminusup<<"\t"
    <<(boost::math::pdf<double>(hg3_dist, A_boundminusup))<<"\t"<<setprecision(4)<<(expected_upminus)<<"\t"<<upminusED<<"\tup-"<<endl<<"{"<<min_fold_boundary<<", "<<max_fold_boundary<<"}"<< std::endl;  
    
   
    fcp_table<<"\t"<<setprecision(3)<<(boost::math::pdf<double>(hg_dist, A_boundplusdown))
             <<"\t"<<setprecision(3)<<(boost::math::pdf<double>(hg2_dist, A_boundminusdown))
             <<"\t"<<min_fold_boundary
              <<"\t"<<setprecision(5)<<-log10((boost::math::pdf<double>(hg_dist, A_boundplusdown)))<<"\t"<<downplusED
              <<"\t"<<setprecision(5)<<-log10((boost::math::pdf<double>(hg2_dist, A_boundminusdown)))<<"\t"<<downminusED
             // <<"\t"<<min_fold_boundary
             <<"\t"<<setprecision(3)<<(boost::math::pdf<double>(hg1_dist, A_boundplusup))
             <<"\t"<<setprecision(3)<<(boost::math::pdf<double>(hg3_dist, A_boundminusup))
           <<"\t"<<max_fold_boundary
             <<"\t"<<setprecision(5)<<-log10((boost::math::pdf<double>(hg1_dist, A_boundplusup)))<<"\t"<<upplusED
             <<"\t"<<setprecision(5)<<-log10((boost::math::pdf<double>(hg3_dist, A_boundminusup)))<<"\t"<<upminusED
           
             <<endl;
    }

    
  



