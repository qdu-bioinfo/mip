// Updated at Apr 18, 2019
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Code by: Gongchao Jing, Xiaoquan Su
// Microbial Index of Pathogen
// PM3 version 3.4 or higher

#include "mip.h"

char Ref_db;

//single input
string Infilename;

//list input
string Listfilename;
string Listprefix;

//table input
string Tablefilename;

string Outprefix;

int Coren = 0;

bool Is_cp_correct = true;
int Mode = 0; //0: single, 1: multi_list, 2: multi_table

int printhelp(){
    
    cout << "Parse the MiP version: " << Version << endl;
    cout << "\tCompute the Microbal Index of Pathogen (16S only)" << endl;
    cout << "Usage:" << endl;
    cout << "PM-parse-mip [Option] Value" << endl;
    cout << "Options: " << endl;
    
    cout << "\t[Input options, required]" << endl;
    cout << "\t  -i Input single file" << endl;
    cout << "\tor" << endl;
    cout << "\t  -l Input files list" << endl;
    cout << "\t  -p List file path prefix [Optional for -l]" << endl;
    cout << "\tor" << endl;
    cout << "\t  -T (upper) Input OTU count table (*.OTU.Count)" << endl;
    
    cout << "\t[Output options]" << endl;
    cout << "\t  -o Output prefix, default is \"MiP\"" << endl;
    
    cout << "\t[Other options]" << endl;
    cout << "\t  -r rRNA copy number correction, T(rue) or F(alse), default is T" << endl; 
    //cout << "\t  -t Cpu core number, default is auto" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
    return 0;
    }

int Parse_Para(int argc, char * argv[]){
    //init
	Outprefix = "MiP";
	Is_cp_correct = true;
	Mode = 0;
	    
    if (argc ==1) 
		printhelp();
    
    int i = 1;
    
    while(i<argc){
         if (argv[i][0] != '-') {
                           cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
                           exit(0);
                           };           
         switch(argv[i][1]){
                            //case 'D': Ref_db = argv[i+1][0]; break;
                            //single input
                            case 'i': Infilename = argv[i+1]; Mode = 0; break;
                            //list input
							case 'l': Listfilename = argv[i+1]; Mode = 1; break;
                            case 'p': Listprefix = argv[i+1]; break;
                            //table inpuut
							case 'T': Tablefilename = argv[i+1]; Mode = 2; break;
                            
							case 'o': Outprefix = argv[i+1]; break;
                            
							case 'r': if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F')) Is_cp_correct = false; break;
                            //case 't': Coren = atoi(argv[i+1]); break;         
                            case 'h': printhelp(); break;
                            default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break; 
                            }
         i+=2;
         }
    
	return 0;
    }

int main(int argc, char * argv[]){
	
	Parse_Para(argc, argv);
	
	_mip_parser mip_parser;
	vector <_mip> mips;
	vector <string> sam_name;
	
	cout << "Pathogen analysis starts" << endl;
	
	switch (Mode){
		case 0: {
				sam_name.push_back("sample_1");
				mips.push_back(mip_parser.Get_pathogen(Infilename.c_str(), Is_cp_correct));
				break;
				}
		case 1: {
				vector <string> file_list;
				int file_count = Load_List(Listfilename.c_str(), file_list, sam_name, Listprefix);
				for (int i = 0; i < file_count; i ++)
					mips.push_back(mip_parser.Get_pathogen(file_list[i].c_str(), Is_cp_correct));
				break;
				}
		case 2: {
				_Table_Format tablefile(Tablefilename.c_str());
				sam_name = tablefile.Get_Sample_Names();
				for (int i = 0; i < tablefile.Get_Sample_Size(); i ++)
					mips.push_back(mip_parser.Get_pathogen(&tablefile, i, Is_cp_correct));
				break;
				}		
		default: break;
	}
	
	mip_parser.Output_pathogen_summary((Outprefix + ".summary.out").c_str(), mips, sam_name);
	mip_parser.Output_pathogen_taxa((Outprefix + ".taxa.out").c_str(), mips, sam_name);
	mip_parser.Output_pathogen_site((Outprefix + ".site.out").c_str(), mips, sam_name);
	mip_parser.Output_pathogen_infection((Outprefix + ".infection.out").c_str(), mips, sam_name);
	mip_parser.Output_pathogen_otu((Outprefix + ".OTU.Abd.out").c_str(), mips, sam_name);
	
	cout << "Pathogen analysis finished" << endl;
	
	return 0;
}
