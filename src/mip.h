// Updated at Oct 14, 2021
// Bioinformatics Group, Qingdao University
// Code by: Gongchao Jing, Xiaoquan Su
// Microbial Index of Pathogen
// PM3 version 3.4 or higher

#include <iostream>
#include <fstream> 
#include <sstream>

#include "db.h"
#include "otu_parser.h"
#include "table_format.h"
#include "utility.h"
#include "version.h"

using namespace std;

class _mip{
      
      public: friend class _mip_parser;
             _mip(){
                  Total_otu_count = 0;
                  Pathogen_otu_count = 0;
                  Pathogen_otu_abd = 0;
                  }
      
      private:
            int Total_otu_count;
            int Pathogen_otu_count;
            float Pathogen_otu_abd;
            
            hash_map <string, float, std_string_hash> taxa_info;
			hash_map <string, float, std_string_hash> site_info;
            hash_map <string, float, std_string_hash> infection_info;
			hash_map <string, float, std_string_hash> otu_info;   
                   
      };

class _mip_parser{
      
      public:
      		 _mip_parser(){
			   			string path = Check_Env();
			   			path += "/databases/mip_16s/pathogen/pathogen.dat";
                    
			   			_PMDB ref_db('P'); //default
						Otu_parser = _OTU_Parser(ref_db);
			   			Init(path.c_str());
			   			};
      		
      		_mip Get_pathogen (const char * classfilename, bool is_cp);
      		_mip Get_pathogen (_Table_Format * tablefile, unsigned int n, bool is_cp);
      		
      		void Output_pathogen_summary(const char * summary_file_name, vector <_mip> mips, vector <string> sample_names);
      		void Output_pathogen_taxa(const char * taxa_file_name, vector <_mip> mips, vector <string> sample_names);
      		void Output_pathogen_site(const char * site_file_name, vector <_mip> mips, vector <string> sample_names);
      		void Output_pathogen_infection(const char * infection_file_name, vector <_mip> mips, vector <string> sample_names);      		
      		void Output_pathogen_otu(const char * infection_file_name, vector <_mip> mips, vector <string> sample_names); 
      		
      private:
              _OTU_Parser Otu_parser;
			  hash_set <string, std_string_hash> Pathogen_otus;
			  hash_map <string, string, std_string_hash> Otu_taxa_info;
			  hash_map <string, vector <int>, std_string_hash> Otu_site_info;
              hash_map <string, vector <string>, std_string_hash> Otu_infection_info;
              vector <string> Pathogen_sites; //title
              
              int Init(const char * infilename);
              _mip Get_pathogen (hash_map<string, int, std_string_hash> otu_count, bool is_cp);
      		  _mip Get_pathogen (hash_map<string, int, std_string_hash> otu_count);
      };

int _mip_parser::Init(const char * infilename){
	
	ifstream infile(infilename);
	if (!infile){
		cerr << "Error: Cannot open pathogen otu information, please check then environment variable" << endl;
		exit(0);
		//return 0;
	}
	
	string buffer;
	string str_temp;
	//title
	getline(infile, buffer);
	stringstream strin(buffer);
	strin >> str_temp >> str_temp;
	while(strin >> str_temp)
		Pathogen_sites.push_back(str_temp);
	
	//table	
	while(getline(infile, buffer)){
		stringstream strin(buffer);
		string a_otu, a_taxon;
		vector <int> a_otu_site_info;
		vector <string> a_otu_infection_info;
		strin >> str_temp;
		a_otu = Check_OTU(str_temp);
		strin >> a_taxon;
		
		for (int i = 0; i < Pathogen_sites.size(); i ++){
			strin >> str_temp;
			if (str_temp[0] == '0') 
               a_otu_site_info.push_back(0);
			else{
                 a_otu_site_info.push_back(1);
			     //decode
			     string a_infection = "";
			     for (int j = 0; j <= str_temp.size(); j ++){
                     if ((j >= str_temp.size()) || (str_temp[j] == '|')){
                        if (a_infection.size() > 0) 
                           a_otu_infection_info.push_back("(" + Pathogen_sites[i] + ")_" + a_infection);  
                        a_infection = "";
                     }
                     else 
                          a_infection += str_temp[j];  
                     }
                 
                 }
		}
			
		Pathogen_otus.insert(a_otu);
		Otu_taxa_info[a_otu] = a_taxon;
		Otu_site_info[a_otu] = a_otu_site_info;
		Otu_infection_info[a_otu] = a_otu_infection_info;
	}
		
	infile.close();
	infile.clear();
		
	return Pathogen_otus.size();
	}

_mip _mip_parser:: Get_pathogen (const char * classfilename, bool is_cp){
	
	hash_map<string, int, std_string_hash> otu_count;
	Otu_parser.Load_file_to_hash(classfilename, otu_count);
	
	return Get_pathogen(otu_count, is_cp);
}

_mip _mip_parser:: Get_pathogen (_Table_Format * tablefile, unsigned int n, bool is_cp){
	
	if (n >= tablefile->Get_Sample_Size()) return (_mip());
	
	 vector <string> otus = tablefile->Get_Feature_Names();
	 vector <float> count = tablefile->Get_Abd(n);
	 
	 hash_map<string, int, std_string_hash> otu_count;
	 for (int i = 0; i < otus.size(); i ++)
	 	if (count[i] > 0) otu_count[otus[i]] = (int)count[i];
	 
	 return Get_pathogen(otu_count, is_cp);
}

_mip _mip_parser:: Get_pathogen (hash_map<string, int, std_string_hash> otu_count, bool is_cp){
	
	float total_otu_abd = 0;
	float pathogen_otu_abd = 0;
	
	_mip a_mip;
	
	for (hash_map<string, int, std_string_hash>::iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){
		
		string a_otu = Check_OTU(miter->first);
		int a_otu_seq_count = miter->second;
		
		if (a_otu_seq_count <= 0) continue;
		a_mip.Total_otu_count ++;
		
		float cp_number = 1.0;
		//cp number correction 
		if (is_cp) cp_number = Otu_parser.Get_cp_by_OTU(a_otu);
		total_otu_abd += (float)a_otu_seq_count / cp_number;
		
		if (Pathogen_otus.count(a_otu) != 0){ //a pathogen otu
			a_mip.Pathogen_otu_count ++;
			pathogen_otu_abd += (float)a_otu_seq_count / cp_number;
			
			//parse taxa
			string taxon = Otu_taxa_info[a_otu];
			if (a_mip.taxa_info.count(taxon) == 0) a_mip.taxa_info[taxon] = 0;
			a_mip.taxa_info[taxon] += (float)a_otu_seq_count / cp_number;
			
			//parse site
			vector <int> a_otu_site_info = Otu_site_info[a_otu];
			for (int i = 0; i < a_otu_site_info.size(); i ++){                
                if (a_otu_site_info[i] == 0) continue;
                string site = Pathogen_sites[i];
                if (a_mip.site_info.count(site) == 0) a_mip.site_info[site] = 0;
                a_mip.site_info[site] += (float)a_otu_seq_count / cp_number;
                }
			
			//parse infection
			vector <string> a_otu_infection_info = Otu_infection_info[a_otu];
			for (int i = 0; i < a_otu_infection_info.size(); i ++){
				string infection = a_otu_infection_info[i];
				if (a_mip.infection_info.count(infection) == 0) a_mip.infection_info[infection] = 0;
				a_mip.infection_info[infection] += (float)a_otu_seq_count / cp_number;
				}
			//parse otu
			if (a_mip.otu_info.count(a_otu) == 0) a_mip.otu_info[a_otu] = 0;
			a_mip.otu_info[a_otu] += (float)a_otu_seq_count / cp_number; 
		}	
	}
	
	//normalization
	a_mip.Pathogen_otu_abd = pathogen_otu_abd / total_otu_abd;
	for (hash_map <string, float, std_string_hash> :: iterator miter = a_mip.taxa_info.begin(); miter != a_mip.taxa_info.end(); miter ++)
		miter->second /= total_otu_abd; 
	for (hash_map <string, float, std_string_hash> :: iterator miter = a_mip.site_info.begin(); miter != a_mip.site_info.end(); miter ++)
		miter->second /= total_otu_abd; 
	for (hash_map <string, float, std_string_hash> :: iterator miter = a_mip.infection_info.begin(); miter != a_mip.infection_info.end(); miter ++)
		miter->second /= total_otu_abd; 	
	for (hash_map <string, float, std_string_hash> :: iterator miter = a_mip.otu_info.begin(); miter != a_mip.otu_info.end(); miter ++)
		miter->second /= total_otu_abd; 
	
	return a_mip;
}

_mip _mip_parser:: Get_pathogen(hash_map<string, int, std_string_hash> otu_count){
		return Get_pathogen(otu_count, true);	
}

void _mip_parser:: Output_pathogen_summary(const char * summary_file_name, vector <_mip> mips, vector <string> sample_names){
	
	//summary
	ofstream summary_out(summary_file_name, ios::out);
	if (!summary_out){
		cerr << "Error: Cannot open output file: " << summary_file_name << endl;
		return;
	}
	
	summary_out << "Sample_ID\tMiP\t#_of_OTU\t#_of_Pathogen_OTU" << endl;
	for (int i = 0; i < mips.size(); i ++)
		summary_out << sample_names[i] << "\t" << mips[i].Pathogen_otu_abd << "\t" << mips[i].Total_otu_count << "\t" << mips[i].Pathogen_otu_count << endl;
		
	summary_out.close();
	summary_out.clear();
}

void _mip_parser:: Output_pathogen_taxa(const char * taxa_file_name, vector <_mip> mips, vector <string> sample_names){
	
	//taxa
	vector <string> output_taxa;
	hash_set <string, std_string_hash> output_taxa_hash;
	
	//make taxa table title
	for (int i = 0; i < mips.size(); i ++){
		for (hash_map <string, float, std_string_hash> ::iterator miter = mips[i].taxa_info.begin(); miter !=  mips[i].taxa_info.end(); miter ++)
			if (output_taxa_hash.count(miter->first) == 0) output_taxa_hash.insert(miter->first);
	}
	
	for (hash_set <string, std_string_hash> :: iterator siter = output_taxa_hash.begin(); siter != output_taxa_hash.end(); siter ++)
		output_taxa.push_back(*siter);
 	
 	_Table_Format taxa_table(output_taxa);
 	for (int i = 0; i < mips.size(); i ++){
		vector <float> taxa_abd;
		for (int j = 0; j < output_taxa.size(); j ++)
			if (mips[i].taxa_info.count(output_taxa[j]) != 0) taxa_abd.push_back(mips[i].taxa_info[output_taxa[j]]);
			else taxa_abd.push_back(0);
		taxa_table.Add_Abd(taxa_abd, sample_names[i]);
	}
	taxa_table.Output_Table(taxa_file_name);
}

void _mip_parser:: Output_pathogen_site(const char * site_file_name, vector <_mip> mips, vector <string> sample_names){
	
	//site 
	_Table_Format site_table(Pathogen_sites);
	for (int i = 0; i < mips.size(); i ++){
		vector <float> site_abd;
		for (int j = 0; j < Pathogen_sites.size(); j ++)
			if (mips[i].site_info.count(Pathogen_sites[j]) != 0) site_abd.push_back(mips[i].site_info[Pathogen_sites[j]]);
			else site_abd.push_back(0);
		site_table.Add_Abd(site_abd, sample_names[i]);
		
		//debug
		/* 
		cout << "debug" << endl; 
		for (hash_map <string, float, std_string_hash> :: iterator miter = mips[i].site_info.begin(); miter != mips[i].site_info.end(); miter ++)
			cout << miter->first << "\t" << miter->second << endl; 
		*/ 
	}
	site_table.Output_Table(site_file_name);
}

void _mip_parser:: Output_pathogen_infection(const char * infection_file_name, vector <_mip> mips, vector <string> sample_names){
	
	//infection
	vector <string> output_infections;
	hash_set <string, std_string_hash> output_infection_hash;
	
	//make infection table title
	for (int i = 0; i < mips.size(); i ++){
		for (hash_map <string, float, std_string_hash> ::iterator miter = mips[i].infection_info.begin(); miter !=  mips[i].infection_info.end(); miter ++)
			if (output_infection_hash.count(miter->first) == 0) output_infection_hash.insert(miter->first);
	}
	
	for (hash_set <string, std_string_hash> :: iterator siter = output_infection_hash.begin(); siter != output_infection_hash.end(); siter ++)
		output_infections.push_back(*siter);
 	
 	_Table_Format infection_table(output_infections);
 	for (int i = 0; i < mips.size(); i ++){
		vector <float> infection_abd;
		for (int j = 0; j < output_infections.size(); j ++)
			if (mips[i].infection_info.count(output_infections[j]) != 0) infection_abd.push_back(mips[i].infection_info[output_infections[j]]);
			else infection_abd.push_back(0);
		infection_table.Add_Abd(infection_abd, sample_names[i]);
	}
	infection_table.Output_Table(infection_file_name);
}

void _mip_parser:: Output_pathogen_otu(const char * otu_file_name, vector <_mip> mips, vector <string> sample_names){
	
	//otu
	vector <string> output_otus;
	hash_set <string, std_string_hash> output_otu_hash;
	
	//make otu table title
	for (int i = 0; i < mips.size(); i ++){
		for (hash_map <string, float, std_string_hash> ::iterator miter = mips[i].otu_info.begin(); miter !=  mips[i].otu_info.end(); miter ++)
			if (output_otu_hash.count(miter->first) == 0) output_otu_hash.insert(miter->first);
	}
	
	for (hash_set <string, std_string_hash> :: iterator siter = output_otu_hash.begin(); siter != output_otu_hash.end(); siter ++)
		output_otus.push_back("otu_" + (*siter));
 	
 	_Table_Format otu_table(output_otus);
 	for (int i = 0; i < mips.size(); i ++){
		vector <float> otu_abd;
		for (int j = 0; j < output_otus.size(); j ++){
		    string a_otu = Check_OTU(output_otus[j]);
			if (mips[i].otu_info.count(a_otu) != 0) otu_abd.push_back(mips[i].otu_info[a_otu]);
			else otu_abd.push_back(0);
		}
		otu_table.Add_Abd(otu_abd, sample_names[i]);
	}
	otu_table.Output_Table(otu_file_name);
}
