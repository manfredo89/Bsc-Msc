//
//  input_translator.cpp
//
//
//  Created by Manfredo di Porcia on 19/11/14.
//  Programma per tradurre le configurazione di oxDNA
//  in quelle CG sia top che conf

#include "input_translator.h"


using namespace std;

bool is_arm (int i){
	
	if (i%49 == 42) return true;
	else return false;
}

bool is_sticky (int i){
	
	if (i%49 == 43 || i%49 == 44 || i%49 == 45 || i%49 == 46 || i%49 == 47 || i%49 == 48) return true;
	else return false;
}

bool is_center (int i){
	
	if (i%49 == 20 || i%49 == 21) return true;
	else return false;
}

ifstream& GotoLine(ifstream& file, unsigned int num){
	file.seekg(0, ios::beg);
	for(int i=0; i < num; ++i){
		file.ignore(numeric_limits<streamsize>::max(),'\n');
	}
	return file;
}


int main (int argc, char *argv[]){
	
	
	/******************************************************************************/
	/*** topology file management ***/
	/*** usage is executable file.conf file.top ***/
	
	/* Per costruire la topologia scommentare la prima parte e attenzione alle cose che vengono
	 ridefinite
	 */
	
	/*ifstream top_file (argv[2]);
	try{
		
		if (!top_file.good()) throw string ("Error in opening top file");
		
	}
	catch (const string s){
		
		cout << s.c_str() << endl;
	}
	
	
	char line [64];
	int N_nucleotides, N_strands;
	top_file.getline(line,64);
	sscanf (line, "%d %d\n", &N_nucleotides, &N_strands);
	
	string s = argv[2];
	size_t position = s.find(".top");
	s.insert (position, "_new");
	ofstream new_top_file (s.c_str());
	
	int N_tetramer = N_strands / 4;
	int N_CG_tetramer = 1 + 4 + 4 * 6;
	int N_CG_nucleotides = N_tetramer * N_CG_tetramer;
	new_top_file << N_CG_nucleotides << ' ' << N_tetramer << endl;
	
	char arm_type[9] = "CACGATCG";
	for (int i=0; i<N_tetramer; i++){
		
		new_top_file <<  4 * i + 1 << ' ' << arm_type[0] << ' ' << -1 << ' ' << N_CG_tetramer * i + 1 << endl;
		
		for (int j=0; j<4; j++){
			
			for (int k=0; k<7; k++){
				
				int strand_id = 4 * i + j + 1;
				
				int n3;
				if (k==0 && j!=0) n3 = -1;
				else n3 = (i * N_CG_tetramer) + (j * 7) + k;
				
				int n5;
				if (k==6) n5 = -1;
				else n5 = (i * N_CG_tetramer) + (j * 7) + k + 2;
				
				
				new_top_file << strand_id << ' ' << arm_type[k+1] << ' ' << n3 << ' ' << n5 << endl;
			}
		}
		
	}
	
	top_file.close();
	new_top_file.close(); */
	
	/******************************************************************************/
	/******************************************************************************/
	/*** configuration file management ***/
	
	int N_nucleotides = 19600;
	int N_strands = 400;
	int N_tetramer = N_strands / 4;
	ifstream conf_file (argv[1]);
	
	try{
		
		if (!conf_file.good()) throw string ("Error in opening conf file");
		
	}
	catch (const string s){
		
		cout << s.c_str() << endl;
	}
	
	string s = argv[1];
	size_t position = s.find(".dat");
	s.insert (position, "_new");
	ofstream new_conf_file (s.c_str());
	
	int nlines = 0;
	string file_line;
	while (!conf_file.eof()) {nlines++; getline(conf_file, file_line);}
	
	conf_file.clear();
	conf_file.seekg(0, ios::beg);
	file_line.clear();
	
	for (int i=0; i<3; i++){
		
		getline (conf_file, file_line, '\n');
		new_conf_file << file_line << endl;
	}
	
	double values[15];
	
	int N_in_tetramer = N_nucleotides / N_tetramer;
	for (int k=0; k < N_tetramer; k++){
		
		for (int i=0; i<15; i++) values[i] = 0.;
		conf_file.clear();
		GotoLine (conf_file, 3 + k * N_in_tetramer);
		//center coordinates
		for (int i=0; i<N_in_tetramer; i++){
			for (int j=0; j<15; j++){
				
				conf_file >> file_line;
				if (is_center(i)) values [j] += atof(file_line.c_str());
				
			}
		}
		
		file_line.clear();
		
		for (int i=0; i<15; i++){
			
			file_line.append(to_string(values[i]/8.));
			file_line.append (" ");
		}
		
		new_conf_file << file_line.c_str() << endl;
		
		// arm+sticky coordinates
		
		file_line.clear();
		conf_file.clear();
		GotoLine (conf_file, 3 + k * N_in_tetramer);
		
		for (int i=0; i<N_in_tetramer; i++){
			
			getline (conf_file, file_line, '\n');
			
			if (is_arm(i) || is_sticky(i)) new_conf_file << file_line.c_str() << endl;

			else continue;
		}
	}
	
	conf_file.close();
	new_conf_file.close();
	
}