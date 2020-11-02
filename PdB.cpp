/*
 * PdB.cpp
 *
 *  Created on: Dec 6, 2014
 *      Author: manfredo
 */

#include "PdB.h"

using namespace std;

template<typename number>
PdB<number>::PdB() {
	_nconf = 0;
}

template<typename number>
PdB<number>::~PdB() {
	
}

template<typename number>
int PdB<number>::_get_p_type(BaseParticle<number> *p) {
	int N_per_tetramer = 29;
	int rel_ind = p->index % N_per_tetramer;
	if(rel_ind == 0) return CENTRE;
	if(((rel_ind-1) % 7) == 0) return ARM;
	
	return STICKY;
}

/*template<typename number>
int PdB<number>::_get_p_arm(BaseParticle<number> *p) {
	int N_per_tetramer = 29;
	int rel_ind = p->index % N_per_tetramer;
	if((rel_ind > 0) && (rel_ind < 8)) return 1;
	if((rel_ind > 7) && (rel_ind < 15)) return 2;
	if((rel_ind > 14) && (rel_ind < 22)) return 3;
	if((rel_ind > 21) && (rel_ind < 29)) return 4;
		
	return -1;
}*/

template<typename number>
bool PdB<number>::_are_bonded_arms(int p, int q){
	
	int cont = 0;
	for (int i=p; i<(p+6); i++){
		BaseParticle<number> *part1 = this->_config_info.particles[i];
		
		for (int j=q; j<(q+6); j++){
			BaseParticle<number> *part2 = this->_config_info.particles[j];
			
			number E = this->_config_info.interaction->pair_interaction(part1,part2);
			if (E<0.) cont++;
		}
	}
	
	if (cont > 2) return true;
	else return false;
}



template<typename number>
std::string PdB<number>::get_output_string(llint step) {
	int N = *this->_config_info.N;
	_nconf += 1;
	
	int N_bonds = 0;
	for (int i = 0; i < N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		int p_type = _get_p_type(p);
		for (int j = 0; j < i; j++) {
			BaseParticle<number> *q = this->_config_info.particles[j];
			int q_type = _get_p_type(q);
			
			if ((_both(p_type, q_type, ARM)) && ((i/29)!=(j/29))) {
				
				if (_are_bonded_arms(i+1,j+1)) N_bonds++;
				
			}
		}
	}
	
	stringstream ret;
	ret.precision(9);
	
	ret << _nconf << " " << N_bonds << endl;
	
	return ret.str();
}

template<typename number>
void PdB<number>::get_settings (input_file &my_inp, input_file &sim_inp) {
	
}

template class PdB<float>;
template class PdB<double>;
