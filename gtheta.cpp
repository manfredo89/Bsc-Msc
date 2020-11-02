/*
 * gtheta.cpp
 *
 *  Created on: Dec 6, 2014
 *      Author: manfredo
 */

#include "gtheta.h"
#include <iostream>

using namespace std;

template<typename number>
gtheta<number>::gtheta() {
	_nconf = 0;
	_n_bins = 500;
	for (int i=0; i<_n_bins; i++) _g_theta[i] = 0;

}

template<typename number>
gtheta<number>::~gtheta() {
	
}

template<typename number>
int gtheta<number>::_get_p_type(BaseParticle<number> *p) {
	int N_per_tetramer = 29;
	int rel_ind = p->index % N_per_tetramer;
	if(rel_ind == 0) return CENTRE;
	if(((rel_ind-1) % 7) == 0) return ARM;
	
	return STICKY;
}

template <typename number>
int gtheta<number>::_get_rel_center_index(BaseParticle<number> *p){
	int N_per_tetramer = 29;
	int strand_ind = p->index / N_per_tetramer;
	return strand_ind * 29;
}

template<typename number>
bool gtheta<number>::_are_bonded_arms(int p, int q){
	
	double cont = 0.;
	for (int i=p; i<(p+6); i++){
		BaseParticle<number> *part1 = this->_config_info.particles[i];
		
		for (int j=q; j<(q+6); j++){
			BaseParticle<number> *part2 = this->_config_info.particles[j];
			
			number E = this->_config_info.interaction->pair_interaction(part1,part2);
			if (E < 0.) cont += E;
		}
	}
	
	if (cont < -2.){return true;}
	else return false;
}



template<typename number>
std::string gtheta<number>::get_output_string(llint step) {
	int N = *this->_config_info.N;
	_nconf += 1;
	
stringstream ret;
	//number delg = 2. * M_PI / _n_bins;
	for (int i = 0; i < N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		int p_type = _get_p_type(p);
		for (int j = 0; j < i; j++) {
			BaseParticle<number> *q = this->_config_info.particles[j];
			int q_type = _get_p_type(q);
			
			if ((_both(p_type, q_type, ARM)) && ((i/29)!=(j/29))) {
				
				if (_are_bonded_arms(i+1,j+1)){
					
					LR_vector<number> AA_vector = p->pos.minimum_image(q->pos,*this->_config_info.box_side);
					
				
					
					int c_index = _get_rel_center_index (q);
					BaseParticle<number> *c = this->_config_info.particles[c_index];
					LR_vector<number> AC_vector = (q->get_abs_pos(*this->_config_info.box_side) - c->get_abs_pos(*this->_config_info.box_side));
					
					printf ("%lf %lf %d %d %d %d \n", AA_vector.module(), AC_vector.module(), (i/29), (j/29), i, j);
					AA_vector.normalize();

					AC_vector.normalize();
					
					number val = (AA_vector * AC_vector);
					//number val = acos (AA_vector * AC_vector);
					
					//int bin_number = val / delg;
					//_g_theta [bin_number] += 2;
					
					ret << val << endl;
					
				}
				
			}
		}
	}
	
	ret.precision(9);
	
	/*for (int i=0; i<_n_bins; i++){
		
		double alpha = ((double)i + 0.5) * delg;
		double factor = 1. / ((double) _nconf * delg);
		alpha *= 180. / M_PI;
		ret << alpha << " " << _g_theta[i] * factor << endl;
	}*/
	
	return ret.str();
}

template<typename number>
void gtheta<number>::get_settings (input_file &my_inp, input_file &sim_inp) {
	
}

template class gtheta<float>;
template class gtheta<double>;
