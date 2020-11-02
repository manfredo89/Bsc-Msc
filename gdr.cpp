/*
 * GdR.cpp
 *
 *  Created on: Dec 4, 2014
 *      Author: manfredo
 *
 *  Ultimo file usato per calcolare le g(r) in oxDNA
 *
 *
 */

#include "gdr.h"

using namespace std;

template<typename number>
GdR<number>::GdR() {
	_nconf = 0;
	_nbins = -1;
	_bin_size = (number) -1.;
}

template<typename number>
GdR<number>::~GdR() {
	
}

template<typename number>
int GdR<number>::_get_p_type(BaseParticle<number> *p) {
	int N_per_tetramer = 29;
	int rel_ind = p->index % N_per_tetramer;
	if(rel_ind == 0) return CENTRE;
	if(((rel_ind-1) % 7) == 0) return ARM;
	
	return -10000000;
}

template<typename number>
int GdR<number>::_get_index(BaseParticle<number> *p) {
	int N_per_tetramer = 29;
	int rel_ind = p->index / N_per_tetramer;
	
	return rel_ind;
}

template<typename number>
std::string GdR<number>::get_output_string(llint step) {
	int N = *this->_config_info.N;
	number box_side = *this->_config_info.box_side;
	_nconf += 1;
	
if (_max_value > box_side / 2. && _nconf == 1) OX_LOG(Logger::LOG_WARNING, "Observable Rdf: computing profile with max_value > box_size/2. (%g > %g/2.)", _max_value, box_side);
	
	for (int i = 0; i < N; i ++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		int p_type = _get_p_type(p);
		for (int j = 0; j < i; j ++) {
			BaseParticle<number> *q = this->_config_info.particles[j];
			int q_type = _get_p_type(q);
			
			if (_get_index(p) != _get_index(q)){
				
				if (_both(p_type, q_type, CENTRE)) {
					LR_vector <number> dr = p->pos.minimum_image(q->pos, box_side);
					dr = LR_vector<number> (dr.x, dr.y, dr.z);
					number drmod = dr.module();
					if (drmod < _max_value){
						int mybin = (int) (0.01 + floor (drmod / _bin_size));
						_gr_center_center[mybin] += 2;
					}
				}
				if (_any(p_type, q_type, CENTRE) && _any(p_type, q_type, ARM)) {
					LR_vector <number> dr = p->pos.minimum_image(q->pos, box_side);
					dr = LR_vector<number> (dr.x, dr.y, dr.z);
					number drmod = dr.module();
					if (drmod < _max_value){
						int mybin = (int) (0.01 + floor (drmod / _bin_size));
						_gr_arm_center[mybin] += 2;
					}
				}
				if (_both(p_type, q_type, ARM)) {
					LR_vector <number> dr = p->pos.minimum_image(q->pos, box_side);
					dr = LR_vector<number> (dr.x, dr.y, dr.z);
					number drmod = dr.module();
					if (drmod < _max_value){
						int mybin = (int) (0.01 + floor (drmod / _bin_size));
						_gr_arm_arm[mybin] += 2;
					}
				}
			}
		}
	}
	
	stringstream ret;
	ret.precision(9);
	int N_tetramer = N/29;
	number fact = box_side * box_side * box_side / (4 * M_PI * _bin_size * _nconf * N_tetramer * (N_tetramer - 1.));
	double myx = _bin_size / 2.;
	for (int i=0; i<_nbins; i++){
		
		double norm =  fact / (myx * myx);
		ret << myx << " " << _gr_center_center[i] * norm << " " << _gr_arm_center[i] * norm / 8. << " " << _gr_arm_arm[i] * norm / 16. << endl;
		myx += _bin_size;
		
	}
	ret << endl;
	
	return ret.str();
}

template<typename number>
void GdR<number>::get_settings (input_file &my_inp, input_file &sim_inp) {
	
	float tmpf;
	getInputFloat(&my_inp, "max_value", &tmpf, 1);
	_max_value = (number) tmpf;
	getInputFloat(&my_inp, "bin_size", &tmpf, 1);
	_nbins = (int) (floor(_max_value / tmpf) + 0.01);
	_bin_size = (number) _max_value / _nbins;
	
	OX_LOG(Logger::LOG_INFO, "Observable GdR initialized with  bin_size %g (%g), nbins %d", _bin_size, tmpf, _nbins);
	
	_gr_center_center.resize(_nbins);
	_gr_arm_center.resize(_nbins);
	_gr_arm_arm.resize(_nbins);
}

template class GdR<float>;
template class GdR<double>;
