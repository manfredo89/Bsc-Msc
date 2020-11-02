/*
 * gtheta.h
 *
 *  Created on: Dec 6, 2014
 *      Author: manfredo
 */
#ifndef gtheta_H_
#define gtheta_H_

#include "/home/manfredo/programmi/oxDNA/src/Observables/BaseObservable.h"
#include <sstream>


template<typename number>
class gtheta : public BaseObservable<number> {
private:
	enum {
		CENTRE,
		ARM,
		STICKY
	};
	long int _nconf;
	int _n_bins;
	int _g_theta[500];
	int _get_p_type(BaseParticle<number> *p);
	int _get_rel_center_index(BaseParticle<number> *p);
	
	bool _both(int p_type, int q_type, int to_test) { return (p_type == to_test && q_type == to_test); }
	bool _are_bonded_arms(int, int);
	
public:
	gtheta();
	virtual ~gtheta();
	
	virtual std::string get_output_string(llint curr_step);
	void get_settings (input_file &my_inp, input_file &sim_inp);
};

extern "C" BaseObservable<float> *make_float() { return new gtheta<float>(); }
extern "C" BaseObservable<double> *make_double() { return new gtheta<double>(); }

#endif /* gtheta_H_ */