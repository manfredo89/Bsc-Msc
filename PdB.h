/*
 * PdB.h
 *
 *  Created on: Dec 6, 2014
 *      Author: manfredo
 */
#ifndef PDB_H_
#define PDB_H_

#include "/home/manfredo/programmi/oxDNA/src/Observables/BaseObservable.h"
#include <sstream>


template<typename number>
class PdB : public BaseObservable<number> {
private:
	enum {
		CENTRE,
		ARM,
		STICKY
	};
	long int _nconf;
	int _get_p_type(BaseParticle<number> *p);
	
	bool _both(int p_type, int q_type, int to_test) { return (p_type == to_test && q_type == to_test); }
	bool _are_bonded_arms(int, int);
	
public:
	PdB();
	virtual ~PdB();
	
	virtual std::string get_output_string(llint curr_step);
	void get_settings (input_file &my_inp, input_file &sim_inp);
};

extern "C" BaseObservable<float> *make_float() { return new PdB<float>(); }
extern "C" BaseObservable<double> *make_double() { return new PdB<double>(); }

#endif /* PDB_H_ */