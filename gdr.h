/*
 * GdR.h
 *
 *  Created on: Dec 4, 2014
 *      Author: manfredo
 */
#ifndef GDR_H_
#define GDR_H_

#include "/home/manfredo/programmi/oxDNA/src/Observables/BaseObservable.h"
#include <sstream>


template<typename number>
class GdR : public BaseObservable<number> {
private:
	enum {
		CENTRE,
		ARM,
		STICKY
	};
	long int _nconf;
	number _max_value;
	number _bin_size;
	int _nbins;
	std::vector<long int> _gr_center_center;
	std::vector<long int> _gr_arm_center;
	std::vector<long int> _gr_arm_arm;
	int _get_p_type(BaseParticle<number> *p);
	
	bool _any(int p_type, int q_type, int to_test) { return (p_type == to_test || q_type == to_test); }
	bool _both(int p_type, int q_type, int to_test) { return (p_type == to_test && q_type == to_test); }
	
public:
	GdR();
	virtual ~GdR();
	
	virtual std::string get_output_string(llint curr_step);
	void get_settings (input_file &my_inp, input_file &sim_inp);
};

extern "C" BaseObservable<float> *make_float() { return new GdR<float>(); }
extern "C" BaseObservable<double> *make_double() { return new GdR<double>(); }

#endif /* GDR_H_ */