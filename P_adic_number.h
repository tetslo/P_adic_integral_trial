/*
 *  P_adic_number.h
 *  P_adic_integral
 *
 *  Created by Olga Stetyukhina on 3/11/13.
 *  Copyright 2013 RusBigHedgehog. All rights reserved.
 *
 */

#include <vector>
#include <complex>

#include <iostream>

typedef std::complex<double> complex_d;

const int p = 3;

const int gamma_min = -3;
const int gamma_max = 3;
const int num_size = gamma_max - gamma_min + 1;

// canonical expansion of p_adic number
struct P_adic_number {
	std::vector<int> coef;
	
	P_adic_number();
	P_adic_number(const int*);
	
	bool next_number();
	void print_number();
};


int log_norm_dif(const P_adic_number, const P_adic_number);

complex_d indicator(const P_adic_number, const P_adic_number, const int);

complex_d wavelet(const P_adic_number, const P_adic_number, const int);

complex_d integral();


