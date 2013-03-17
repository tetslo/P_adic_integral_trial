/*
 *  P_adic_number.cpp
 *  P_adic_integral
 *
 *  Created by Olga Stetyukhina on 3/11/13.
 *  Copyright 2013 RusBigHedgehog. All rights reserved.
 *
 */

#include "P_adic_number.h"
#include <math.h>

using namespace std;

// initializes canonical expansion of p-adic number by nulls

P_adic_number :: P_adic_number() {
	coef.resize(num_size);
}

// initializes canonical expansion of p-adic number by array

P_adic_number :: P_adic_number(const int* arr) {
	coef.resize(num_size);
	
	vector<int>::iterator it = coef.begin();
	int i = 0;
	while (it != coef.end()) {
		*it = arr[i];
		it++;
		i++;
	}
}

// gets the next p-adic number in subsetting of integration domain

bool P_adic_number :: next_number() {
	vector<int>::reverse_iterator it;
	it = coef.rbegin();
	while (it != coef.rend()) {
		if (*it < p - 1) {
			*it += 1;
			return true;
		}
		*it = 0;
		it++;
	}
	
	return false;
}

// outputs the number in console (can be used for verification)

void P_adic_number :: print_number() {
	vector<int>::iterator it;
	int i = 0;
	for (it = coef.begin(); it != coef.end(); it++) {
		if (i == gamma_max) {
			cout << ".";
		}
		cout << *it;
		i++;
	}
	cout << endl;
}

// returns the logarithm of norm of difference (gamma in B_gamma)

int log_norm_dif(const P_adic_number a, const P_adic_number b) {
	int gamma = gamma_max;
	vector<int>::const_iterator a_it;
	vector<int>::const_iterator b_it;
	a_it = a.coef.begin();
	b_it = b.coef.begin();
	
	while (a_it != a.coef.end()) {
		if (*a_it != *b_it) {
			return gamma;
		}
		gamma--;
		a_it++;
		b_it++;
	}
	
	return gamma;
}

// returns indicator-function value

complex_d indicator(const P_adic_number x, const P_adic_number c, const int gamma) {
	if (log_norm_dif(x, c) <= gamma) {
		return 1;
	}
	return 0;
}

// returns wavelet-function value

complex_d wavelet(const P_adic_number x, const P_adic_number c, const int gamma) {
	int g_index = gamma_max - gamma;
	P_adic_number n;
	for (int k = 0; k < p; k++) {
		int i;
		for (i = 0; i < g_index; i++) {
			n.coef[i] = c.coef[i];
		}
		n.coef[i] = k;
		if (log_norm_dif(x, n) <= gamma - 1) {
			double s_arg = 2 * M_PI * k / p;
			complex_d s = polar(1.0, s_arg);
			return s;
		}
	}
	return 0;
}

// returns integral value

complex_d integral() {
	int center_arr[] = {0, 2, 0, 0, 2, 0, 0};
	
	P_adic_number center (center_arr);		// for indicator or wavelet
	P_adic_number x;						// variable
	
	int gamma_ind = 3;						// for indicator or wavelet 
			
	complex_d s = 0;
	double ball_measure = pow((double) p, (double) (gamma_min - 1));
	do {		
//		s += ball_measure * indicator(x, center, gamma_ind);
		s += ball_measure * wavelet(x, center, gamma_ind);
	} while (x.next_number());
	
	if (abs(s) < 0.00000001) {
		return 0;
	}
	
	return s;
}
