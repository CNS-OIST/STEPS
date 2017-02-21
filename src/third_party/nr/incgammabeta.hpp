#ifndef _IGB_H_
#define _IGB_H_ 1

#include "nr3.hpp"
#include "gamma.hpp"


class Gauleg18 {
	public:

		Gauleg18 (void);

        Int ngau;
        Doub y[18];
        Doub w[18];
};

class Gamma : Gauleg18 {
	public:

		Gamma(void);

        Int ASWITCH;
        Doub EPS;
        Doub FPMIN;
        Doub gln;

        Doub gammp(const Doub a, const Doub x);

        Doub gammq(const Doub a, const Doub x);

        Doub gser(const Doub a, const Doub x);

        Doub gcf(const Doub a, const Doub x);

        Doub gammpapprox(Doub a, Doub x, Int psig);

        Doub invgammp(Doub p, Doub a);

};

class Beta : Gauleg18 {
	public:

		Beta(void);

        Int SWITCH;
        Doub EPS;
        Doub FPMIN;

        Doub betai(const Doub a, const Doub b, const Doub x);

        Doub betacf(const Doub a, const Doub b, const Doub x);

        Doub betaiapprox(Doub a, Doub b, Doub x);

        Doub invbetai(Doub p, Doub a, Doub b) ;

};

class Gammadist : Gamma {
	public:

        Gammadist(Doub aalph, Doub bbet = 1.);

		Doub alph;
		Doub bet;
		Doub fac;

        Doub p(Doub x);

        Doub cdf(Doub x);

        Doub invcdf(Doub p);

};

class Betadist : Beta {

	public:

		Betadist(Doub aalph, Doub bbet);
        Doub alph;
        Doub bet;
        Doub fac;
        Doub p(Doub x) ;
        Doub cdf(Doub x);
        Doub invcdf(Doub p);
};

class Studenttdist : Beta {

	public:

		Studenttdist(Doub nnu, Doub mmu = 0., Doub ssig = 1.);

        Doub nu;
        Doub mu;
        Doub sig;
        Doub np;
        Doub fac;

        Doub p(Doub t);

        Doub cdf(Doub t);

        Doub invcdf(Doub p);

        Doub aa(Doub t);

        Doub invaa(Doub p);

};

class Poissondist : Gamma {
	public:
        Poissondist(Doub llam);

        Doub lam;

        Doub p(Int n);

        Doub cdf(Int n);

        Int invcdf(Doub p);
};

class Binomialdist : Beta {
	public:
        Binomialdist(Int nn, Doub ppe);

        Int n;
        Doub pe;
        Doub fac;

        Doub p(Int k);
        Doub cdf(Int k);
        Int invcdf(Doub p);
};

class Chisqdist : Gamma {

public:
	Chisqdist(Doub nnu);

        Doub nu;
        Doub fac;

        Doub p(Doub x2);

        Doub cdf(Doub x2);

        Doub invcdf(Doub p);

};
class Fdist : Beta {

public:
    Fdist(Doub nnu1, Doub nnu2);

		Doub nu1;
		Doub nu2;
        Doub fac;

        Doub p(Doub f);

        Doub cdf(Doub f) ;

        Doub invcdf(Doub p);

};


#endif /* _IGB_H_ */

