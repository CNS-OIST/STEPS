
#include "incgammabeta.hpp"


Gauleg18::Gauleg18 (void)
: ngau(18)
, y()
, w()

{
	y[0] = 0.0021695375159141994;
	y[1] =  0.011413521097787704;
	y[2] =  0.027972308950302116;
	y[3] =  0.051727015600492421;
	y[4] =  0.082502225484340941;
	y[5] =  0.12007019910960293;
	y[6] =  0.16415283300752470;
	y[7] =  0.21442376986779355;
	y[8] =  0.27051082840644336;
	y[9] =  0.33199876341447887;
	y[10] =  0.39843234186401943;
	y[11] =  0.46931971407375483;
	y[12] =  0.54413605556657973;
	y[13] =  0.62232745288031077;
	y[14] =  0.70331500465597174;
	y[15] =  0.78649910768313447;
	y[16] =  0.87126389619061517;
	y[17] =  0.95698180152629142;

	w[0] = 0.0055657196642445571;
	w[1] = 0.012915947284065419;
	w[2] = 0.020181515297735382;
	w[3] =  0.027298621498568734;
	w[4] =  0.034213810770299537;
	w[5] =  0.040875750923643261;
	w[6] =  0.047235083490265582;
	w[7] = 	0.053244713977759692;
	w[8] = 	0.058860144245324798;
	w[9] = 	0.064039797355015485;
	w[10] = 	0.068745323835736408;
	w[11] = 	0.072941885005653087;
	w[12] = 	0.076598410645870640;
	w[13] = 	0.079687828912071670;
	w[14] = 	0.082187266704339706;
	w[15] = 	0.084078218979661945;
	w[16] = 	0.085346685739338721;
	w[17] = 	0.085983275670394821;


};


Gamma::Gamma (void)
:Gauleg18()
, ASWITCH(100)
, EPS()
, FPMIN()
{
	EPS = numeric_limits<Doub>::epsilon();
	FPMIN = numeric_limits<Doub>::min()/EPS;
}


Doub Gamma::gammp(const Doub a, const Doub x) {
                if (x < 0.0 || a <= 0.0) throw("bad args in gammp");
                if (x == 0.0) return 0.0;
                else if ((Int)a >= ASWITCH) return gammpapprox(a,x,1);
                else if (x < a+1.0) return gser(a,x);
                else return 1.0-gcf(a,x);
        }

Doub Gamma::gammq(const Doub a, const Doub x) {
                if (x < 0.0 || a <= 0.0) throw("bad args in gammq");
                if (x == 0.0) return 1.0;
                else if ((Int)a >= ASWITCH) return gammpapprox(a,x,0);
                else if (x < a+1.0) return 1.0-gser(a,x);
                else return gcf(a,x);
        }

Doub Gamma::gser(const Doub a, const Doub x) {
                Doub sum,del,ap;
                gln=gammln(a);
                ap=a;
                del=sum=1.0/a;
                for (;;) {
                        ++ap;
                        del *= x/ap;
                        sum += del;
                        if (fabs(del) < fabs(sum)*EPS) {
                                return sum*exp(-x+a*log(x)-gln);
                        }
                }
        }

Doub Gamma::gcf(const Doub a, const Doub x) {
                Int i;
                Doub an,b,c,d,del,h;
                gln=gammln(a);
                b=x+1.0-a;
                c=1.0/FPMIN;
                d=1.0/b;
                h=d;
                for (i=1;;i++) {
                        an = -i*(i-a);
                        b += 2.0;
                        d=an*d+b;
                        if (fabs(d) < FPMIN) d=FPMIN;
                        c=b+an/c;
                        if (fabs(c) < FPMIN) c=FPMIN;
                        d=1.0/d;
                        del=d*c;
                        h *= del;
                        if (fabs(del-1.0) <= EPS) break;
                }
                return exp(-x+a*log(x)-gln)*h;
        }

Doub Gamma::gammpapprox(Doub a, Doub x, Int psig) {
                Int j;
                Doub xu,t,sum,ans;
                Doub a1 = a-1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
                gln = gammln(a);
                if (x > a1) xu = MAX(a1 + 11.5*sqrta1, x + 6.0*sqrta1);
                else xu = MAX(0.,MIN(a1 - 7.5*sqrta1, x - 5.0*sqrta1));
                sum = 0;
                for (j=0;j<ngau;j++) {
                        t = x + (xu-x)*y[j];
                        sum += w[j]*exp(-(t-a1)+a1*(log(t)-lna1));
                }
                ans = sum*(xu-x)*exp(a1*(lna1-1.)-gln);
                return (psig?(ans>0.0? 1.0-ans:-ans):(ans>=0.0? ans:1.0+ans));
        }

Doub Gamma::invgammp(Doub p, Doub a) {
        Int j;
        Doub x,err,t,u,pp,lna1,afac,a1=a-1;
        const Doub EPS=1.e-8;
        gln=gammln(a);
        if (a <= 0.) throw("a must be pos in invgammap");
        if (p >= 1.) return MAX(100.,a + 100.*sqrt(a));
        if (p <= 0.) return 0.0;
        if (a > 1.) {
                lna1=log(a1);
                afac = exp(a1*(lna1-1.)-gln);
                pp = (p < 0.5)? p : 1. - p;
                t = sqrt(-2.*log(pp));
                x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
                if (p < 0.5) x = -x;
                x = MAX(1.e-3,a*pow(1.-1./(9.*a)-x/(3.*sqrt(a)),3));
        } else {
                t = 1.0 - a*(0.253+a*0.12);
                if (p < t) x = pow(p/t,1./a);
                else x = 1.-log(1.-(p-t)/(1.-t));
        }
        for (j=0;j<12;j++) {
                if (x <= 0.0) return 0.0;
                err = gammp(a,x) - p;
                if (a > 1.) t = afac*exp(-(x-a1)+a1*(log(x)-lna1));
                else t = exp(-x+a1*log(x)-gln);
                u = err/t;
                x -= (t = u/(1.-0.5*MIN(1.,u*((a-1.)/x - 1))));
                if (x <= 0.) x = 0.5*(x + t);
                if (fabs(t) < EPS*x ) break;
        }
        return x;
}

Beta::Beta (void)
:Gauleg18()
, SWITCH(3000)
, EPS()
, FPMIN()
{
	EPS = numeric_limits<Doub>::epsilon();
	FPMIN = numeric_limits<Doub>::min()/EPS;
}


Doub Beta::betai(const Doub a, const Doub b, const Doub x) {
                Doub bt;
                if (a <= 0.0 || b <= 0.0) throw("Bad a or b in routine betai");
                if (x < 0.0 || x > 1.0) throw("Bad x in routine betai");
                if (x == 0.0 || x == 1.0) return x;
                if (a > SWITCH && b > SWITCH) return betaiapprox(a,b,x);
                bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
                if (x < (a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
                else return 1.0-bt*betacf(b,a,1.0-x)/b;
        }

Doub Beta::betacf(const Doub a, const Doub b, const Doub x) {
                Int m,m2;
                Doub aa,c,d,del,h,qab,qam,qap;
                qab=a+b;
                qap=a+1.0;
                qam=a-1.0;
                c=1.0;
                d=1.0-qab*x/qap;
                if (fabs(d) < FPMIN) d=FPMIN;
                d=1.0/d;
                h=d;
                for (m=1;m<10000;m++) {
                        m2=2*m;
                        aa=m*(b-m)*x/((qam+m2)*(a+m2));
                        d=1.0+aa*d;
                        if (fabs(d) < FPMIN) d=FPMIN;
                        c=1.0+aa/c;
                        if (fabs(c) < FPMIN) c=FPMIN;
                        d=1.0/d;
                        h *= d*c;
                        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
                        d=1.0+aa*d;
                        if (fabs(d) < FPMIN) d=FPMIN;
                        c=1.0+aa/c;
                        if (fabs(c) < FPMIN) c=FPMIN;
                        d=1.0/d;
                        del=d*c;
                        h *= del;
                        if (fabs(del-1.0) <= EPS) break;
                }
                return h;
        }

Doub Beta::betaiapprox(Doub a, Doub b, Doub x) {
                Int j;
                Doub xu,t,sum,ans;
                Doub a1 = a-1.0, b1 = b-1.0, mu = a/(a+b);
                Doub lnmu=log(mu),lnmuc=log(1.-mu);
                t = sqrt(a*b/(SQR(a+b)*(a+b+1.0)));
                if (x > a/(a+b)) {
                        if (x >= 1.0) return 1.0;
                        xu = MIN(1.,MAX(mu + 10.*t, x + 5.0*t));
                } else {
                        if (x <= 0.0) return 0.0;
                        xu = MAX(0.,MIN(mu - 10.*t, x - 5.0*t));
                }
                sum = 0;
                for (j=0;j<18;j++) {
                        t = x + (xu-x)*y[j];
                        sum += w[j]*exp(a1*(log(t)-lnmu)+b1*(log(1-t)-lnmuc));
                }
                ans = sum*(xu-x)*exp(a1*lnmu-gammln(a)+b1*lnmuc-gammln(b)+gammln(a+b));
                return ans>0.0? 1.0-ans : -ans;
        }

Doub Beta::invbetai(Doub p, Doub a, Doub b) {
                const Doub EPS = 1.e-8;
                Doub pp,t,u,err,x,al,h,w,afac,a1=a-1.,b1=b-1.;
                Int j;
                if (p <= 0.) return 0.;
                else if (p >= 1.) return 1.;
                else if (a >= 1. && b >= 1.) {
                        pp = (p < 0.5)? p : 1. - p;
                        t = sqrt(-2.*log(pp));
                        x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
                        if (p < 0.5) x = -x;
                        al = (SQR(x)-3.)/6.;
                        h = 2./(1./(2.*a-1.)+1./(2.*b-1.));
                        w = (x*sqrt(al+h)/h)-(1./(2.*b-1)-1./(2.*a-1.))*(al+5./6.-2./(3.*h));
                        x = a/(a+b*exp(2.*w));
                } else {
                        Doub lna = log(a/(a+b)), lnb = log(b/(a+b));
                        t = exp(a*lna)/a;
                        u = exp(b*lnb)/b;
                        w = t + u;
                        if (p < t/w) x = pow(a*w*p,1./a);
                        else x = 1. - pow(b*w*(1.-p),1./b);
                }
                afac = -gammln(a)-gammln(b)+gammln(a+b);
                for (j=0;j<10;j++) {
                        if (x == 0. || x == 1.) return x;
                        err = betai(a,b,x) - p;
                        t = exp(a1*log(x)+b1*log(1.-x) + afac);
                        u = err/t;
                        x -= (t = u/(1.-0.5*MIN(1.,u*(a1/x - b1/(1.-x)))));
                        if (x <= 0.) x = 0.5*(x + t);
                        if (x >= 1.) x = 0.5*(x + t + 1.);
                        if (fabs(t) < EPS*x && j > 0) break;
                }
                return x;
        }



Gammadist::Gammadist(Doub aalph, Doub bbet)
: Gamma()
, alph(aalph)
, bet(bbet)
{
    if (alph <= 0. || bet <= 0.) throw("bad alph,bet in Gammadist");
    fac = alph*log(bet)-gammln(alph);
}

Doub Gammadist::p(Doub x) {
                if (x <= 0.) throw("bad x in Gammadist");
                return exp(-bet*x+(alph-1.)*log(x)+fac);
        }

Doub Gammadist::cdf(Doub x) {
                if (x < 0.) throw("bad x in Gammadist");
                return gammp(alph,bet*x);
        }

Doub Gammadist::invcdf(Doub p) {
                if (p < 0. || p >= 1.) throw("bad p in Gammadist");
                return invgammp(p,alph)/bet;
        }


Betadist::Betadist(Doub aalph, Doub bbet)
: Beta()
, alph(aalph)
, bet(bbet)
{
    if (alph <= 0. || bet <= 0.) throw("bad alph,bet in Betadist");
    fac = gammln(alph+bet)-gammln(alph)-gammln(bet);
}

Doub Betadist::p(Doub x) {
                if (x <= 0. || x >= 1.) throw("bad x in Betadist");
                return exp((alph-1.)*log(x)+(bet-1.)*log(1.-x)+fac);
        }

Doub Betadist::cdf(Doub x) {
                if (x < 0. || x > 1.) throw("bad x in Betadist");
                return betai(alph,bet,x);
        }

Doub Betadist::invcdf(Doub p) {
                if (p < 0. || p > 1.) throw("bad p in Betadist");
                return invbetai(p,alph,bet);
        }

Studenttdist::Studenttdist(Doub nnu, Doub mmu, Doub ssig)
: Beta()
, nu(nnu)
, mu(mmu)
, sig(ssig)
{
                if (sig <= 0. || nu <= 0.) throw("bad sig,nu in Studentdist");
                np = 0.5*(nu + 1.);
                fac = gammln(np)-gammln(0.5*nu);
        }

Doub Studenttdist::p(Doub t) {
                return exp(-np*log(1.+SQR((t-mu)/sig)/nu)+fac)
                        /(sqrt(3.14159265358979324*nu)*sig);
        }

Doub Studenttdist::cdf(Doub t) {
                Doub p = 0.5*betai(0.5*nu, 0.5, nu/(nu+SQR((t-mu)/sig)));
                if (t >= mu) return 1. - p;
                else return p;
        }

Doub Studenttdist::invcdf(Doub p) {
                if (p <= 0. || p >= 1.) throw("bad p in Studentdist");
                Doub x = invbetai(2.*MIN(p,1.-p), 0.5*nu, 0.5);
                x = sig*sqrt(nu*(1.-x)/x);
                return (p >= 0.5? mu+x : mu-x);
        }

Doub Studenttdist::aa(Doub t) {
                if (t < 0.) throw("bad t in Studentdist");
                return 1.-betai(0.5*nu, 0.5, nu/(nu+SQR(t)));
        }

Doub Studenttdist::invaa(Doub p) {
                if (p < 0. || p >= 1.) throw("bad p in Studentdist");
                Doub x = invbetai(1.-p, 0.5*nu, 0.5);
                return sqrt(nu*(1.-x)/x);
        }

Poissondist::Poissondist(Doub llam)
: Gamma()
, lam(llam)
{
                if (lam <= 0.) throw("bad lam in Poissondist");
        }

Doub Poissondist::p(Int n) {
                if (n < 0) throw("bad n in Poissondist");
                return exp(-lam + n*log(lam) - gammln(n+1.));
        }

Doub Poissondist::cdf(Int n) {
                if (n < 0) throw("bad n in Poissondist");
                if (n == 0) return 0.;
                return gammq((Doub)n,lam);
        }

Int Poissondist::invcdf(Doub p) {
                Int n,nl,nu,inc=1;
                if (p <= 0. || p >= 1.) throw("bad p in Poissondist");
                if (p < exp(-lam)) return 0;
                n = (Int)MAX(sqrt(lam),5.);
                if (p < cdf(n)) {
                        do {
                                n = MAX(n-inc,0);
                                inc *= 2;
                        } while (p < cdf(n));
                        nl = n; nu = n + inc/2;
                } else {
                        do {
                                n += inc;
                                inc *= 2;
                        } while (p > cdf(n));
                        nu = n; nl = n - inc/2;
                }
                while (nu-nl>1) {
                        n = (nl+nu)/2;
                        if (p < cdf(n)) nu = n;
                        else nl = n;
                }
                return nl;
        }

Binomialdist::Binomialdist(Int nn, Doub ppe)
: Beta()
, n(nn)
, pe(ppe) {
                if (n <= 0 || pe <= 0. || pe >= 1.) throw("bad args in Binomialdist");
                fac = gammln(n+1.);
        }

Doub Binomialdist::p(Int k) {
                if (k < 0) throw("bad k in Binomialdist");
                if (k > n) return 0.;
                return exp(k*log(pe)+(n-k)*log(1.-pe)
                        +fac-gammln(k+1.)-gammln(n-k+1.));
        }

Doub Binomialdist::cdf(Int k) {
                if (k < 0) throw("bad k in Binomialdist");
                if (k == 0) return 0.;
                if (k > n) return 1.;
                return 1. - betai((Doub)k,n-k+1.,pe);
        }

Int Binomialdist::invcdf(Doub p) {
                Int k,kl,ku,inc=1;
                if (p <= 0. || p >= 1.) throw("bad p in Binomialdist");
                k = MAX(0,MIN(n,(Int)(n*pe)));
                if (p < cdf(k)) {
                        do {
                                k = MAX(k-inc,0);
                                inc *= 2;
                        } while (p < cdf(k));
                        kl = k; ku = k + inc/2;
                } else {
                        do {
                                k = MIN(k+inc,n+1);
                                inc *= 2;
                        } while (p > cdf(k));
                        ku = k; kl = k - inc/2;
                }
                while (ku-kl>1) {
                        k = (kl+ku)/2;
                        if (p < cdf(k)) ku = k;
                        else kl = k;
                }
                return kl;
        }

Chisqdist::Chisqdist(Doub nnu)
: Gamma()
, nu(nnu)
{
                if (nu <= 0.) throw("bad nu in Chisqdist");
                fac = 0.693147180559945309*(0.5*nu)+gammln(0.5*nu);
        }

Doub Chisqdist::p(Doub x2) {
                if (x2 <= 0.) throw("bad x2 in Chisqdist");
                return exp(-0.5*(x2-(nu-2.)*log(x2))-fac);
        }

Doub Chisqdist::cdf(Doub x2) {
                if (x2 < 0.) throw("bad x2 in Chisqdist");
                return gammp(0.5*nu,0.5*x2);
        }

Doub Chisqdist::invcdf(Doub p) {
                if (p < 0. || p >= 1.) throw("bad p in Chisqdist");
                return 2.*invgammp(p,0.5*nu);
        }

Fdist::Fdist(Doub nnu1, Doub nnu2)
: Beta()
, nu1(nnu1)
, nu2(nnu2)
{
                if (nu1 <= 0. || nu2 <= 0.) throw("bad nu1,nu2 in Fdist");
                fac = 0.5*(nu1*log(nu1)+nu2*log(nu2))+gammln(0.5*(nu1+nu2))
                        -gammln(0.5*nu1)-gammln(0.5*nu2);
        }

Doub Fdist::p(Doub f) {
                if (f <= 0.) throw("bad f in Fdist");
                return exp((0.5*nu1-1.)*log(f)-0.5*(nu1+nu2)*log(nu2+nu1*f)+fac);
        }

Doub Fdist::cdf(Doub f) {
                if (f < 0.) throw("bad f in Fdist");
                return betai(0.5*nu1,0.5*nu2,nu1*f/(nu2+nu1*f));
        }

Doub Fdist::invcdf(Doub p) {
                if (p <= 0. || p >= 1.) throw("bad p in Fdist");
                Doub x = invbetai(p,0.5*nu1,0.5*nu2);
                return nu2*x/(nu1*(1.-x));
        }



