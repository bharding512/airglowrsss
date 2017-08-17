// SunAndMoon.cpp: implementation of the CSunAndMoon class.
//
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "SunAndMoon.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CSunAndMoon::CSunAndMoon()
{
	m_MoonAngle = STDHZN;
	m_SunAngle = TWILIGHT;
}

CSunAndMoon::~CSunAndMoon()
{

}

/* given the mjd, find the geocentric ecliptic longitude, lam, and latitude,
 * bet, and horizontal parallax, hp for the moon.
 * N.B. series for long and lat are good to about 10 and 3 arcseconds. however,
 *   math errors cause up to 100 and 30 arcseconds error, even if use double.
 *   why?? suspect highly sensitive nature of difference used to get m1..6.
 * N.B. still need to correct for nutation. then for topocentric location
 *   further correct for parallax and refraction.
 */
void CSunAndMoon::moon (double mjd, double *lam, double *bet, double *hp)
{
	double t, t2;
	double ld;
	double ms;
	double md;
	double de;
	double f;
	double n;
	double a, sa, sn, b, sb, c, sc, e, e2, l, g, w1, w2;
	double m1, m2, m3, m4, m5, m6;

	t = mjd/36525.;
	t2 = t*t;

	m1 = mjd/27.32158213;
	m1 = 360.0*(m1-(long)m1);
	m2 = mjd/365.2596407;
	m2 = 360.0*(m2-(long)m2);
	m3 = mjd/27.55455094;
	m3 = 360.0*(m3-(long)m3);
	m4 = mjd/29.53058868;
	m4 = 360.0*(m4-(long)m4);
	m5 = mjd/27.21222039;
	m5 = 360.0*(m5-(long)m5);
	m6 = mjd/6798.363307;
	m6 = 360.0*(m6-(long)m6);

	ld = 270.434164+m1-(.001133-.0000019*t)*t2;
	ms = 358.475833+m2-(.00015+.0000033*t)*t2;
	md = 296.104608+m3+(.009192+.0000144*t)*t2;
	de = 350.737486+m4-(.001436-.0000019*t)*t2;
	f = 11.250889+m5-(.003211+.0000003*t)*t2;
	n = 259.183275-m6+(.002078+.000022*t)*t2;

	a = degrad(51.2+20.2*t);
	sa = sin(a);
	sn = sin(degrad(n));
	b = 346.56+(132.87-.0091731*t)*t;
	sb = .003964*sin(degrad(b));
	c = degrad(n+275.05-2.3*t);
	sc = sin(c);
	ld = ld+.000233*sa+sb+.001964*sn;
	ms = ms-.001778*sa;
	md = md+.000817*sa+sb+.002541*sn;
	f = f+sb-.024691*sn-.004328*sc;
	de = de+.002011*sa+sb+.001964*sn;
	e = 1-(.002495+7.52e-06*t)*t;
	e2 = e*e;

	ld = degrad(ld);
	ms = degrad(ms);
	n = degrad(n);
	de = degrad(de);
	f = degrad(f);
	md = degrad(md);

	l = 6.28875*sin(md)+1.27402*sin(2*de-md)+.658309*sin(2*de)+
	    .213616*sin(2*md)-e*.185596*sin(ms)-.114336*sin(2*f)+
	    .058793*sin(2*(de-md))+.057212*e*sin(2*de-ms-md)+
	    .05332*sin(2*de+md)+.045874*e*sin(2*de-ms)+.041024*e*sin(md-ms);
	l = l-.034718*sin(de)-e*.030465*sin(ms+md)+.015326*sin(2*(de-f))-
	    .012528*sin(2*f+md)-.01098*sin(2*f-md)+.010674*sin(4*de-md)+
	    .010034*sin(3*md)+.008548*sin(4*de-2*md)-e*.00791*sin(ms-md+2*de)-
	    e*.006783*sin(2*de+ms);
	l = l+.005162*sin(md-de)+e*.005*sin(ms+de)+.003862*sin(4*de)+
	    e*.004049*sin(md-ms+2*de)+.003996*sin(2*(md+de))+
	    .003665*sin(2*de-3*md)+e*.002695*sin(2*md-ms)+
	    .002602*sin(md-2*(f+de))+e*.002396*sin(2*(de-md)-ms)-
	    .002349*sin(md+de);
	l = l+e2*.002249*sin(2*(de-ms))-e*.002125*sin(2*md+ms)-
	    e2*.002079*sin(2*ms)+e2*.002059*sin(2*(de-ms)-md)-
	    .001773*sin(md+2*(de-f))-.001595*sin(2*(f+de))+
	    e*.00122*sin(4*de-ms-md)-.00111*sin(2*(md+f))+.000892*sin(md-3*de);
	l = l-e*.000811*sin(ms+md+2*de)+e*.000761*sin(4*de-ms-2*md)+
	     e2*.000704*sin(md-2*(ms+de))+e*.000693*sin(ms-2*(md-de))+
	     e*.000598*sin(2*(de-f)-ms)+.00055*sin(md+4*de)+.000538*sin(4*md)+
	     e*.000521*sin(4*de-ms)+.000486*sin(2*md-de);
	l = l+e2*.000717*sin(md-2*ms);
	*lam = ld+degrad(l);
	range (lam, 2*PI);

	g = 5.12819*sin(f)+.280606*sin(md+f)+.277693*sin(md-f)+
	    .173238*sin(2*de-f)+.055413*sin(2*de+f-md)+.046272*sin(2*de-f-md)+
	    .032573*sin(2*de+f)+.017198*sin(2*md+f)+.009267*sin(2*de+md-f)+
	    .008823*sin(2*md-f)+e*.008247*sin(2*de-ms-f);
	g = g+.004323*sin(2*(de-md)-f)+.0042*sin(2*de+f+md)+
	    e*.003372*sin(f-ms-2*de)+e*.002472*sin(2*de+f-ms-md)+
	    e*.002222*sin(2*de+f-ms)+e*.002072*sin(2*de-f-ms-md)+
	    e*.001877*sin(f-ms+md)+.001828*sin(4*de-f-md)-e*.001803*sin(f+ms)-
	    .00175*sin(3*f);
	g = g+e*.00157*sin(md-ms-f)-.001487*sin(f+de)-e*.001481*sin(f+ms+md)+
	     e*.001417*sin(f-ms-md)+e*.00135*sin(f-ms)+.00133*sin(f-de)+
	     .001106*sin(f+3*md)+.00102*sin(4*de-f)+.000833*sin(f+4*de-md)+
	     .000781*sin(md-3*f)+.00067*sin(f+4*de-2*md);
	g = g+.000606*sin(2*de-3*f)+.000597*sin(2*(de+md)-f)+
	    e*.000492*sin(2*de+md-ms-f)+.00045*sin(2*(md-de)-f)+
	    .000439*sin(3*md-f)+.000423*sin(f+2*(de+md))+
	    .000422*sin(2*de-f-3*md)-e*.000367*sin(ms+f+2*de-md)-
	    e*.000353*sin(ms+f+2*de)+.000331*sin(f+4*de);
	g = g+e*.000317*sin(2*de+f-ms+md)+e2*.000306*sin(2*(de-ms)-f)-
	    .000283*sin(md+3*f);
	w1 = .0004664*cos(n);
	w2 = .0000754*cos(c);
	*bet = degrad(g)*(1-w1-w2);

	*hp = .950724+.051818*cos(md)+.009531*cos(2*de-md)+.007843*cos(2*de)+
	      .002824*cos(2*md)+.000857*cos(2*de+md)+e*.000533*cos(2*de-ms)+
	      e*.000401*cos(2*de-md-ms)+e*.00032*cos(md-ms)-.000271*cos(de)-
	      e*.000264*cos(ms+md)-.000198*cos(2*f-md);
	*hp = *hp+.000173*cos(3*md)+.000167*cos(4*de-md)-e*.000111*cos(ms)+
	     .000103*cos(4*de-2*md)-.000084*cos(2*md-2*de)-
	     e*.000083*cos(2*de+ms)+.000079*cos(2*de+2*md)+.000072*cos(4*de)+
	     e*.000064*cos(2*de-ms+md)-e*.000063*cos(2*de+ms-md)+
	     e*.000041*cos(ms+de);
	*hp = *hp+e*.000035*cos(2*md-ms)-.000033*cos(3*md-2*de)-
	     .00003*cos(md+de)-.000029*cos(2*(f-de))-e*.000029*cos(2*md+ms)+
	     e2*.000026*cos(2*(de-ms))-.000023*cos(2*(f-de)+md)+
	     e*.000019*cos(4*de-ms-md);
	*hp = degrad(*hp);
}

void CSunAndMoon::riset (double ra, double dec, double lat, double dis, 
			double *lstr, double *lsts, double *azr, double *azs, int *status)
{
#define	EPS	(1e-6)	/* math rounding fudge - always the way, eh? */
	double d;	/* angle from pole */
	double h;	/* hour angle */
	double crho;	/* cos hour-angle complement */
	int shemi;	/* flag for southern hemisphere reflection */

	d = PI/2 - dec;

	/* reflect if in southern hemisphere.
	 * (then reflect azimuth back after computation.)
	 */
	if (shemi = lat < 0) {
	    lat = -lat;
	    d = PI - d;
	}

	/* do the easy ones (and avoid CSunAndMoon::violated assumptions) if d arc never
	 * meets horizon. 
	 */
	if (d <= lat + dis + EPS) {
	    *status = -1; /* never sets */
	    return;
	}
	if (d >= PI - lat + dis - EPS) {
	    *status = 1; /* never rises */
	    return;
	}

	/* find rising azimuth and cosine of hour-angle complement */
	if (lat > EPS) {
	    double d2, d1; /* polr arc to ideal hzn, and corrctn for apparent */
	    double z2, z1; /* azimuth to ideal horizon, and " */
	    double a;	   /* intermediate temp */
	    double sdis, slat, clat, cz2, cd2;	/* trig temps */
	    sdis = sin(dis);
	    slat = sin(lat);
	    a = sdis*sdis + slat*slat + 2*cos(d)*sdis*slat;
	    if (a <= 0) {
		*status = 2; /* can't happen - hah! */
		return;
	    }

	    d1 = asin (sin(d) * sdis / sqrt(a));
	    d2 = d - d1;
	    cd2 = cos(d2);
	    clat = cos(lat);
	    cz2 = cd2/clat;
	    z2 = acos (cz2);
	    z1 = acos (cos(d1)/cos(dis));
	    if (dis < 0)
		z1 = -z1;
	    *azr = z1 + z2;
	    range (azr, PI);
	    crho = (cz2 - cd2*clat)/(sin(d2)*slat);
	} else {
	    *azr = acos (cos(d)/cos(dis));
	    crho = sin(dis)/sin(d);
	}

	if (shemi)
	    *azr = PI - *azr;
        *azs = 2*PI - *azr;
	
	/* find hour angle */
	h = PI - acos (crho);
        *lstr = radhr(ra-h);
	*lsts = radhr(ra+h);
	range (lstr, 24.0);
	range (lsts, 24.0);

	*status = 0;
}


/* insure 0 <= *v < r.
 */
void CSunAndMoon::range (double *v, double r)
{
	while (*v <  0) *v += r;
	while (*v >= r) *v -= r;
}

double CSunAndMoon::mjd_day(double jd)
{
	return (floor(jd-0.5)+0.5);
}

double CSunAndMoon::mjd_hr(double jd)
{
	return ((jd-mjd_day(jd))*24.0);
}

void CSunAndMoon::now_lst (Now *np, double *lst)
{
	utc_gst (mjd_day(np->n_mjd), mjd_hr(np->n_mjd), lst);
	*lst += radhr(np->n_lng);
	range (lst, 24.0);
}

/* given a modified julian date, mjd, and a universally coordinated time, utc,
 * return greenwich mean siderial time, *gst.
 */
void CSunAndMoon::utc_gst (double dmjd, double utc, double *gst)
{
	static double lastmjd = -10000;
	static double t0;

	if (dmjd != lastmjd) {
	    t0 = tnaught (dmjd);
	    lastmjd = dmjd;
	}
	*gst = (1.0/SIDRATE)*utc + t0;
	range (gst, 24.0);
}

double CSunAndMoon::tnaught(double ddmjd)
/* julian days since 1900 jan 0.5 */
{
	double dmjd;
	int m, y;
	double d;
	double t, t0;

	mjd_cal (ddmjd, &m, &d, &y);
	cal_mjd (1, 0., y, &dmjd);
	t = dmjd/36525;
	t0 = 6.57098e-2 * (ddmjd - dmjd) - 
	     (24 - (6.6460656 + (5.1262e-2 + (t * 2.581e-5))*t) -
		   (2400 * (t - (((double)y - 1900)/100))));
	return (t0);
}


/* given the modified Julian date (number of days elapsed since 1900 jan 0.5,),
 * mjd, return the calendar date in months, *mn, days, *dy, and years, *yr.
 */
void CSunAndMoon::mjd_cal (double ddmjd, int *mn, double *dy, int *yr)
{
	double d, f;
	double i, a, b, ce, g;

	d = ddmjd + 0.5;
	i = floor(d);
	f = d-i;
	if (f == 1) {
	    f = 0;
	    i += 1;
	}

	if (i > -115860.0) {
	    a = floor((i/36524.25)+.9983573)+14;
	    i += 1 + a - floor(a/4.0);
	}

	b = floor((i/365.25)+.802601);
	ce = i - floor((365.25*b)+.750001)+416;
	g = floor(ce/30.6001);
	*mn = (int) g - 1;
	*dy = ce - floor(30.6001*g)+f;
	*yr = (int) b + 1899;

	if (g > 13.5)
	    *mn = (int) g - 13;
	if (*mn < 2.5)
	    *yr = (int) b + 1900;
	if (*yr < 1)
	    *yr -= 1;
}

/* given a date in months, mn, days, dy, years, yr,
 * return the modified Julian date (number of days elapsed since 1900 jan 0.5),
 * *mjd.
 */
void CSunAndMoon::cal_mjd (int mn, double dy, int yr, double *ddmjd)
{
	int b, d, m, y;
	long c;

	m = mn;
	y = (yr < 0) ? yr + 1 : yr;
	if (mn < 3) {
	    m += 12;
	    y -= 1;
	}

	if (yr < 1582 || yr == 1582 && (mn < 10 || mn == 10 && dy < 15)) 
	    b = 0;
	else {
	    int a;
	    a = y/100;
	    b = 2 - a + a/4;
	}

	if (y < 0)
	    c = (long)((365.25*y) - 0.75) - 694025L;
	else
	    c = (long)(365.25*y) - 694025L;

	d = (int) (30.6001*(m+1));

	*ddmjd = b + c + d + dy - 0.5;
}

/* find body p's circumstances now.
 * to save some time the caller may specify a desired accuracy, in arc seconds.
 * if, based on its mean motion, it would not have moved this much since the
 * last time we were called we only recompute altitude and azimuth and avoid
 * recomputing the planet's heliocentric position. use 0.0 for best possible.
 * we always recompute the user-defined objects' position regardless.
 * return 0 if only alt/az changes, else 1 if all other stuff updated too.
 * TODO: beware of fast retrograde motions.
 */
void CSunAndMoon::body_cir (int p, double as, Now *np, Sky *sp)
{
	typedef struct {
	    double l_dpas;	/* mean days per arc second */
	    Now l_now;		/* when l_sky was found */
	    double l_ra, l_dec;	/* the eod, ie, unprecessed, ra/dec values */
	    Sky l_sky;
	} Last;
	/* must be in same order as the astro.h object #define's */
	static Last last[8] = {
	    {.000068, {NOMJD}},
	    {.00017, {NOMJD}},
	    {.00053, {NOMJD}},
	    {.0034, {NOMJD}},
	    {.0092, {NOMJD}},
	    {.027, {NOMJD}},
	    {.046, {NOMJD}},
	    {.069, {NOMJD}}
	};
	//Last objxlast, objylast;
	double lst, alt, az;
	double ehp, ha, dec;	/* ehp: angular dia of earth from body */
	Last *lp;
    double lpd0, psi0;	/* heliocentric ecliptic long and lat */
    
	    double rp0;		/* dist from sun */
	    double rho0;	/* dist from earth */
	    double lam, bet;	/* geocentric ecliptic long and lat */
	    double dia, mag;	/* angular diameter at 1 AU and magnitude */
	    double lsn, rsn;	/* true geoc lng of sun, dist from sn to earth*/
	    double el;	/* elongation */
	    double f;   /* phase from earth */

	switch (p) {
	case SUN: sun_cir (as, np, sp);
				return;
	case MOON: moon_cir (as, np, sp);
				return;
	default: lp = last + p; break;
	}

	
	    lp->l_now = *np;
	    sunpos (np->n_mjd, &lsn, &rsn);

		double deps, dpsi;
		double a;
		plans(np->n_mjd, p, &lpd0, &psi0, &rp0, &rho0, &lam, &bet, &dia,&mag);
		nutation (np->n_mjd, &deps, &dpsi);	/* correct for nutation */
		lam += dpsi;
		a = lsn-lam;			/* and 20.4" aberation */
		lam -= degrad(20.4/3600)*cos(a)/cos(bet);
		bet -= degrad(20.4/3600)*sin(a)*sin(bet);

	    ecl_eq (np->n_mjd, bet, lam, &lp->l_ra, &lp->l_dec);

	    sp->s_ra = lp->l_ra;

	    sp->s_dec = lp->l_dec;
	    if (np->n_epoch != EOD)
		precess (np->n_mjd, np->n_epoch, &sp->s_ra, &sp->s_dec);
	    sp->s_edist = rho0;
	    sp->s_sdist = rp0;
	    elongation (lam, bet, lsn, &el);
	    el = raddeg(el);
	    sp->s_elong = el;
	    f = (rp0 > 0.0)
		? 0.25 * (((rp0+rho0)*(rp0+rho0) - rsn*rsn)/(rp0*rho0)) : 0.0;
	    sp->s_phase = f*100.0; /* percent */

		sp->s_size = dia/rho0;
		sp->s_mag = mag + 5.0*log(rp0*rho0/sqrt(f))/log(10.0);

	    sp->s_hlong = lpd0;
	    sp->s_hlat = psi0;

	/* alt, az; correct for parallax and refraction; use eod ra/dec */
	now_lst (np, &lst);
	ha = hrrad(lst) - lp->l_ra;
	if (sp->s_edist > 0.0) {
	    ehp = (2.0*6378.0/146.0e6) / sp->s_edist;
	    ta_par (ha, lp->l_dec, np->n_lat, np->n_height, ehp, &ha, &dec);
	} else
	    dec = lp->l_dec;
	hadec_aa (np->n_lat, ha, dec, &alt, &az);
	refract (np->n_pressure, np->n_temp, alt, &alt);
	sp->s_alt = alt;
	sp->s_az = az;
	lp->l_sky = *sp;
}

/* corrects ra and dec, both in radians, for precession from epoch 1 to epoch 2.
 * the epochs are given by their modified JDs, mjd1 and mjd2, respectively.
 * N.B. ra and dec are modifed IN PLACE.
 * TODO: find a better algorithm; this one is not even symmetric.
 */
void CSunAndMoon::precess (double mjd1, double mjd2, double *ra, double *dec)
/*double mjd1, mjd2;	/* initial and final epoch modified JDs *
double *ra, *dec;	/* ra/dec for mjd1 in, for mjd2 out */
{
	static double lastmjd1 = -10000, lastmjd2 = -10000;
	static double m, n, nyrs;
	double dra, ddec;	/* ra and dec corrections */

	if (mjd1 != lastmjd1 || mjd2 != lastmjd2) {
	    double t1, t2; /* Julian centuries of 36525 days since Jan .5 1900*/
	    double m1, n1; /* "constants" for t1 */
	    double m2, n2; /* "constants" for t2 */
	    t1 = mjd1/36525.;
	    t2 = mjd2/36525.;
	    m1 = 3.07234+(1.86e-3*t1);
	    n1 = 20.0468-(8.5e-3*t1);
	    m2 = 3.07234+(1.86e-3*t2);
	    n2 = 20.0468-(8.5e-3*t2);
	    m = (m1+m2)/2;	/* average m for the two epochs */
	    n = (n1+n2)/2;	/* average n for the two epochs */
	    nyrs = (mjd2-mjd1)/365.2425;
	    lastmjd1 = mjd1;
	    lastmjd2 = mjd2;
	}

	dra = (m+(n*sin(*ra)*tan(*dec)/15))*7.272205e-5*nyrs;
	ddec = n*cos(*ra)*4.848137e-6*nyrs;
	*ra += dra;
	*dec += ddec;
	/* added by ECD */
	if (*dec > PI/2) {
	    *dec = PI - *dec;
	    *ra += PI;
	} else if (*dec < -PI/2) {
	    *dec = -PI - *dec;
	    *ra += PI;
	}
	range (ra, 2*PI);
}
/* given true ha and dec, tha and tdec, the geographical latitude, phi, the
 * height above sea-level (as a fraction of the earths radius, 6378.16km),
 * ht, and the equatorial horizontal parallax, ehp, find the apparent
 * ha and dec, aha and adec allowing for parallax.
 * all angles in radians. ehp is the angle subtended at the body by the
 * earth's equator.
 */
void CSunAndMoon::ta_par (double tha, double tdec, double phi, double ht, double ehp, double *aha, double *adec)
{
	static double last_phi, last_ht, rsp, rcp;
	double rp;	/* distance to object in Earth radii */
	double ctha;
	double stdec, ctdec;
	double tdtha, dtha;
	double caha;

	/* avoid CSunAndMoon::calcs involving the same phi and ht */
	if (phi != last_phi || ht != last_ht) {
	    double cphi, sphi, u;
	    cphi = cos(phi);
	    sphi = sin(phi);
	    u = atan(9.96647e-1*sphi/cphi);
	    rsp = (9.96647e-1*sin(u))+(ht*sphi);
	    rcp = cos(u)+(ht*cphi);
	    last_phi  =  phi;
	    last_ht  =  ht;
	}

        rp = 1/sin(ehp);

        ctha = cos(tha);
	stdec = sin(tdec);
	ctdec = cos(tdec);
        tdtha = (rcp*sin(tha))/((rp*ctdec)-(rcp*ctha));
        dtha = atan(tdtha);
	*aha = tha+dtha;
	caha = cos(*aha);
	range (aha, 2*PI);
        *adec = atan(caha*(rp*stdec-rsp)/(rp*ctdec*ctha-rcp));
}


/* correct the true altitude, ta, for refraction to the apparent altitude, aa,
 * each in radians, given the local atmospheric pressure, pr, in mbars, and
 * the temperature, tr, in degrees C.
 */
void CSunAndMoon::refract (double pr, double tr, double ta, double *aa)
{
	double r;	/* refraction correction*/

        if (ta >= degrad(15.)) {
	    /* model for altitudes at least 15 degrees above horizon */
            r = 7.888888e-5*pr/((273+tr)*tan(ta));
	} else if (ta > degrad(-5.)) {
	    /* hairier model for altitudes at least -5 and below 15 degrees */
	    double a, b, tadeg = raddeg(ta);
	    a = ((2e-5*tadeg+1.96e-2)*tadeg+1.594e-1)*pr;
	    b = (273+tr)*((8.45e-2*tadeg+5.05e-1)*tadeg+1);
	    r = degrad(a/b);
	} else {
	    /* do nothing if more than 5 degrees below horizon.
	     */
	    r = 0;
	}

	*aa  =  ta + r;
}

/* correct the apparent altitude, aa, for refraction to the true altitude, ta,
 * each in radians, given the local atmospheric pressure, pr, in mbars, and
 * the temperature, tr, in degrees C.
 */
void CSunAndMoon::unrefract (double pr, double tr, double aa, double *ta)
{
	double err;
	double appar;
	double truef;

	/* iterative solution: search for the true that refracts to the
	 *   given apparent.
	 * since refract() is discontinuous at -5 degrees, there is a range
	 *   of apparent altitudes between about -4.5 and -5 degrees that are
	 *   not invertable (the graph of ap vs. true has a vertical step at
	 *   true = -5). thus, the iteration just oscillates if it gets into
	 *   this region. if this happens the iteration is forced to abort.
	 *   of course, this makes unrefract() discontinuous too.
	 */
	truef = aa;
	do {
	    refract (pr, tr, truef, &appar);
	    err = appar - aa;
	    truef -= err;
	} while (fabs(err) >= 1e-6 && truef > degrad(-5));

	*ta = truef;
}

/* find sun's circumstances now.
 * as is the desired accuracy, in arc seconds; use 0.0 for best possible.
 * return 0 if only alt/az changes, else 1 if all other stuff updated too.
 */
void CSunAndMoon::sun_cir (double as, Now *np, Sky *sp)
{
	static Sky last_sky;
	static Now last_now = {NOMJD};
	static double last_ra, last_dec;	/* unprecessed ra/dec */
	double lst, alt, az;
	double ehp, ha, dec;	/* ehp: angular dia of earth from body */


	    double lsn, rsn;
	    double deps, dpsi;

	    last_now = *np;
	    sunpos (np->n_mjd, &lsn, &rsn);		/* sun's true ecliptic long
						 * and dist
						 */
	    nutation (np->n_mjd, &deps, &dpsi);	/* correct for nutation */
	    lsn += dpsi;
	    lsn -= degrad(20.4/3600);		/* and light travel time */

	    sp->s_edist = rsn;
	    sp->s_sdist = 0.0;
	    sp->s_elong = 0.0;
	    sp->s_size = raddeg(4.65242e-3/rsn)*3600*2;
	    sp->s_mag = -26.8;
	    sp->s_hlong = lsn-PI;	/* geo- to helio- centric */
	    range (&sp->s_hlong, 2*PI);
	    sp->s_hlat = 0.0;

	    ecl_eq (np->n_mjd, 0.0, lsn, &last_ra, &last_dec);
	    sp->s_ra = last_ra;
	    sp->s_dec = last_dec;
	    if (np->n_epoch != EOD)
		precess (np->n_mjd, np->n_epoch, &sp->s_ra, &sp->s_dec);
		
	now_lst (np, &lst);
	ha = hrrad(lst) - last_ra;
	ehp = (2.0 * 6378.0 / 146.0e6) / sp->s_edist;
	ta_par (ha, last_dec, np->n_lat, np->n_height, ehp, &ha, &dec);
	hadec_aa (np->n_lat, ha, dec, &alt, &az);
	refract (np->n_pressure, np->n_temp, alt, &alt);
	sp->s_alt = alt;
	sp->s_az = az;
	last_sky = *sp;
}

/* find moon's circumstances now.
 * as is the desired accuracy, in arc seconds; use 0.0 for best possible.
 * return 0 if only alt/az changes, else 1 if all other stuff updated too.
 */
void CSunAndMoon::moon_cir (double as, Now *np, Sky *sp)
{
	static Sky last_sky;
	static Now last_now = {NOMJD};
	static double ehp;
	static double last_ra, last_dec;	/* unprecessed */
	double lst, alt, az;
	double ha, dec;

	    double lam, bet;
	    double deps, dpsi;
	    double lsn, rsn;	/* sun long in rads, earth-sun dist in au */
	    double edistau;	/* earth-moon dist, in au */
	    double el;		/* elongation, rads east */

	    last_now = *np;
	    moon (np->n_mjd, &lam, &bet, &ehp);	/* moon's true ecliptic loc */
	    nutation (np->n_mjd, &deps, &dpsi);	/* correct for nutation */
	    lam += dpsi;
	    range (&lam, 2*PI);

	    sp->s_edist = 6378.14/sin(ehp);	/* earth-moon dist, want km */
	    sp->s_size = 3600*31.22512*sin(ehp);/* moon angular dia, seconds */
	    sp->s_hlong = lam;			/* save geo in helio fields */
	    sp->s_hlat = bet;

	    ecl_eq (np->n_mjd, bet, lam, &last_ra, &last_dec);
	    sp->s_ra = last_ra;
	    sp->s_dec = last_dec;
	    if (np->n_epoch != EOD)
		precess (np->n_mjd, np->n_epoch, &sp->s_ra, &sp->s_dec);

	    sunpos (np->n_mjd, &lsn, &rsn);
	    range (&lsn, 2*PI);
	    elongation (lam, bet, lsn, &el);

	    /* solve triangle of earth, sun, and elongation for moon-sun dist */
	    edistau = sp->s_edist/1.495979e8; /* km -> au */
	    sp->s_sdist =
		sqrt (edistau*edistau + rsn*rsn - 2.0*edistau*rsn*cos(el));

	    /* TODO: improve mag; this is based on a flat moon model. */
	    sp->s_mag = -12.7 + 2.5*(log10(PI) - log10(PI/2*(1+1.e-6-cos(el))));

	    sp->s_elong = raddeg(el);	/* want degrees */
	    sp->s_phase = fabs(el)/PI*100.0;	/* want non-negative % */
	    
	/* show topocentric alt/az by correcting ra/dec for parallax 
	 * as well as refraction.
	 */
	now_lst (np, &lst);
	ha = hrrad(lst) - last_ra;
	ta_par (ha, last_dec, np->n_lat, np->n_height, ehp, &ha, &dec);
	hadec_aa (np->n_lat, ha, dec, &alt, &az);
	refract (np->n_pressure, np->n_temp, alt, &alt);
	sp->s_alt = alt;
	sp->s_az = az;
	last_sky = *sp;
}

/* given the modified JD, mjd, return the true geocentric ecliptic longitude
 *   of the sun for the mean equinox of the date, *lsn, in radians, and the
 *   sun-earth distance, *rsn, in AU. (the true ecliptic latitude is never more
 *   than 1.2 arc seconds and so may be taken to be a constant 0.)
 * if the APPARENT ecliptic longitude is required, correct the longitude for
 *   nutation to the true equinox of date and for aberration (light travel time,
 *   approximately  -9.27e7/186000/(3600*24*365)*2*pi = -9.93e-5 radians).
 */
void CSunAndMoon::sunpos (double mjd, double *lsn, double *rsn)
{
	double t, t2;
	double ls, ms;    /* mean longitude and mean anomoay */
	double s, nu, ea; /* eccentricity, true anomaly, eccentric anomaly */
	double a, b, a1, b1, c1, d1, e1, h1, dl, dr;

	t = mjd/36525.;
	t2 = t*t;
	a = 100.0021359*t;
	b = 360.*(a-(long)a);
	ls = 279.69668+.0003025*t2+b;
	a = 99.99736042000039*t;
	b = 360*(a-(long)a);
	ms = 358.47583-(.00015+.0000033*t)*t2+b;
	s = .016751-.0000418*t-1.26e-07*t2;
	anomaly (degrad(ms), s, &nu, &ea);
	a = 62.55209472000015*t;
	b = 360*(a-(long)a);
	a1 = degrad(153.23+b);
	a = 125.1041894*t;
	b = 360*(a-(long)a);
	b1 = degrad(216.57+b);
	a = 91.56766028*t;
	b = 360*(a-(long)a);
	c1 = degrad(312.69+b);
	a = 1236.853095*t;
	b = 360*(a-(long)a);
	d1 = degrad(350.74-.00144*t2+b);
	e1 = degrad(231.19+20.2*t);
	a = 183.1353208*t;
	b = 360*(a-(long)a);
	h1 = degrad(353.4+b);
	dl = .00134*cos(a1)+.00154*cos(b1)+.002*cos(c1)+.00179*sin(d1)+
								.00178*sin(e1);
	dr = 5.43e-06*sin(a1)+1.575e-05*sin(b1)+1.627e-05*sin(c1)+
					    3.076e-05*cos(d1)+9.27e-06*sin(h1);
	*lsn = nu+degrad(ls-ms+dl);
	*rsn = 1.0000002*(1-s*cos(ea))+dr;
	range (lsn, 2*PI);
}
static char *dbfile;			/* !0 if set by -d option */
static char dbfdef[] = "ephem.db"; 	/* default database file name */

/* structures to describe objects of various types.
 */
#define	MAXNM		16	/* longest allowed object name, inc \0 */
typedef struct {
    double f_ra;	/* ra, rads, at given epoch */
    double f_dec;	/* dec, rads, at given epoch */
    double f_mag;	/* visual magnitude */
    double f_epoch;	/* the given epoch, as an mjd */
    char f_name[MAXNM];	/* name */
} ObjF;			/* fixed object */
typedef struct {
    double e_inc;	/* inclination, degrees */
    double e_Om;	/* longitude of ascending node, degrees */
    double e_om;	/* argument of perihelion, degress */
    double e_a;		/* mean distance, aka, semi-maj axis, in AU */
    double e_n;		/* daily motion, degrees/day */
    double e_e;		/* eccentricity */
    double e_M;		/* mean anomaly, ie, degrees from perihelion at... */
    double e_cepoch;	/* epoch date (M reference), as an mjd */
    double e_epoch;	/* equinox year (inc/Om/om reference), as an mjd */
    double e_m1, e_m2;	/* magnitude model coefficients: H/G or g/k */
    int e_whichm;	/* MAG_HG (default) or MAG_gk */
    char e_name[MAXNM];	/* name */
} ObjE;			/* object in heliocentric elliptical orbit */
typedef struct {
    double p_ep;	/* epoch of perihelion, as an mjd */
    double p_inc;	/* inclination, degs */
    double p_qp;	/* perihelion distance, AU */
    double p_ap;	/* argument of perihelion, degs. */
    double p_om;	/* longitude of ascending node, degs */
    double p_epoch;	/* reference epoch, as an mjd */
    double p_g, p_k;	/* magnitude model coefficients */
    char p_name[MAXNM];	/* name */
} ObjP;			/* object in heliocentric parabolic trajectory */

typedef struct {
    int  o_type;	/* current object type; see flags, below */
    int  o_on;		/* !=0 while current object is active */
    ObjF o_f;		/* the fixed object */
    ObjE o_e;		/* the elliptical orbit object */
    ObjP o_p;		/* the parabolic orbit object */
} Obj;
#define	FIXED		1
#define	ELLIPTICAL	2
#define	PARABOLIC	3
#define	MAG_HG		0	/* using 0 makes HG the initial default */
#define	MAG_gk		1

static Obj objx;
static Obj objy;

#define	DY	0		/* decimal year flag for set_year() */
#define	YMD	1		/* year/mon/day flag for set_year() */

/* fill in info about object x or y.
 * most arguments and conditions are the same as for plans().
 * only difference is that mag is already apparent, not absolute magnitude.
 * this is called by body_cir() for object x and y just like plans() is called
 * for the planets.
 */
void CSunAndMoon::obj_cir (double jd, int p, double *lpd0, double *psi0, double *rp0, 
		double *rho0, double *lam, double *bet, double *mag)
/*double jd;	/* mjd now *
int p;		/* OBJX or OBJY *
double *lpd0;	/* heliocentric longitude, or NOHELIO  *
double *psi0;	/* heliocentric latitude, or 0 if *lpd0 set to NOHELIO *
double *rp0;	/* distance from the sun, or 0 *
double *rho0;	/* true distance from the Earth, or 0 *
double *lam;	/* apparent geocentric ecliptic longitude *
double *bet;	/* apparent geocentric ecliptic latitude *
double *mag;	/* APPARENT magnitude */
{
	Obj *op =  &objy;

	switch (op->o_type) {
	case FIXED: {
	    double xr, xd;
	    xr = op->o_f.f_ra;
	    xd = op->o_f.f_dec;
	    if (op->o_f.f_epoch != jd)
		precess (op->o_f.f_epoch, jd, &xr, &xd);
	    eq_ecl (jd, xr, xd, bet, lam);

	    *lpd0 = NOHELIO;
	    *psi0 = *rp0 = *rho0 = 0.0;
	    *mag = op->o_f.f_mag;
	    }
	    break;
	case PARABOLIC: {
	    double inc, ap, om;
	    double lpd, psi, rp, rho;
	    double dt;
	    int pass;
	    /* two passes to correct lam and bet for light travel time. */
	    dt = 0.0;
	    for (pass = 0; pass < 2; pass++) {
		reduce_elements (op->o_p.p_epoch, jd-dt, degrad(op->o_p.p_inc),
		    degrad(op->o_p.p_ap), degrad(op->o_p.p_om), &inc, &ap, &om);
		comet (jd-dt, op->o_p.p_ep, inc, ap, op->o_p.p_qp, om,
					&lpd, &psi, &rp, &rho, lam, bet);
		if (pass == 0) {
		    *lpd0 = lpd;
		    *psi0 = psi;
		    *rp0 = rp;
		    *rho0 = rho;
		}
		dt = rho*5.775518e-3;	/* au to light-days */
	    }
	    *mag = op->o_p.p_g + 5*log10(*rho0) + 2.5*op->o_p.p_k*log10(*rp0);
	    }
	    break;
	case ELLIPTICAL: {
	    /* this is basically the same code as pelement() and plans()
	     * combined and simplified for the special case of osculating
	     * (unperturbed) elements.
	     * inputs have been changed to match the Astronomical Almanac.
	     * we have added reduction of elements using reduce_elements().
	     */
	    double dt, lg, lsn, rsn;
	    double nu, ea;
	    double ma, rp, lo, slo, clo;
	    double inc, psi, spsi, cpsi;
	    double y, lpd, rpd, ll, rho, sll, cll;
	    double om;		/* arg of perihelion */
	    double Om;		/* long of ascending node. */
	    double e;
	    int pass;

	    dt = 0;
	    sunpos (jd, &lsn, &rsn);
	    lg = lsn + PI;
	    e = op->o_e.e_e;

	    for (pass = 0; pass < 2; pass++) {

		reduce_elements (op->o_e.e_epoch, jd-dt, degrad(op->o_e.e_inc),
				degrad (op->o_e.e_om), degrad (op->o_e.e_Om),
				&inc, &om, &Om);

		ma = degrad (op->o_e.e_M
				+ (jd - op->o_e.e_cepoch - dt) * op->o_e.e_n);
		anomaly (ma, e, &nu, &ea);
		rp= op->o_e.e_a * (1-e*e) / (1+e*cos(nu));
		lo = nu + om;
		slo = sin(lo);
		clo = cos(lo);
		spsi = slo*sin(inc);
		y = slo*cos(inc);
		psi = asin(spsi);
		lpd = atan(y/clo)+Om;
		if (clo<0) lpd += PI;
		range (&lpd, 2*PI);
		cpsi = cos(psi);
		rpd = rp*cpsi;
		ll = lpd-lg;
		rho = sqrt(rsn*rsn+rp*rp-2*rsn*rp*cpsi*cos(ll));
		dt = rho*5.775518e-3;	/* light travel time, in days */
		if (pass == 0) {
		    *lpd0 = lpd;
		    *psi0 = psi;
		    *rp0 = rp;
		    *rho0 = rho;
		}
	    }

	    sll = sin(ll);
	    cll = cos(ll);
	    if (rpd < rsn)
		*lam = atan(-1*rpd*sll/(rsn-rpd*cll))+lg+PI;
	    else
		*lam = atan(rsn*sll/(rpd-rsn*cll))+lpd;
	    range (lam, 2*PI);
	    *bet = atan(rpd*spsi*sin(*lam-lpd)/(cpsi*rsn*sll));

	    if (op->o_e.e_whichm == MAG_HG) {
		/* this is for the H and G parameters from the Astro. Almanac.
		 */
		double psi_t, Psi_1, Psi_2, beta;
		beta = acos((rp*rp + rho*rho - rsn*rsn)/ (2*rp*rho));
		psi_t = exp(log(tan(beta/2.0))*0.63);
		Psi_1 = exp(-3.33*psi_t);
		psi_t = exp(log(tan(beta/2.0))*1.22);
		Psi_2 = exp(-1.87*psi_t);
		*mag = op->o_e.e_m1 + 5.0*log10(rp*rho)
		    - 2.5*log10((1-op->o_e.e_m2)*Psi_1 + op->o_e.e_m2*Psi_2);

	    } else {
		/* this uses the g/k model of comets */
		*mag =
		  op->o_e.e_m1 + 5*log10(*rho0) + 2.5*op->o_e.e_m2*log10(*rp0);
	    }
	    }
	    break;
	default:
	    break;
	}
}

#define	EQtoECL	1
#define	ECLtoEQ	(-1)

/* given the modified Julian date, mjd, and an equitorial ra and dec, each in
 * radians, find the corresponding geocentric ecliptic latitude, *lat, and
 * longititude, *lng, also each in radians.
 * correction for the effect on the angle of the obliquity due to nutation is
 * included.
 */
void CSunAndMoon::eq_ecl (double mjd, double ra, double dec, double *lat, double *lng)
{
	ecleq_aux (EQtoECL, mjd, ra, dec, lng, lat);
}

/* given the modified Julian date, mjd, and a geocentric ecliptic latitude,
 * *lat, and longititude, *lng, each in radians, find the corresponding
 * equitorial ra and dec, also each in radians.
 * correction for the effect on the angle of the obliquity due to nutation is
 * included.
 */
void CSunAndMoon::ecl_eq (double mjd, double lat, double lng, double *ra, double *dec)
{
	ecleq_aux (ECLtoEQ, mjd, lng, lat, ra, dec);
}

void CSunAndMoon::ecleq_aux (int sw, double mjd, double x, double y, double *p, double *q)
/*int sw;			/* +1 for eq to ecliptic, -1 for vv. 
double mjd, x, y;	/* sw==1: x==ra, y==dec.  sw==-1: x==lng, y==lat. 
double *p, *q;		/* sw==1: p==lng, q==lat. sw==-1: p==ra, q==dec. */
{
	double lastmjd = -10000;	/* last mjd calculated */
	double seps, ceps;	/* sin and cos of mean obliquity */
	double sx, cx, sy, cy, ty;

	if (mjd != lastmjd) {
	    double eps;
	    double deps, dpsi;
	    obliquity (mjd, &eps);		/* mean obliquity for date */
	    nutation (mjd, &deps, &dpsi);
	    eps += deps;
    	    seps = sin(eps);
	    ceps = cos(eps);
	    lastmjd = mjd;
	}

	sy = sin(y);
	cy = cos(y);				/* always non-negative */
        if (fabs(cy)<1e-20) cy = 1e-20;		/* insure > 0 */
        ty = sy/cy;
	cx = cos(x);
	sx = sin(x);
        *q = asin((sy*ceps)-(cy*seps*sx*sw));
        *p = atan(((sx*ceps)+(ty*seps*sw))/cx);
        if (cx<0) *p += PI;		/* account for atan quad ambiguity */
	range (p, 2*PI);
}

/* given the modified Julian date, mjd, find the obliquity of the
 * ecliptic, *eps, in radians.
 */
void CSunAndMoon::obliquity (double mjd, double *eps)
{
	static double lastmjd = -10000, lasteps;

	if (mjd != lastmjd) {
	    double t;
	    t = mjd/36525.;
	    lasteps = degrad(2.345229444E1
			- ((((-1.81E-3*t)+5.9E-3)*t+4.6845E1)*t)/3600.0);
	    lastmjd = mjd;
	}
	*eps = lasteps;
}


/* convert those orbital elements that change from epoch mjd0 to epoch mjd.
 */
void CSunAndMoon::reduce_elements (double mjd0, double mjd, double inc0, double ap0, double om0, 
				double *inc, double *ap, double *om)
/*double mjd0;	/* initial epoch 
double mjd;	/* desired epoch 
double inc0;	/* initial inclination, rads 
double ap0;	/* initial argument of perihelion, as an mjd 
double om0;	/* initial long of ascending node, rads 
double *inc;	/* desired inclination, rads 
double *ap;	/* desired epoch of perihelion, as an mjd 
double *om;	/* desired long of ascending node, rads */
{
	double t0, t1;
	double tt, tt2, t02, tt3;
	double eta, th, th0;
	double a, b;
	double dap;
	double cinc, sinc;
	double ot, sot, cot, ot1;
	double seta, ceta;

	t0 = mjd0/365250.0;
	t1 = mjd/365250.0;

	tt = t1-t0;
	tt2 = tt*tt;
        t02 = t0*t0;
	tt3 = tt*tt2;
        eta = (471.07-6.75*t0+.57*t02)*tt+(.57*t0-3.37)*tt2+.05*tt3;
        th0 = 32869.0*t0+56*t02-(8694+55*t0)*tt+3*tt2;
        eta = degrad(eta/3600.0);
        th0 = degrad((th0/3600.0)+173.950833);
        th = (50256.41+222.29*t0+.26*t02)*tt+(111.15+.26*t0)*tt2+.1*tt3;
        th = th0+degrad(th/3600.0);
	cinc = cos(inc0);
        sinc = sin(inc0);
	ot = om0-th0;
	sot = sin(ot);
        cot = cos(ot);
	seta = sin(eta);
        ceta = cos(eta);
	a = sinc*sot;
        b = ceta*sinc*cot-seta*cinc;
	ot1 = atan(a/b);
        if (b<0) ot1 += PI;
        b = sinc*ceta-cinc*seta*cot;
        a = -1*seta*sot;
	dap = atan(a/b);
        if (b<0) dap += PI;

        *ap = ap0+dap;
	range (ap, 2*PI);
        *om = ot1+th;
	range (om, 2*PI);

        if (inc0<.175)
	    *inc = asin(a/sin(dap));
	else
	    *inc = 1.570796327-asin((cinc*ceta)+(sinc*seta*cot));
}


/* given a modified Julian date, mjd, and a set of heliocentric parabolic
 * orbital elements referred to the epoch of date (mjd):
 *   ep:   epoch of perihelion,
 *   inc:  inclination,
 *   ap:   argument of perihelion (equals the longitude of perihelion minus the
 *	   longitude of ascending node)
 *   qp:   perihelion distance,
 *   om:   longitude of ascending node;
 * find:
 *   lpd:  heliocentric longitude, 
 *   psi:  heliocentric latitude,
 *   rp:   distance from the sun to the planet, 
 *   rho:  distance from the Earth to the planet,
 *   lam:  geocentric ecliptic longitude, 
 *   bet:  geocentric ecliptic latitude,
 *         none are corrected for light time, ie, they are the true values for
 *	   the given instant.
 *
 * all angles are in radians, all distances in AU.
 * mutual perturbation corrections with other solar system objects are not
 * applied. corrections for nutation and abberation must be made by the caller.
 * The RA and DEC calculated from the fully-corrected ecliptic coordinates are
 * then the apparent geocentric coordinates. Further corrections can be made,
 * if required, for atmospheric refraction and geocentric parallax.
 */
void CSunAndMoon::comet (double mjd, double ep, double inc, double ap, double qp, double om, 
		double *lpd, double *psi, double *rp, double *rho, double *lam, double *bet)
{
	double w, s, s2;
	double l, sl, cl, y;
	double spsi, cpsi;
	double rd, lsn, rsn;
	double lg, re, ll;
	double cll, sll;
	double nu;

#define	ERRLMT	0.0001
        w = ((mjd-ep)*3.649116e-02)/(qp*sqrt(qp));
        s = w/3;
	while (1) {
	    double d;
	    s2 = s*s;
	    d = (s2+3)*s-w;
	    if (fabs(d) <= ERRLMT)
		break;
	    s = ((2*s*s2)+w)/(3*(s2+1));
	}

        nu = 2*atan(s);
	*rp = qp*(1+s2);
	l = nu+ap;
        sl = sin(l);
	cl = cos(l);
	spsi = sl*sin(inc);
        *psi = asin(spsi);
	y = sl*cos(inc);
        *lpd = atan(y/cl)+om;
	cpsi = cos(*psi);
        if (cl<0) *lpd += PI;
	range (lpd, 2*PI);
        rd = *rp * cpsi;
	sunpos (mjd, &lsn, &rsn);
	lg = lsn+PI;
        re = rsn;
	ll = *lpd - lg;
        cll = cos(ll);
	sll = sin(ll);
        *rho = sqrt((re * re)+(*rp * *rp)-(2*re*rd*cll));
        if (rd<re) 
            *lam = atan((-1*rd*sll)/(re-(rd*cll)))+lg+PI;
	else
	    *lam = atan((re*sll)/(rd-(re*cll)))+*lpd;
	range (lam, 2*PI);
        *bet = atan((rd*spsi*sin(*lam-*lpd))/(cpsi*re*sll));
}


#define	TWOPI	(2*PI)

/* given the mean anomaly, ma, and the eccentricity, s, of elliptical motion,
 * find the true anomaly, *nu, and the eccentric anomaly, *ea.
 * all angles in radians.
 */
void CSunAndMoon::anomaly (double ma, double s, double *nu, double *ea)
{
	double m, dla, fea;

	m = ma-TWOPI*(long)(ma/TWOPI);
	fea = m;
	while (1) {
	    dla = fea-(s*sin(fea))-m;
	    if (fabs(dla)<1e-6)
		break;
	    dla /= 1-(s*cos(fea));
	    fea -= dla;
	}

	*nu = 2*atan(sqrt((1+s)/(1-s))*tan(fea/2));
	*ea = fea;
}

#define	mod2PI(x)	((x) - (long)((x)/TWOPI)*TWOPI)

/* given a modified Julian date, mjd, and a planet, p, find:
 *   lpd0: heliocentric longitude, 
 *   psi0: heliocentric latitude,
 *   rp0:  distance from the sun to the planet, 
 *   rho0: distance from the Earth to the planet,
 *         none corrected for light time, ie, they are the true values for the
 *         given instant.
 *   lam:  geocentric ecliptic longitude, 
 *   bet:  geocentric ecliptic latitude,
 *         each corrected for light time, ie, they are the apparent values as
 *	   seen from the center of the Earth for the given instant.
 *   dia:  angular diameter in arcsec at 1 AU, 
 *   mag:  visual magnitude when 1 AU from sun and earth at 0 phase angle.
 *
 * all angles are in radians, all distances in AU.
 * the mean orbital elements are found by calling pelement(), then mutual
 *   perturbation corrections are applied as necessary.
 *
 * corrections for nutation and abberation must be made by the caller. The RA 
 *   and DEC calculated from the fully-corrected ecliptic coordinates are then
 *   the apparent geocentric coordinates. Further corrections can be made, if
 *   required, for atmospheric refraction and geocentric parallax although the
 *   intrinsic error herein of about 10 arcseconds is usually the dominant
 *   error at this stage.
 * TODO: combine the several intermediate expressions when get a good compiler.
 */
void CSunAndMoon::plans (double mjd, int p, double *lpd0, double *psi0, double *rp0, 
		double *rho0, double *lam, double *bet, double *dia, double *mag)
{
	static double plan[8][9];
	static double lastmjd = -10000;
	double dl;	/* perturbation correction for longitude */
	double dr;	/*  "   orbital radius */
	double dml;	/*  "   mean longitude */
	double ds;	/*  "   eccentricity */
	double dm;	/*  "   mean anomaly */
	double da;	/*  "   semi-major axis */
	double dhl;	/*  "   heliocentric longitude */
	double lsn, rsn;/* true geocentric longitude of sun and sun-earth rad */
	double mas;	/* mean anomaly of the sun */
	double re;	/* radius of earth's orbit */
	double lg;	/* longitude of earth */
	double map[8];	/* array of mean anomalies for each planet */
	double lpd, psi, rp, rho;
	double ll, sll, cll;
	double t;
	double dt;
	int pass;
	int j;
	double s, ma;
	double nu, ea;
	double lp, om;
	double lo, slo, clo;
	double inc, y;
	double spsi, cpsi;
	double rpd;

	/* only need to fill in plan[] once for a given mjd */
	if (mjd != lastmjd) {
	    pelement (mjd, plan);
	    lastmjd = mjd;
	}

	dt = 0;
	t = mjd/36525.;
	sunpos (mjd, &lsn, &rsn);
	masun (mjd, &mas);
        re = rsn;
	lg = lsn+PI;

	/* first find the true position of the planet at mjd.
	 * then repeat a second time for a slightly different time based
	 * on the position found in the first pass to account for light-travel
	 * time.
	 */
	for (pass = 0; pass < 2; pass++) {

	    for (j = 0; j < 8; j++)
		map[j] = degrad(plan[j][0]-plan[j][2]-dt*plan[j][1]);

	    /* set initial corrections to 0.
	     * then modify as necessary for the planet of interest.
	     */
	    dl = 0;
	    dr = 0;
	    dml = 0;
	    ds = 0;
	    dm = 0;
	    da = 0;
	    dhl = 0;

	    s = plan[p][3]+ds;
	    ma = map[p]+dm;
	    anomaly (ma, s, &nu, &ea);
	    rp = (plan[p][6]+da)*(1-s*s)/(1+s*cos(nu));
	    lp = raddeg(nu)+plan[p][2]+raddeg(dml-dm);
	    lp = degrad(lp);
	    om = degrad(plan[p][5]);
	    lo = lp-om;
	    slo = sin(lo);
	    clo = cos(lo);
	    inc = degrad(plan[p][4]);
	    rp = rp+dr;
	    spsi = slo*sin(inc);
	    y = slo*cos(inc);
	    psi = asin(spsi)+dhl;
	    spsi = sin(psi);
	    lpd = atan(y/clo)+om+degrad(dl);
	    if (clo<0) lpd += PI;
	    range (&lpd, TWOPI);
	    cpsi = cos(psi);
	    rpd = rp*cpsi;
	    ll = lpd-lg;
	    rho = sqrt(re*re+rp*rp-2*re*rp*cpsi*cos(ll));

	    /* when we view a planet we see it in the position it occupied
	     * dt days ago, where rho is the distance between it and earth,
	     * in AU. use this as the new time for the next pass.
	     */
	    dt = rho*5.775518e-3;

	    if (pass == 0) {
		/* save heliocentric coordinates after first pass since, being
		 * true, they are NOT to be corrected for light-travel time.
		 */
		*lpd0 = lpd;
		range (lpd0, TWOPI);
		*psi0 = psi;
		*rp0 = rp;
		*rho0 = rho;
	    }
	}

        sll = sin(ll);
	cll = cos(ll);
        if (p < MARS) 
	    *lam = atan(-1*rpd*sll/(re-rpd*cll))+lg+PI;
	else
	    *lam = atan(re*sll/(rpd-re*cll))+lpd;
	range (lam, TWOPI);
        *bet = atan(rpd*spsi*sin(*lam-lpd)/(cpsi*re*sll));
	*dia = plan[p][7];
	*mag = plan[p][8];
}

/* this array contains polynomial coefficients to find the various orbital
 *   elements for the mean orbit at any instant in time for each major planet.
 *   the first five elements are in the form a0 + a1*t + a2*t**2 + a3*t**3,
 *   where t is the number of Julian centuries of 36525 Julian days since 1900
 *   Jan 0.5. the last three elements are constants.
 *
 * the orbital element (column) indeces are:
 *   [ 0- 3]: coefficients for mean longitude, in degrees;
 *   [ 4- 7]: coefficients for longitude of the perihelion, in degrees;
 *   [ 8-11]: coefficients for eccentricity;
 *   [12-15]: coefficients for inclination, in degrees;
 *   [16-19]: coefficients for longitude of the ascending node, in degrees;
 *      [20]: semi-major axis, in AU;
 *      [21]: angular diameter at 1 AU, in arcsec;
 *      [22]: standard visual magnitude, ie, the visual magnitude of the planet
 *	      when at a distance of 1 AU from both the Sun and the Earth and
 *	      with zero phase angle.
 *
 * the planent (row) indeces are:
 *   [0]: Mercury; [1]: Venus;   [2]: Mars;  [3]: Jupiter; [4]: Saturn;
 *   [5]: Uranus;  [6]: Neptune; [7]: Pluto.
 */
#define	NPELE	(5*4 + 3)	/* 4 coeffs for ea of 5 elems, + 3 constants */
static double elements[8][NPELE] = {

	{   /*     mercury... */

	    178.179078,	415.2057519,	3.011e-4,	0.0,
	    75.899697,	1.5554889,	2.947e-4,	0.0,
	    .20561421,	2.046e-5,	3e-8,		0.0,
	    7.002881,	1.8608e-3,	-1.83e-5,	0.0,
	    47.145944,	1.1852083,	1.739e-4,	0.0,
	    .3870986,	6.74, 		-0.42
	},

	{   /*     venus... */

	    342.767053,	162.5533664,	3.097e-4,	0.0,
	    130.163833,	1.4080361,	-9.764e-4,	0.0,
	    6.82069e-3,	-4.774e-5,	9.1e-8,		0.0,
	    3.393631,	1.0058e-3,	-1e-6,		0.0,
	    75.779647,	.89985,		4.1e-4,		0.0,
	    .7233316,	16.92,		-4.4
	},

	{   /*     mars... */

	    293.737334,	53.17137642,	3.107e-4,	0.0,
	    3.34218203e2, 1.8407584,	1.299e-4,	-1.19e-6,
	    9.33129e-2,	9.2064e-5,	7.7e-8,		0.0,
	    1.850333,	-6.75e-4,	1.26e-5,	0.0,
	    48.786442,	.7709917,	-1.4e-6,	-5.33e-6,
	    1.5236883,	9.36,		-1.52
	},

	{   /*     jupiter... */

	    238.049257,	8.434172183,	3.347e-4,	-1.65e-6,
	    1.2720972e1, 1.6099617,	1.05627e-3,	-3.43e-6,
	    4.833475e-2, 1.6418e-4,	-4.676e-7,	-1.7e-9,
	    1.308736,	-5.6961e-3,	3.9e-6,		0.0,
	    99.443414,	1.01053,	3.5222e-4,	-8.51e-6,
	    5.202561,	196.74,		-9.4
	},

	{   /*     saturn... */

	    266.564377,	3.398638567,	3.245e-4,	-5.8e-6,
	    9.1098214e1, 1.9584158,	8.2636e-4,	4.61e-6,
	    5.589232e-2, -3.455e-4,	-7.28e-7,	7.4e-10,
	    2.492519,	-3.9189e-3,	-1.549e-5,	4e-8,
	    112.790414,	.8731951,	-1.5218e-4,	-5.31e-6,
	    9.554747,	165.6,		-8.88
	},

	{   /*     uranus... */

	    244.19747,	1.194065406,	3.16e-4,	-6e-7,
	    1.71548692e2, 1.4844328,	2.372e-4,	-6.1e-7,
	    4.63444e-2,	-2.658e-5,	7.7e-8,		0.0,
	    .772464,	6.253e-4,	3.95e-5,	0.0,
	    73.477111,	.4986678,	1.3117e-3,	0.0,
	    19.21814,	65.8,		-7.19
	},

	{   /*     neptune... */

	    84.457994,	.6107942056,	3.205e-4,	-6e-7,
	    4.6727364e1, 1.4245744,	3.9082e-4,	-6.05e-7,
	    8.99704e-3,	6.33e-6,	-2e-9,		0.0,
	    1.779242,	-9.5436e-3,	-9.1e-6,	0.0,
	    130.681389,	1.098935,	2.4987e-4,	-4.718e-6,
	    30.10957,	62.2,		-6.87
	},

	{   /*     pluto...(osculating 1984 jan 21) */

	    95.3113544,	.3980332167,	0.0,		0.0,
	    224.017,	0.0,		0.0,		0.0,
	    .25515,	0.0,		0.0,		0.0,
	    17.1329,	0.0,		0.0,		0.0,
	    110.191,	0.0,		0.0,		0.0,
	    39.8151,	8.2,		-1.0
	}
};

/* given a modified Julian date, mjd, return the elements for the mean orbit
 *   at that instant of all the major planets, together with their
 *   mean daily motions in longitude, angular diameter and standard visual
 *   magnitude.
 * plan[i][j] contains all the values for all the planets at mjd, such that
 *   i = 0..7: mercury, venus, mars, jupiter, saturn, unranus, neptune, pluto;
 *   j = 0..8: mean longitude, mean daily motion in longitude, longitude of 
 *     the perihelion, eccentricity, inclination, longitude of the ascending
 *     node, length of the semi-major axis, angular diameter from 1 AU, and
 *     the standard visual magnitude (see elements[][] comment, above).
 */
void CSunAndMoon::pelement (double mjd, double plan[8][9])
{
	register double *ep, *pp;
	register double t = mjd/36525.;
	double aa;
	int planet, i;


	for (planet = 0; planet < 8; planet++) {
	    ep = elements[planet];
	    pp = plan[planet];
	    aa = ep[1]*t;
	    pp[0] = ep[0] + 360.*(aa-(long)aa) + (ep[3]*t + ep[2])*t*t;
	    range (pp, 360.);
	    pp[1] = (ep[1]*9.856263e-3) + (ep[2] + ep[3])/36525;

	    for (i = 4; i < 20; i += 4)
		pp[i/4+1] = ((ep[i+3]*t + ep[i+2])*t + ep[i+1])*t + ep[i+0];

	    pp[6] = ep[20];
	    pp[7] = ep[21];
	    pp[8] = ep[22];
	}
}

/* find the mean anomaly of the sun at mjd.
 * this is the same as that used in sun() but when it was converted to C it
 * was not known it would be required outside that routine.
 * TODO: add an argument to sun() to return mas and eliminate this routine.
 */
void CSunAndMoon::masun (double mjd, double *mas)
{
	double t, t2;
	double a, b;

	t = mjd/36525;
	t2 = t*t;
	a = 9.999736042e1*t;
	b = 360.*(a-(long)a);
	*mas = degrad (3.5847583e2-(1.5e-4+3.3e-6*t)*t2+b);
}

/* given the modified JD, mjd, find the nutation in obliquity, *deps, and
 * the nutation in longitude, *dpsi, each in radians.
 */
void CSunAndMoon::nutation (double mjd, double *deps, double *dpsi)
{
	static double lastmjd = -10000, lastdeps, lastdpsi;
	double ls, ld;	/* sun's mean longitude, moon's mean longitude */
	double ms, md;	/* sun's mean anomaly, moon's mean anomaly */
	double nm;	/* longitude of moon's ascending node */
	double t, t2;	/* number of Julian centuries of 36525 days since
			 * Jan 0.5 1900.
			 */
	double tls, tnm, tld;	/* twice above */
	double a, b;	/* temps */

	if (mjd == lastmjd) {
	    *deps = lastdeps;
	    *dpsi = lastdpsi;
	    return;
	}
	    
	t = mjd/36525.;
	t2 = t*t;

	a = 100.0021358*t;
	b = 360.*(a-(long)a);
	ls = 279.697+.000303*t2+b;

	a = 1336.855231*t;
	b = 360.*(a-(long)a);
	ld = 270.434-.001133*t2+b;

	a = 99.99736056000026*t;
	b = 360.*(a-(long)a);
	ms = 358.476-.00015*t2+b;

	a = 13255523.59*t;
	b = 360.*(a-(long)a);
	md = 296.105+.009192*t2+b;

	a = 5.372616667*t;
	b = 360.*(a-(long)a);
	nm = 259.183+.002078*t2-b;

	/* convert to radian forms for use with trig functions.
	 */
	tls = 2*degrad(ls);
	nm = degrad(nm);
	tnm = 2*degrad(nm);
	ms = degrad(ms);
	tld = 2*degrad(ld);
	md = degrad(md);

	/* find delta psi and eps, in arcseconds.
	 */
	lastdpsi = (-17.2327-.01737*t)*sin(nm)+(-1.2729-.00013*t)*sin(tls)
		   +.2088*sin(tnm)-.2037*sin(tld)+(.1261-.00031*t)*sin(ms)
		   +.0675*sin(md)-(.0497-.00012*t)*sin(tls+ms)
		   -.0342*sin(tld-nm)-.0261*sin(tld+md)+.0214*sin(tls-ms)
		   -.0149*sin(tls-tld+md)+.0124*sin(tls-nm)+.0114*sin(tld-md);
	lastdeps = (9.21+.00091*t)*cos(nm)+(.5522-.00029*t)*cos(tls)
		   -.0904*cos(tnm)+.0884*cos(tld)+.0216*cos(tls+ms)
		   +.0183*cos(tld-nm)+.0113*cos(tld+md)-.0093*cos(tls-ms)
		   -.0066*cos(tls-nm);

	/* convert to radians.
	 */
	lastdpsi = degrad(lastdpsi/3600);
	lastdeps = degrad(lastdeps/3600);

	lastmjd = mjd;
	*deps = lastdeps;
	*dpsi = lastdpsi;
}

/* given geocentric ecliptic longitude and latitude, lam and bet, of some object
 * and the longitude of the sun, lsn, find the elongation, el. this is the
 * actual angular separation of the object from the sun, not just the difference
 * in the longitude. the sign, however, IS set simply as a test on longitude
 * such that el will be >0 for an evening object <0 for a morning object.
 * to understand the test for el sign, draw a graph with lam going from 0-2*PI
 *   down the vertical axis, lsn going from 0-2*PI across the hor axis. then
 *   define the diagonal regions bounded by the lines lam=lsn+PI, lam=lsn and
 *   lam=lsn-PI. the "morning" regions are any values to the lower left of the
 *   first line and bounded within the second pair of lines.
 * all angles in radians.
 */
void CSunAndMoon::elongation (double lam, double bet, double lsn, double *el)
{
	*el = acos(cos(bet)*cos(lam-lsn));
	if (lam>lsn+PI || lam>lsn-PI && lam<lsn) *el = - *el;
}

/* given latitude (n+, radians), lat, altitude (up+, radians), alt, and
 * azimuth (angle round to the east from north+, radians),
 * return hour angle (radians), ha, and declination (radians), dec.
 */
void CSunAndMoon::aa_hadec (double lat, double alt, double az, double *ha, double *dec)
{
	aaha_aux (lat, az, alt, ha, dec);
}

/* given latitude (n+, radians), lat, hour angle (radians), ha, and declination
 * (radians), dec,
 * return altitude (up+, radians), alt, and
 * azimuth (angle round to the east from north+, radians),
 */
void CSunAndMoon::hadec_aa (double lat, double ha, double dec, double *alt, double *az)
{
	aaha_aux (lat, ha, dec, az, alt);
}

/* the actual formula is the same for both transformation directions so
 * do it here once for each way.
 * N.B. all arguments are in radians.
 */
void CSunAndMoon::aaha_aux (double lat, double x, double y, double *p, double *q)
{
	static double lastlat = -1000.;
	static double sinlastlat, coslastlat;
	double sy, cy;
	double sx, cx;
	double sq, cq;
	double a;
	double cp;

	/* latitude doesn't change much, so try to reuse the sin and cos evals.
	 */
	if (lat != lastlat) {
	    sinlastlat = sin (lat);
	    coslastlat = cos (lat);
	    lastlat = lat;
	}

	sy = sin (y);
	cy = cos (y);
	sx = sin (x);
	cx = cos (x);

/* define GOODATAN2 if atan2 returns full range -PI through +PI.
 */
#ifdef GOODATAN2
	*q = asin ((sy*sinlastlat) + (cy*coslastlat*cx));
	*p = atan2 (-cy*sx, -cy*cx*sinlastlat + sy*coslastlat);
#else
#define	EPS2	(1e-20)
	sq = (sy*sinlastlat) + (cy*coslastlat*cx);
	*q = asin (sq);
	cq = cos (*q);
	a = coslastlat*cq;
	if (a > -EPS2 && a < EPS2)
	    a = a < 0 ? -EPS2 : EPS2; /* avoid CSunAndMoon::/ 0 */
	cp = (sy - (sinlastlat*sq))/a;
	if (cp >= 1.0)	/* the /a can be slightly > 1 */
	    *p = 0.0;
	else if (cp <= -1.0)
	    *p = PI;
	else
	    *p = acos ((sy - (sinlastlat*sq))/a);
	if (sx>0) *p = 2.0*PI - *p;
#endif
}

/* find where and when a body, p, will rise and set and
 *   its transit circumstances. all times are local, angles rads e of n.
 * status is set from the RS_* #defines in circum.h.
 * also used to find astro twilight by calling with dis of 18 degrees.
 */
void CSunAndMoon::riset_cir (int p, Now *np, int hzn,double *ltr, 
							 double *lts, double *ltt, double *azr, 
							 double *azs, double *altt, int *status)
{
	typedef struct {
	    Now l_now;
	    double l_ltr, l_lts, l_ltt, l_azr, l_azs, l_altt;
	    int l_hzn;
	    int l_status;
	} Last;
	/* must be in same order as the astro.h/screen.h #define's */
	Last last[12] = {
	    {NOMJD}, {NOMJD}, {NOMJD}, {NOMJD}, {NOMJD}, {NOMJD},
	    {NOMJD}, {NOMJD}, {NOMJD}, {NOMJD}, {NOMJD}, {NOMJD}
	};
	Last *lp;

	lp = last + p;
	
    *status = 0;
    iterative_riset (p, np, hzn, ltr, lts, ltt, azr, azs, altt, status);
    lp->l_ltr = *ltr;
    lp->l_lts = *lts;
    lp->l_ltt = *ltt;
    lp->l_azr = *azr;
    lp->l_azs = *azs;
    lp->l_altt = *altt;
    lp->l_status = *status;
    lp->l_hzn = hzn;
    lp->l_now = *np;
	    
	return;
}


#define	STDREF	degrad(34./60.)	/* nominal horizon refraction amount */
#define	TWIREF	degrad(18.)	/* twilight horizon displacement */
#define	TMACC	(15./3600.)	/* convergence accuracy, hours */

void CSunAndMoon::iterative_riset (int p, Now *np, int hzn, double *ltr, 
								   double *lts, double *ltt, double *azr, 
								   double *azs, double *altt, int *status)
{
#define	MAXPASSES	6
	double lstr, lsts, lstt; /* local sidereal times of rising/setting */
	double mjd0;		/* mjd estimates of rise/set event */
	double lnoon;		/* mjd of local noon */
	double x;		/* discarded tmp value */
	Now n;			/* just used to call now_lst() */
	double lst;		/* lst at local noon */
	double diff, lastdiff;	/* iterative improvement to mjd0 */
	int pass;
	int rss;

	/* first approximation is to find rise/set times of a fixed object
	 * in its position at local noon.
	 */
	lnoon = mjd_day(np->n_mjd - np->n_tz/24.0) + (12.0 + np->n_tz)/24.0; /*mjd of local noon*/
	n.n_mjd = lnoon;
	n.n_lng = np->n_lng;
	now_lst (&n, &lst);	/* lst at local noon */
	mjd0 = lnoon;
	stationary_riset (p,mjd0,np,hzn,&lstr,&lsts,&lstt,&x,&x,&x,&rss);
    chkrss:
	switch (rss) {
	case  0:  break;
	case  1: *status = RS_NEVERUP; return;
	case -1: *status = RS_CIRCUMPOLAR; goto transit;
	default: *status = RS_ERROR; return;
	}

	/* find a better approximation to the rising circumstances based on
	 * more passes, each using a "fixed" object at the location at
	 * previous approximation of the rise time.
	 */
	lastdiff = 1000.0;
	for (pass = 1; pass < MAXPASSES; pass++) {
	    diff = (lstr - lst)*SIDRATE; /* next guess at rise time wrt noon */
	    if (diff > 12.0)
		diff -= 24.0*SIDRATE;	/* not tomorrow, today */
	    else if (diff < -12.0)
		diff += 24.0*SIDRATE;	/* not yesterday, today */
	    mjd0 = lnoon + diff/24.0;	/* next guess at mjd of rise */
	    stationary_riset (p,mjd0,np,hzn,&lstr,&x,&x,azr,&x,&x,&rss);
	    if (rss != 0) goto chkrss;
	    if (fabs (diff - lastdiff) < TMACC)
		break;
	    lastdiff = diff;
	}
	if (pass == MAXPASSES)
	    *status |= RS_NORISE;	/* didn't converge - no rise today */
	else {
	    *ltr = 12.0 + diff;
	    if (p != MOON &&
		    (*ltr <= 24.0*(1.0-SIDRATE) || *ltr >= 24.0*SIDRATE))
		*status |= RS_2RISES;
	}

	/* find a better approximation to the setting circumstances based on
	 * more passes, each using a "fixed" object at the location at
	 * previous approximation of the set time.
	 */
	lastdiff = 1000.0;
	for (pass = 1; pass < MAXPASSES; pass++) {
	    diff = (lsts - lst)*SIDRATE; /* next guess at set time wrt noon */
	    if (diff > 12.0)
		diff -= 24.0*SIDRATE;	/* not tomorrow, today */
	    else if (diff < -12.0)
		diff += 24.0*SIDRATE;	/* not yesterday, today */
	    mjd0 = lnoon + diff/24.0;	/* next guess at mjd of set */
	    stationary_riset (p,mjd0,np,hzn,&x,&lsts,&x,&x,azs,&x,&rss);
	    if (rss != 0) goto chkrss;
	    if (fabs (diff - lastdiff) < TMACC)
		break;
	    lastdiff = diff;
	}
	if (pass == MAXPASSES)
	    *status |= RS_NOSET;	/* didn't converge - no set today */
	else {
	    *lts = 12.0 + diff;
	    if (p != MOON &&
		    (*lts <= 24.0*(1.0-SIDRATE) || *lts >= 24.0*SIDRATE))
		*status |= RS_2SETS;
	}

    transit:
	/* find a better approximation to the transit circumstances based on
	 * more passes, each using a "fixed" object at the location at
	 * previous approximation of the transit time.
	 */
	lastdiff = 1000.0;
	for (pass = 1; pass < MAXPASSES; pass++) {
	    diff = (lstt - lst)*SIDRATE; /*next guess at transit time wrt noon*/
	    if (diff > 12.0)
		diff -= 24.0*SIDRATE;	/* not tomorrow, today */
	    else if (diff < -12.0)
		diff += 24.0*SIDRATE;	/* not yesterday, today */
	    mjd0 = lnoon + diff/24.0;	/* next guess at mjd of transit */
	    stationary_riset (p,mjd0,np,hzn,&x,&x,&lstt,&x,&x,altt,&rss);
	    if (fabs (diff - lastdiff) < TMACC)
		break;
	    lastdiff = diff;
	}
	if (pass == MAXPASSES)
	    *status |= RS_NOTRANS;	/* didn't converge - no transit today */
	else {
	    *ltt = 12.0 + diff;
	    if (p != MOON &&
		    (*ltt <= 24.0*(1.0-SIDRATE) || *ltt >= 24.0*SIDRATE))
		*status |= RS_2TRANS;
	}
}



void CSunAndMoon::stationary_riset (int p, double mjd0, Now *np, int hzn, 
					   double *lstr, double *lsts, double *lstt, 
					   double *azr, double *azs, double *altt, int *status)
{
	/*extern void bye();*/
	double dis;
	Now n;
	Sky s;

	/* find object p's topocentric ra/dec at mjd0
	 * (this must include parallax)
	 */
	n = *np;
	n.n_mjd = mjd0;
	(void) body_cir (p, 0.0, &n, &s);
	if (np->n_epoch != EOD)
	    precess (np->n_epoch, mjd0, &s.s_ra, &s.s_dec);
	if (s.s_edist > 0) {
	    /* parallax, if we can */
	    double ehp, lst, ha;
	    if (p == MOON)
		ehp = asin (6378.14/s.s_edist);
	    else
		ehp = (2.*6378./146e6)/s.s_edist;
	    now_lst (&n, &lst);
	    ha = hrrad(lst) - s.s_ra;
	    ta_par (ha, s.s_dec, np->n_lat, np->n_height, ehp, &ha, &s.s_dec);
	    s.s_ra = hrrad(lst) - ha;
	    range (&s.s_ra, 2*PI);
	}

	switch (hzn) {
	case STDHZN:
	    /* nominal atmospheric refraction.
	     * then add nominal moon or sun semi-diameter, as appropriate.
	     * other objects assumes to be negligibly small.
	     */
	    dis = STDREF;
	    if (p == MOON || p == SUN)
		dis += degrad (32./60./2.);
	    break;
	case TWILIGHT:
	    if (p != SUN) {
		/*bye();*/
	    }
	    dis = TWIREF;
	    break;
	case MANUAL:
		if (p == MOON)
			dis = m_MoonAngle + 2.*degrad(32./60./2.);
		else if (p == SUN)
			dis = m_SunAngle;
		break;
	case ADPHZN:
	    /* adaptive includes actual refraction conditions and also
	     * includes object's semi-diameter.
	     */
	    unrefract (np->n_pressure, np->n_temp, 0.0, &dis);
	    dis = -dis;
	    dis += degrad(s.s_size/3600./2.0);
	    break;
	}

	riset (s.s_ra, s.s_dec, np->n_lat, dis, lstr, lsts, azr, azs, status);
	transit (s.s_ra, s.s_dec, np, lstt, altt);
}


/* find when and how hi object at (r,d) is when it transits. */
void CSunAndMoon::transit (double r, double d, Now *np, double *lstt, double *altt)
{
	*lstt = radhr(r);
	*altt = PI/2 - np->n_lat + d;
	if (*altt > PI/2)
	    *altt = PI - *altt;
	refract (np->n_pressure, np->n_temp, *altt, altt);
}


void CSunAndMoon::SetSunAngle(double val)
{
	m_SunAngle = degrad(val);
}

void CSunAndMoon::SetMoonAngle(double val)
{
	m_MoonAngle = degrad(val);
}
