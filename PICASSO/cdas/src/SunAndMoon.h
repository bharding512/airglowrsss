// SunAndMoon.h: interface for the CSunAndMoon class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SUNANDMOON_H
#define SUNANDMOON_H

/*Moon Stuff */

#ifndef PI
#define	PI		3.141592653589793
#endif

/* conversions among hours (of ra), degrees and radians. */
#define	degrad(x)	((x)*PI/180.)
#define	raddeg(x)	((x)*180./PI)
#define	hrdeg(x)	((x)*15.)
#define	deghr(x)	((x)/15.)
#define	hrrad(x)	degrad(hrdeg(x))
#define	radhr(x)	deghr(raddeg(x))

/* ratio of from synodic (solar) to sidereal (stellar) rate */
#define	SIDRATE		.9972695677

/* manifest names for planets.
 * N.B. must cooincide with usage in pelement.c and plans.c.
 */
#define	MERCURY	0
#define	VENUS	1
#define	MARS	2
#define	JUPITER	3
#define	SATURN	4
#define	URANUS	5
#define	NEPTUNE	6
#define	PLUTO	7

#define	SUN	    8
#define	MOON    9

#define	OBJX	(PLUTO+3)	/* the user-defined object */
#define	OBJY	(PLUTO+4)	/* the user-defined object */
#define	NOBJ	(OBJY+1)	/* total number of objects */

/* More moon stuff */
#define	SPD	(24.0*3600.0)	/* seconds per day */

#define	EOD	(-9786)		/* special epoch flag: use epoch of date */
#define	RTC	(-1324)		/* special tminc flag: use rt clock */
#define	NOMJD	(-58631.)	/* an unlikely mjd for initing static mjd's */
#define	NOHELIO	(-2314)		/* special s_hlong flag: means it and s_hlat are
				 * undefined
				 */

#define	STDHZN		0	/* rise/set times based on nominal conditions */
#define	ADPHZN		1	/* rise/set times based on exact current " */
#define	TWILIGHT	2	/* rise/set times for sun 18 degs below hor */
#define MANUAL		3   /* use values in m_MoonAngle and m_SunAngle */

/* info about our local observing circumstances */
typedef struct {
	double n_mjd;	/* modified Julian date, ie, days since
			 * Jan 0.5 1900 (== 12 noon, Dec 30, 1899), utc.
			 * enough precision to get well better than 1 second.
			 * N.B. if not first member, must move NOMJD inits.
			 */
	double n_lat;	/* latitude, >0 north, rads */
	double n_lng;	/* longitude, >0 east, rads */
	double n_tz;	/* time zone, hrs behind UTC */
	double n_temp;	/* atmospheric temp, degrees C */
	double n_pressure; /* atmospheric pressure, mBar */
	double n_height;	/* height above sea level, earth radii */
	double n_epoch;	/* desired precession display epoch as an mjd, or EOD */
	char n_tznm[4];	/* time zone name; 3 chars or less, always 0 at end */
} Now;
extern double	mjd_day(), mjd_hr();

/* info about where and how we see something in the sky */
typedef struct {
	double s_ra;	/* ra, rads (precessed to n_epoch) */
	double s_dec;	/* dec, rads (precessed to n_epoch) */
	double s_az;	/* azimuth, >0 e of n, rads */
	double s_alt;	/* altitude above topocentric horizon, rads */
	double s_sdist;	/* dist from object to sun, au */
	double s_edist;	/* dist from object to earth, au */
	double s_elong;	/* angular sep between object and sun, >0 if east */
	double s_hlong;	/* heliocentric longitude, rads */
	double s_hlat;	/* heliocentric latitude, rads */
	double s_size;	/* angular size, arc secs */
	double s_phase;	/* phase, % */
	double s_mag;	/* visual magnitude */
} Sky;

/* flags for riset_cir() status */
#define	RS_NORISE	0x001	/* object does not rise as such today */
#define	RS_2RISES	0x002	/* object rises more than once today */
#define	RS_NOSET	0x004	/* object does not set as such today */
#define	RS_2SETS	0x008	/* object sets more than once today */
#define	RS_CIRCUMPOLAR	0x010	/* object stays up all day today */
#define	RS_2TRANS	0x020	/* transits twice in one day */
#define	RS_NEVERUP	0x040	/* object never rises today */
#define	RS_NOTRANS	0x080	/* doesn't transit today */
#define	RS_ERROR	0x100	/* can't figure out times... */

#define MY_ERROR -20
#define MY_UPALL -21
#define MY_DOWNALL -22
#define MY_NOSET -2
#define MY_NORISE -1

class CSunAndMoon  
{
public:	CSunAndMoon();
	void SetMoonAngle(double val);
	void SetSunAngle(double val);
	virtual ~CSunAndMoon();
	double m_SunAngle;
	double m_MoonAngle;	
	void riset_cir (int p, Now *np, int hzn,double *ltr, 
							 double *lts, double *ltt, double *azr, 
							 double *azs, double *altt, int *status);
	void cal_mjd (int mn, double dy, int yr, double *dmjd);


protected:
	void moon (double mjd, double *lam, double *bet, double *hp);
	void riset (double ra, double dec, double lat, double dis, 
			double *lstr, double *lsts, double *azr, double *azs, int *status);
	void range(double *v, double r);
	double mjd_day(double jd);
	double mjd_hr(double jd);
	void now_lst(Now *np, double *lst);
	void utc_gst(double dmjd, double dutc, double *dgst);
	double tnaught(double dmjd);
	void mjd_cal(double dmjd, int *mn, double *dy, int *yr);
	void body_cir (int p, double as, Now *np, Sky *sp);
	void precess (double mjd1, double mjd2, double *ra, double *dec);
	void ta_par (double tha, double tdec, double phi, double ht, double ehp, double *aha, double *adec);
	void refract (double pr, double tr, double ta, double *aa);
	void unrefract (double pr, double tr, double aa, double *ta);
	void sun_cir (double as, Now *np, Sky *sp);
	void moon_cir (double as, Now *np, Sky *sp);
	void sunpos (double mjd, double *lsn, double *rsn);
	void obj_cir (double jd, int p, double *lpd0, double *psi0, double *rp0, 
			double *rho0, double *lam, double *bet, double *mag);
	void eq_ecl (double mjd, double ra, double dec, double *lat, double *lng);
	void ecl_eq (double mjd, double lat, double lng, double *ra, double *dec);
	void ecleq_aux (int sw, double mjd, double x, double y, double *p, double *q);
	void obliquity (double mjd, double *eps);
	void reduce_elements (double mjd0, double mjd, double inc0, double ap0, double om0, 
					double *inc, double *ap, double *om);
	void comet (double mjd, double ep, double inc, double ap, double qp, double om, 
			double *lpd, double *psi, double *rp, double *rho, double *lam, double *bet);
	void anomaly (double ma, double s, double *nu, double *ea);
	void plans (double mjd, int p, double *lpd0, double *psi0, double *rp0, 
			double *rho0, double *lam, double *bet, double *dia, double *mag);
	void pelement (double mjd, double plan[8][9]);
	void masun (double mjd, double *mas);
	void nutation (double mjd, double *deps, double *dpsi);
	void elongation (double lam, double bet, double lsn, double *el);
	void aa_hadec (double lat, double alt, double az, double *ha, double *dec);
	void hadec_aa (double lat, double ha, double dec, double *alt, double *az);
	void aaha_aux (double lat, double x, double y, double *p, double *q);

							 
	void iterative_riset (int p, Now *np, int hzn, double *ltr, 
								   double *lts, double *ltt, double *azr, 
								   double *azs, double *altt, int *status);

	void stationary_riset (int p, double mjd0, Now *np, int hzn, 
					   double *lstr, double *lsts, double *lstt, 
					   double *azr, double *azs, double *altt, int *status);
	void transit (double r, double d, Now *np, double *lstt, double *altt);

};

#endif
