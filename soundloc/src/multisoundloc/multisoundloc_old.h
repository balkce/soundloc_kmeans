/*
	Source coded by Caleb Rascon and Luis Miguel Gato
	IIMAS, UNAM
	México
	
	This library, using JACK, obtains the angle of different sources in the range of [-180 180].
	
	It requires three channels of independent directional microphones on top of the robot, set at 0 dB gain.
	The following is a simple diagram of their assumed setup:
	
                                                    (front mic)
                                                         ^
                                                         |
                                                         |
                                                        / \
                                                       /   \
                                                      ˇ     ˇ
                                             (left mic)     (right mic)


	The microphones should be set up in a triangle, with an equal distance between their ends and
	a 60 degree angle between left-front, front-right, and right-left.
	
	The angle of one source detected from each pair is compared with the other two to decide which
	provides the best resolution (e.g. the one most in front of a pair, e.g. the smallest angle of all).
	Then, a redundancy check is carried out to see if all the pairs are detecting the source coming from
	the same direction.

	The angle given starts from the line in which the front mic lies. A negative angle is to the
	left; a positive angle is to the right, all the way to the back of the front mic.
	
	If more than source is active, because the sample size is small, it is likely that one sample is able
	to get through redundancy check and provide an adequate angle of one source, and the next redundant
	angle of another source. For this purpose, an expectation system is built on top of the DOA estimator
	to mantain in memory the number of sources and each of their angles, which provides a basic CASA system.

	When running, it will output to standard out its estimate of the angle, and the current angle of the robot
	simulation until it reaches a good enough angle.
	
	Requirements:
		- A soundcard with at least three channels of input, and set to be the system's default.
		- A running JACK Server.
*/
#ifndef _LIBMULTISOUNDLOC_H_
#define _LIBMULTISOUNDLOC_H_


#include "max.h"
#include "unwrap.h"
#include "phat.h"
#include "angleManipulation.h"
#include "kmeans.h"
#include "kalman.h"


#define TDE_CC 1
#define TDE_GCCPHAT 2
#define TRACK_CLUSTER 1
#define TRACK_MCMCDA 2
#define TRACK_KALMAN 3
#define INVALID 999999
#define PI 3.141592653589793238462643383279502884

#ifndef M_PI
#define M_PI 3.1415926535897932384
#endif

#ifndef ROUND
#define ROUND(a) (PLINT)((a)<0. ? ((a)-0.5) : ((a)+0.5))
#endif

#ifdef PL_USE_NAMESPACE
using namespace std;
#endif

extern int SOUNDLOC_CAPTURING;
extern int SOUNDLOC_PAUSE;

// User description and registry variable
typedef struct{
	double value;
	time_t reg_time;
} doa;

typedef struct{
	double main_doa;
	double last_doa;
	time_t main_doa_time;
	std::vector < doa > doas;
	float confidence;
	bool real_source;
} source;

typedef struct{
	std::vector < double > doas;
	time_t reg_time;
	char audio_file_path[100];
} vad_results;

extern std::vector < source > sources;
extern std::vector< doa > doas;
extern vad_results vad_res;

extern void millisleep(int milli);
extern void soundloc_finalize();
extern void soundloc_clear();
extern void soundloc_init(double distance_between_mics_in, int number_of_samples_of_noise_in = 10, double noise_threshold_in = 0.002, double noise_peak_change_in = 0.004, double coherence_threshold_in = 30, int tde_method = TDE_GCCPHAT, int track_method = TRACK_CLUSTER, bool reverse_angle_in = true, bool graph_out_in = true, bool connect_ports = true, bool write_vad_wavs_in = true);


#endif
