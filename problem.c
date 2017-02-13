/**
 * Solar System
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "rebound.h"
#include "reboundx.h"
#include <time.h>


double ss_pos[10][3] = 
{

	{-7.139143380212697E-03, -2.792019770161695E-03,  2.061838852554664E-04}, // Sun
	{-1.478672233442572E-01, -4.466929775364947E-01, -2.313937582786785E-02}, // Mercury
	{-7.257693602841776E-01, -2.529582082587794E-02,  4.137802526208009E-02}, // Venus
	{-1.756637922977122E-01,  9.659912850526895E-01,  2.020629118443605E-04}, // Earth
	{1.383221922520998E+00, -2.380174081741852E-02, -3.441183028447500E-02}, // Mars
	{3.996321310086093E+00,  2.932561197358908E+00, -1.016170544300634E-01}, // Jupiter
	{6.401416890663500E+00,  6.565250734685104E+00, -3.689211141720000E-01}, // Saturn
	{1.442337843936191E+01, -1.373845030140273E+01, -2.379221201389048E-01}, // Uranus
	{1.680361764335730E+01, -2.499544328458694E+01,  1.274772016011350E-01}, // Neptune
	{-9.884006595855096E+00, -2.796081320603301E+01,  5.851020841275061E+00}, // Pluto
};

double ss_vel[10][3] = 
{
	{5.374260940168566E-06, -7.410965396701423E-06, -9.422862838391440E-08},
	{2.117424563261189E-02, -7.105386404267509E-03, -2.522925180072137E-03},	
	{5.189070188671265E-04, -2.031355258779472E-02, -3.072687386494688E-04},
 	{-1.722857156974862E-02, -3.015071224668472E-03, -5.859931223618532E-08},
	{7.533013850513376E-04,  1.517888771209419E-02,  2.996589710207392E-04},
	{-4.558376590671486E-03,  6.439863246141724E-03,  7.537593486203765E-05},
	{-4.285166238539475E-03,  3.884579926659154E-03,  1.025155282571916E-04},
	{2.683840344076701E-03,  2.665016541217002E-03, -2.484232267336756E-05},
	{2.584589572083709E-03,  1.768943546348827E-03, -9.629380362804233E-05},
	{3.044390641302692E-03, -1.537290075737863E-03, -7.173359336658472E-04},
};

double ss_mass[10] =
{
	1.988544e30, // Sun
	3.302e23, // Mercury
	48.685e23, // Venus
	5.97219e24, // Earth
	6.4185e23, // Mars
	1898.13e24, // Jupiter
	5.68319e26, // Saturn
	86.8103e24, // Uranus
	102.41e24, // Neptune
	1.307e22, // Pluto
};


void heartbeat(struct reb_simulation* r);
int main(int argc, char* argv[]) {

	double cm_to_au = 6.68459e-14;

	srand ( time(NULL) );
	double random_number=((double)rand()/(double)RAND_MAX);

//	double kappa = 0.75 + random_number*0.45;
	double kappa = 1.0;

	double k = 1.0;

	char filename[512];
	char rebx_filename[512];

	time_t t = time(NULL);
	struct tm tm = *localtime(&t);

	sprintf(filename,"/data_local/arachkov/solar_system/kappas/kappa_%e_date_%d_%d_%d_time_%dh_%dmin_%dsec.bin",kappa,tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
	sprintf(rebx_filename,"/data_local/arachkov/solar_system/kappas/rebx_date_%d_%d_%d_time_%dh_%dmin_%dsec.bin",tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);

	// Trying to restart from the Simulation Archive.
	struct reb_simulation* r = reb_create_simulation_from_simulationarchive(filename);

	// Check if that was successful
	if (r==NULL){
		printf("No simulation archive found. Creating new simulation.\n");
		r= reb_create_simulation();
		r->G			= 1.4880826e-34;		// in AU^3 / kg / day^2.

		// Initial conditions
		for (int i=0;i<10;i++){
			struct reb_particle p = {0};
			p.x  = ss_pos[i][0]; 		p.y  = ss_pos[i][1];	 	p.z  = ss_pos[i][2];
			p.vx = ss_vel[i][0]; 		p.vy = ss_vel[i][1];	 	p.vz = ss_vel[i][2];
			p.m  = ss_mass[i];
			reb_add(r, p); 
		}

		reb_move_to_com(r);
		r->dt = 8.;
		r->simulationarchive_interval = 1.e3; // output data every 100000 years
		r->ri_whfast.safe_mode = 0;                      
		r->ri_whfast.corrector = 5;    
		r->integrator = REB_INTEGRATOR_WHFAST;
		r->heartbeat  = heartbeat;
	}else{
		printf("Found simulation archive. Loaded snapshot at t=%.16f.\n",r->t);
	}

	struct rebx_extras* rebx = rebx_create_extras_from_binary(r, rebx_filename);
	if (rebx==NULL){
		printf("No reboundx simulation archive found. Creating new rebx effects.\n");
		rebx = rebx_init(r);
		struct rebx_effect* gr_params = rebx_add(rebx, "gr_potential");
		double* c = rebx_add_param(gr_params, "c", REBX_TYPE_DOUBLE);
		*c = REBX_C;
		rebx_output_binary(rebx, rebx_filename);
	}else{
		printf("Found rebx effects archive.");
	}

	r->simulationarchive_filename = filename;
//	reb_integrate(r, 2.*M_PI*5e9); //time
	reb_integrate(r, 5e2); //time
	rebx_free(rebx);
	reb_free_simulation(r);

}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 1.e3*2.*M_PI)){
        reb_output_timing(sim, 2.*M_PI*5e9);
    }
}
