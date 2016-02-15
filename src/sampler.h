#ifndef SAMPLER_H
#define SAMPLER_H

#include <vector>
#include <math.h>

#include "math/point.h"

class Sampler 
{	
	public:	
		
		Sampler(void);								
		
		Sampler(const int num);						
		
		Sampler(const int num, const int num_sets);		

		Sampler(const Sampler& s);						

		Sampler& operator= (const Sampler& rhs);				

		virtual Sampler*clone() const = 0;			

		virtual ~Sampler();							
		
		void set_num_sets(const int np);					
		
		virtual void generate_samples() = 0;

		int get_num_samples();							
		
		void shuffle_x_coordinates();
		
		void shuffle_y_coordinates();
			
		void setup_shuffled_indices();
		
		void map_samples_to_unit_disk();
		
		void map_samples_to_hemisphere(const float p);		
		
		void map_samples_to_sphere();					
			
		Point sample_unit_square();							
		
		Point sample_unit_disk();							
		
		Point sample_hemisphere();								
		
		Point sample_sphere();										
	
		Point sample_one_set();												
		
	protected:
	
		int 					num_samples;     		// the number of sample points in a set
		int 					num_sets;				// the number of sample sets
		std::vector<Point>		samples;				// sample points on a unit square
		std::vector<int>		shuffled_indices;		// shuffled samples array indices
		std::vector<Point>		disk_samples;			// sample points on a unit disk
		std::vector<Point> 		hemisphere_samples;		// sample points on a unit hemisphere
		std::vector<Point> 		sphere_samples;			// sample points on a unit sphere
		unsigned long 			count;					// the current number of sample points used
		int 					jump;					// random index jump
};

#endif

