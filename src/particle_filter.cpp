/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    //number of particles
    num_particles = 3;
    
    //create normal Gaussian distribution for x, y and theta
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    
    //initialize all particles
    for(int i = 0; i < num_particles; ++i){
        
        //declare position values
        double sample_x, sample_y, sample_theta;
        
        //initialize a particle to gaussian distribution around first position
        sample_x = dist_x(gen);
        sample_y = dist_y(gen);
        sample_theta = dist_theta(gen);
        
        //declare single particle
        Particle particle_temp;
        
        //set values
        particle_temp.id = i;
        particle_temp.x = sample_x;
        particle_temp.y = sample_y;
        particle_temp.theta = sample_theta;
        particle_temp.weight = 1.0;
        
        //add to list of particles
        particles.push_back(particle_temp);
        weights.push_back(particle_temp.weight);
    }
    
    //set isinitialized true
    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    
    //create normal distribution for x, y, and yaw_rate
    normal_distribution<double> dist_x(0.0, std_pos[0]);
    normal_distribution<double> dist_y(0.0, std_pos[1]);
    normal_distribution<double> dist_yaw_rate(0.0, std_pos[2]);
    
    //Predict for each particle
    for(int i = 0; i < num_particles; ++i){
        
        //when yaw_rate is zero or very close to zero
        //if(fabs(yaw_rate) < 1e-3){
        if(fabs(yaw_rate) == 0.0){
            
            //bicycle motion model when yaw rate is zero
            particles[i].x = particles[i].x + dist_x(gen) + velocity * delta_t * cos(particles[i].theta);
            particles[i].y = particles[i].y + dist_y(gen) + velocity * delta_t * sin(particles[i].theta);
            
        }
        //when yaw_rate is not zero
        else{
            
            //generate a random yaw rate adding random guassian noise
            double yaw_rate_temp = dist_yaw_rate(gen);
            double yaw_temp = particles[i].theta + yaw_rate_temp * delta_t;
            
            //bicycle motion model when yaw rate is not zero
            particles[i].x = particles[i].x + dist_x(gen) + velocity * (sin(yaw_temp) - sin(particles[i].theta)) / yaw_rate_temp;
            particles[i].y = particles[i].y + dist_y(gen) + velocity * (cos(particles[i].theta) - cos(yaw_temp)) / yaw_rate_temp;
            particles[i].theta = yaw_temp;
        }
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    
    //observations size
    int num_observations = observations.size();
    
    //landmark size
    int num_landmarks = map_landmarks.landmark_list.size();
    
    //calculate normalization term
    double gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
    
    //data association and calculate weight for each particle
    for(int i = 0; i < num_particles; ++i){
        
        //values of particle location
        double x_temp = particles[i].x;
        double y_temp = particles[i].y;
        double theta_temp =  particles[i].theta;
        double weight_temp = 1.0;
        
        //transform coordinates and associations of each observation
        for(int j = 0; j < num_observations; ++j){
            
            //values of single observation
            double observation_x_temp = observations[j].x;
            double observation_y_temp = observations[j].y;
            
            //defind a single particle observation
            LandmarkObs particle_obsevn;
            
            //perform transformation
            particle_obsevn.x = x_temp + (cos(theta_temp) * observation_x_temp) - (sin(theta_temp) * observation_y_temp);
            particle_obsevn.y = y_temp + (sin(theta_temp) * observation_x_temp) + (cos(theta_temp) * observation_y_temp);
            
            std::cout << "trans_observation_x: " << particle_obsevn.x << "\n";
            std::cout << "trans_observation_y: " << particle_obsevn.y << "\n";
            
            //initialize minimum distance between observation and landmark
            double min_dist_obsevn_landmark = sensor_range;
            
            for(int k = 0; k < num_landmarks; ++k){
                
                double dist_particle_landmark = dist(x_temp,y_temp,map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
                
                //if ((fabs(map_landmarks.landmark_list[k].x_f - x_temp) <= sensor_range)
                //    && (fabs(map_landmarks.landmark_list[k].y_f - y_temp) <= sensor_range)) {
                if (dist_particle_landmark <= sensor_range) {
                    
                    double dist_obsevn_landmark = dist(particle_obsevn.x, particle_obsevn.y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
                    
                    //find the true associated landmark by nearest neighbor
                    if (dist_obsevn_landmark < min_dist_obsevn_landmark){
                        
                        particle_obsevn.id = (map_landmarks.landmark_list[k].id_i - 1);
                        min_dist_obsevn_landmark = dist_obsevn_landmark;
                        
                    }
                }
            }
            
            std::cout << "asso_observation_x: " << map_landmarks.landmark_list[particle_obsevn.id].x_f << "\n";
            std::cout << "asso_observation_y: " << map_landmarks.landmark_list[particle_obsevn.id].y_f << "\n";
            
            //calculate exponent
            double exponent_x = pow((map_landmarks.landmark_list[particle_obsevn.id].x_f - particle_obsevn.x),2) / (2 * pow(std_landmark[0],2));
            double exponent_y = pow((map_landmarks.landmark_list[particle_obsevn.id].y_f - particle_obsevn.y),2) / (2 * pow(std_landmark[1],2));
            double exponent = exponent_x + exponent_y;
            
            std::cout << "exponent value: " << exponent << "\n";
            
            //calculate weight using normalization term and exponent
            weight_temp *= (gauss_norm * exp(-exponent));
            std::cout << "weight " << j << ":" << weight_temp << "\n";
        }
        
        //update particle weight
        particles[i].weight = weight_temp;
        weights[i] = weight_temp;
        std::cout << "weight " << i << ":" << weights[i] << "\n";
    }
    
    //weights normalization
    //sum of all particles weights
    double sum_of_weights = accumulate(weights.begin(), weights.end(), 0.0);
    
    for (int i = 0; i < num_particles; ++i){
        
        weights[i] = weights[i] / sum_of_weights;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    //create a discrete distribution of all weights
    std::discrete_distribution<int> d_dist(std::begin(weights), std::end(weights));
    
    //create a temp vector of resampled particles;
    std::vector<Particle> particles_temp;
    
    for(int i=0; i < num_particles; ++i){
        
        particles_temp.push_back(particles[d_dist(gen)]);
    }
    
    particles = particles_temp;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
