classdef pf < handle
    %PF Particle Filter
    
    properties
        particles
        weights
        ci_particles
        ci_weights
        noise_particles
        noise_weights
        system
    end
    
    methods
        function obj = pf(system)
            %PF Construct an instance of this class
            obj.system = system;
        end
        
        function [particles,weights] = sampleGaussianParticles(obj, x0, cov, particle_count)
            %sampleParticles Creates particle set from gaussian
            % distribution
            particles = mvnrnd(x0, cov, particle_count)';
            weights = (1/particle_count).*ones(1,size(particles,2));
        end
        
        function setParticles(obj, particles, weights)
            obj.particles = particles;
            obj.weights = weights;
        end
        
        function setCommonInformation(obj, particles, weights)
            obj.ci_particles = particles;
            obj.ci_weights = weights;
        end
        
        function x_f = predict(obj,dT,u)
            % noise particles
            noise_samples = obj.system.systemNoiseSample(obj.particles,dT);

            % state particles
            obj.particles = obj.system.predict(obj.particles,dT,u)+noise_samples;
            x_f = sum(obj.particles.*obj.weights,2);
                        
            % common information
            if size(obj.ci_particles,2) > 0
                if size(obj.ci_particles,2) == size(noise_samples,2)
                    obj.ci_particles = obj.system.predict(obj.ci_particles,dT,u)+noise_samples;
                else
                    obj.ci_particles = obj.system.predict(obj.ci_particles,dT,u)+obj.system.systemNoiseSample(obj.ci_particles,dT);
                end    
            end
        end
        
        function x_e = update(obj,z)
            likelihoods = obj.system.likelihood(obj.particles,z);            
            obj.weights = obj.weights.*likelihoods;            
            obj.weights = obj.weights ./ sum(obj.weights);
                        
            % state particles
            x_e = sum(obj.particles.*obj.weights,2);        
            
            obj.resampleWhenNecessaire();
        end
        
        function x = getFilterMean(obj)
            x = sum(obj.particles.*obj.weights,2);
        end
        
        function resampleWhenNecessaire(obj)
            % resample
            Ns = size(obj.weights, 2);
            deg = 1 / sum( obj.weights.^2 );
            if deg < Ns / 2
                obj.resample();
            end 
        end
        
        function resample(obj)
            idx = sort( resampleMultinomial( obj.weights ) );
            obj.particles = obj.particles(:,idx);
            obj.weights = ones(size(obj.weights)).*(1/size(obj.weights,2));
        end
        
    end
end

