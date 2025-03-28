#include "spark/particle/pusher.h"
#include "spark/particle/species.h"
#include <thread>

template <>
void spark::particle::move_particles(spark::particle::ChargedSpecies<1, 3>& species,
                                     const core::TMatrix<core::Vec<1>, 1>& force,
                                     const double dt) {
    const size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();
    const auto* f = force.data_ptr();
    const double k = species.q() * dt / species.m();
   
    auto thread_fn = [&](size_t inicio, size_t fim){
	    for(int i = inicio; i < fim ; i++){
	 	v[i].x += f[i].x * k;
	        x[i].x += v[i].x * dt;   
	    }
    };
    int n_threads = 10;
    std::thread th[n_threads];
    int blockD = n/n_threads;
    for(int i = 0; i < n_threads; i++){
    	th[i] = std::thread(thread_fn, i*blockD, ((i*blockD)+blockD)-1);
    }
    for(int i = 0; i < n_threads; i++){
    	th[i].join();
    }
    if(n%n_threads != 0){
    	std::thread thf(thread_fn, (blockD*n_threads)+1, (blockD*n_threads)+1+(n%n_threads));
	thf.join();
    }
}

template <>
void spark::particle::move_particles(spark::particle::ChargedSpecies<2, 3>& species,
                                     const core::TMatrix<core::Vec<2>, 1>& force,
                                     const double dt) {
    const size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();
    const auto* f = force.data_ptr();
    const double k = species.q() * dt / species.m();
    for (size_t i = 0; i < n; i++) {
        v[i].x += f[i].x * k;
        v[i].y += f[i].y * k;

        x[i].x += v[i].x * dt;
        x[i].y += v[i].y * dt;
    }
}
