#include "spark/particle/pusher.h"
#include "spark/threads/pool.h"
#include "spark/particle/species.h"
#include <thread>

template <>
void spark::particle::move_particles(spark::particle::ChargedSpecies<1, 3>& species,
                                     const core::TMatrix<core::Vec<1>, 1>& force,
                                     const double dt,
                                     spark::threads::ThPool& pool) {
    const size_t n = species.n();
    auto* v = species.v();
    auto* x = species.x();
    const auto* f = force.data_ptr();
    const double k = species.q() * dt / species.m();

    // for (size_t i = 0; i < n; i++) {
    //     v[i].x += f[i].x * k;
    //     x[i].x += v[i].x * dt;
    // }
    size_t n_threads = pool.get_num_th();
    int block_id = n/n_threads;
    for(size_t i = 0; i < n_threads; i++){
        pool.enfileira([=](size_t inicio, size_t fim){
	        for(int i = inicio; i < fim ; i++){
	 	        v[i].x += f[i].x * k;
	            x[i].x += v[i].x * dt;
	        }
        }, i*block_id, ((i*block_id)+block_id)-1);

    }

}

template <>
void spark::particle::move_particles(spark::particle::ChargedSpecies<2, 3>& species,
                                     const core::TMatrix<core::Vec<2>, 1>& force,
                                     const double dt,
                                     spark::threads::ThPool& pool) {
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
