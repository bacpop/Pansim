extern crate rand;
extern crate statrs;
extern crate rayon;

use rayon::prelude::*;

use statrs::distribution::Poisson;

use rand::distributions::WeightedIndex;
use rand::Rng;
use rand::seq::IteratorRandom;
use crate::rand::distributions::Distribution;

struct Population {
    pop: Vec<Vec<i8>>,
    core : bool,
    core_vec : Vec<Vec<i8>>,
}

impl Population {
    fn new(size: usize, allele_count: usize, max_variants: i8, core : bool) -> Self {
        let mut pop: Vec<Vec<i8>> = vec![vec![0; allele_count]; size]; // each row is individual, each column is an allele

        let mut start: Vec<i8> = vec![0; allele_count];

        for j in 0..allele_count {
            start[j] = rand::thread_rng().gen_range(0..max_variants);
        }

        for i in 0..size {
            pop[i] = start.clone();
        }

        let core_vec: Vec<Vec<i8>> = vec![vec![1, 2, 3],
                                          vec![0, 2, 3],
                                          vec![0, 1, 3],
                                          vec![0, 1, 2]];

        Self {
            pop,
            core,
            core_vec,
        }
    }

    fn next_generation(&mut self, sample : &Vec<usize>) {

        let next_pop: Vec<Vec<i8>> = sample.iter().map(|&i| self.pop[i].clone()).collect();

        self.pop = next_pop;
    }

    fn mutate_alleles(&mut self, weights : &Vec<i32>, mutations : i32, n_threads : usize) {

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)  // Set the number of threads
            .build()
            .unwrap();

        let weighted_dist = WeightedIndex::new(weights).unwrap();
        if self.core != false {
            pool.install(|| {
                self.pop.par_iter_mut().for_each(|row| {
                    // generate Poisson sampler
                    let mut rng = rand::thread_rng();
                    let poisson = Poisson::new(mutations as f64).unwrap();

                    // sample from Poisson distribution for number of sites to mutate in this isolate
                    let n_sites = rng.sample(poisson) as usize;

                    // iterate for number of mutations required to reach mutation rate
                    for _ in 0..n_sites {
                        // sample new site to mutate
                        let mutant_site = weighted_dist.sample(&mut rng);
                        let value = row.get_mut(mutant_site).unwrap();
                        //let value = self.pop[i][mutant_site];
                        let new_allele : i8 = if *value == 0 as i8 { 1 } else { 0 };

                        // set value in place
                        *value = new_allele;
                    }
                });
            });
        } else {
            pool.install(|| {
                self.pop.par_iter_mut().for_each(|row| {
                    // generate Poisson sampler
                    let mut rng = rand::thread_rng();
                    let poisson = Poisson::new(mutations as f64).unwrap();

                    // sample from Poisson distribution for number of sites to mutate in this isolate
                    let n_sites = rng.sample(poisson) as usize;

                    // iterate for number of mutations required to reach mutation rate
                    for _ in 0..n_sites {
                        // sample new site to mutate
                        let mutant_site = weighted_dist.sample(&mut rng);

                        // get possible values to mutate to, must be different from current value
                        let value = row.get_mut(mutant_site).unwrap();
                        let values = &self.core_vec[*value as usize];

                        // sample new allele
                        let new_allele = values.iter().choose_multiple(&mut rng, 1)[0];

                        // set value in place
                        *value = *new_allele;
                    }
                });
            });
        }
    }
}

fn main() {
    let n_threads = 4;

    let pop_size = 100;
    let core_size = 1000000;
    let pan_size = 5000;
    let n_gen = 100;

    // core and pangenome mutation rates
    let core_mu = 0.05;
    let pan_mu = 0.05;

    // calculate number of mutations per genome per generation
    let n_core_mutations = ((core_size as f64 * core_mu) / n_gen as f64) / 2.0 ;
    let n_pan_mutations = ((pan_size as f64 * pan_mu) / n_gen as f64) / 2.0;

    let core_weights : Vec<i32> = vec![1; core_size];
    let pan_weights : Vec<i32> = vec![1; pan_size];

    let mut core_genome = Population::new(pop_size, core_size, 4, true); // core genome alignment
    let mut pan_genome = Population::new(pop_size, pan_size, 2, false); // pangenome alignment

    let mut rng = rand::thread_rng();

    for j in 0..n_gen { // Run for n_gen generations
        // mutate core genome
        println!("started {}", j);
        core_genome.mutate_alleles(&core_weights, n_core_mutations as i32, n_threads);
        println!("finished mutating core genome {}", j);
        pan_genome.mutate_alleles(&pan_weights, n_pan_mutations as i32, n_threads);
        println!("finished mutating pangenome {}", j);



        // sample new individuals
        let sampled_individuals: Vec<usize> = (0..pop_size).map(|_| rng.gen_range(0..pop_size)).collect();
        core_genome.next_generation(& sampled_individuals);
        println!("finished copying core genome {}", j);
        pan_genome.next_generation(& sampled_individuals);
        println!("finished copying pangenome {}", j);
    }


}