extern crate rand;
extern crate statrs;
extern crate rayon;
extern crate rand_core;
extern crate rand_xoshiro;

use rayon::prelude::*;

use statrs::distribution::Poisson;

use rand::rngs::StdRng;
use rand::distributions::WeightedIndex;
use rand::{Rng, SeedableRng};
use rand::seq::IteratorRandom;
use crate::rand::distributions::Distribution;

fn hamming_distance(row1: &Vec<i8>, row2: &Vec<i8>) -> f64 {
    assert_eq!(row1.len(), row2.len(), "Rows must have the same length");

    let mut distance = 0;
    for (x, y) in row1.iter().zip(row2.iter()) {
        if x != y {
            distance += 1;
        }
    }

    distance as f64 / row1.len() as f64
}

fn pairwise_hamming_distances(matrix: &Vec<Vec<i8>>) -> Vec<f64> {
    let mut distances = Vec::with_capacity(matrix.len() * (matrix.len() - 1) / 2); // Capacity for pairwise combinations

    for i in 0..matrix.len() {
        for j in i + 1..matrix.len() { // Start from i+1 to avoid duplicate pairs and comparing row with itself
            let distance = hamming_distance(&matrix[i], &matrix[j]);
            distances.push(distance);
        }
    }

    distances
}

fn jaccard_distance(row1: &Vec<i8>, row2: &Vec<i8>) -> f64 {
    assert_eq!(row1.len(), row2.len(), "Rows must have the same length");

    let intersection: usize = row1.iter().zip(row2.iter()).filter(|&(x, y)| *x == 1 && *y == 1).count();
    let union: usize = row1.iter().zip(row2.iter()).filter(|&(x, y)| *x == 1 || *y == 1).count();

    1.0 - (intersection as f64 / union as f64)
}

fn pairwise_jaccard_distances(matrix: &Vec<Vec<i8>>) -> Vec<f64> {
    let mut distances = Vec::with_capacity(matrix.len() * (matrix.len() - 1) / 2); // Capacity for pairwise combinations

    for i in 0..matrix.len() {
        for j in i + 1..matrix.len() { // Start from i+1 to avoid duplicate pairs and comparing row with itself
            let distance = jaccard_distance(&matrix[i], &matrix[j]);
            distances.push(distance);
        }
    }

    distances
}

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

    fn mutate_alleles(&mut self, weights : &Vec<i32>, mutations : i32, n_threads : usize, rng : &StdRng) {

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)  // Set the number of threads
            .build()
            .unwrap();

        // Initialize xoshiro RNG

        let weighted_dist = WeightedIndex::new(weights).unwrap();
        if self.core != false {
            pool.install(|| {
                self.pop.par_iter_mut().enumerate().for_each(|(index, row),  |  {
                    // thread-specific random number generator
                    let mut thread_rng = rng.clone();
                    
                    // Jump the state of the generator for this thread
                    for _ in 0..index {
                        thread_rng.gen::<u64>(); // Discard some numbers to mimic jumping
                    }
                    
                    // generate Poisson sampler
                    let poisson = Poisson::new(mutations as f64).unwrap();

                    // sample from Poisson distribution for number of sites to mutate in this isolate
                    let n_sites = poisson.sample(&mut thread_rng) as usize;

                    // iterate for number of mutations required to reach mutation rate
                    for _ in 0..n_sites {
                        // sample new site to mutate
                        let mutant_site = weighted_dist.sample(&mut thread_rng);
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
                self.pop.par_iter_mut().enumerate().for_each(|(index, row),  |  {
                    // thread-specific random number generator
                    let mut thread_rng = rng.clone();
                    
                    // Jump the state of the generator for this thread
                    for _ in 0..index {
                        thread_rng.gen::<u64>(); // Discard some numbers to mimic jumping
                    }
                    
                    // generate Poisson sampler
                    let poisson = Poisson::new(mutations as f64).unwrap();

                    // sample from Poisson distribution for number of sites to mutate in this isolate
                    let n_sites = thread_rng.sample(poisson) as usize;

                    // iterate for number of mutations required to reach mutation rate
                    for _ in 0..n_sites {
                        // sample new site to mutate
                        let mutant_site = weighted_dist.sample(&mut thread_rng);

                        // get possible values to mutate to, must be different from current value
                        let value = row.get_mut(mutant_site).unwrap();
                        let values = &self.core_vec[*value as usize];

                        // sample new allele
                        let new_allele = values.iter().choose_multiple(&mut thread_rng, 1)[0];

                        // set value in place
                        *value = *new_allele;
                    }
                });
            });
        }
    }
}

fn main() {
    use std::time::Instant;
    let now = Instant::now();
    
    let n_threads = 1;

    let pop_size = 100;
    let core_size = 2000000;
    let pan_size = 6000;
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

    let seed: u64 = 0;
    let mut rng: StdRng = StdRng::seed_from_u64(seed);

    for j in 0..n_gen { // Run for n_gen generations
        
        // sample new individuals if not at first generation
        if j > 1 {
            let sampled_individuals: Vec<usize> = (0..pop_size).map(|_| rng.gen_range(0..pop_size)).collect();
            core_genome.next_generation(& sampled_individuals);
            //println!("finished copying core genome {}", j);
            pan_genome.next_generation(& sampled_individuals);
            //println!("finished copying pangenome {}", j);
        }
        
        // if at final generation, just sample, otherwise mutate
        if j < (n_gen - 1) {
                    // mutate core genome
            //println!("started {}", j);

            core_genome.mutate_alleles(&core_weights, n_core_mutations as i32, n_threads, &rng);

            // Jump the state of the generator
            for _ in 0..pop_size {
                rng.gen::<u64>(); // Discard some numbers to mimic jumping
            }

            //println!("finished mutating core genome {}", j);
            pan_genome.mutate_alleles(&pan_weights, n_pan_mutations as i32, n_threads, &rng);
            //println!("finished mutating pangenome {}", j);

            // Jump the state of the generator
            for _ in 0..pop_size {
                rng.gen::<u64>(); // Discard some numbers to mimic jumping
            }
        } else {
            // else calculate hamming and jaccard distances
            let core_distances = pairwise_hamming_distances(&core_genome.pop);
            let acc_distances = pairwise_jaccard_distances(&pan_genome.pop);
        }
    }
    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);


}