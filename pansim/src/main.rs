extern crate rand;
extern crate statrs;

use statrs::distribution::{Poisson, Discrete};
use statrs::statistics::Distribution;
use statrs::prec;

use rand::distributions::WeightedIndex;
use rand::prelude::*;
use rand::Rng;
use rand::{seq::IteratorRandom, thread_rng};

fn sample_with_weighted_replacement(weights: &Vec<f64>) -> usize {
    let dist = WeightedIndex::new(weights).unwrap();
    let mut rng = thread_rng();
    dist.sample(&mut rng)
}

fn normalize_vector(vec: &mut Vec<f64>) {
    let sum: f64 = vec.iter().sum();

    for val in vec.iter_mut() {
        *val /= sum;
    }
}

struct Population {
    size: usize,
    allele_count: usize,
    pop: Vec<Vec<i32>>,
}

impl Population {
    fn new(size: usize, allele_count: usize, max_variants: i32) -> Self {
        let mut pop: Vec<Vec<i32>> = vec![vec![0; allele_count]; size]; // each row is individual, each column is an allele

        for i in 0..size {
            for j in 0..allele_count {
                pop[i][j] = rand::thread_rng().gen_range(0..=(max_variants - 1));
            }
        }

        Self {
            size,
            allele_count,
            pop,
        }
    }

    fn next_generation(&mut self, sample : &Vec<usize>) {
        
        let mut next_pop: Vec<Vec<i32>> = sample.iter().map(|&i| self.pop[i].clone()).collect();

        self.pop = next_pop;
    }

    fn mutate_alleles(&self, weights : Vec<f64>, mutations : i32, max_variants : i32) {
        let mut rng = rand::thread_rng();
        
        // iterate over each genome
        for i in 0..self.size {
            // sample from Poisson distribution for number of sites to mutate in this isolate
            let poisson = Poisson::new(mutations as f64).unwrap();
            let n_sites = rng.sample(poisson) as usize;
            // iterate for number of mutations required to reach mutation rate
            for j in 0..n_sites {
                // sample new site to mutate
                let mutant_site : usize = sample_with_weighted_replacement(&weights);
                
                // get possible values to mutate to, must be different from current value
                let mut values = Vec::from_iter(0..(max_variants - 1));
                let value = self.pop[i][j];
                values.retain(|&x| x != value);

                // sample new allele
                let mut new_allele = values.iter().choose_multiple(&mut rng, 1)[0];

                // set value in place
                self.pop[i][j] = *new_allele;
            }
        }
    }
}

fn main() {
    let pop_size = 100;
    let core_size = 100;
    let pan_size = 100;
    let n_gen = 10;

    // core and pangenome mutation rates
    let core_mu = 0.01;
    let pan_mu = 0.001;

    // calculate number of mutations per genome per generation
    let n_core_mutations = ((core_size as f64 * core_mu) / n_gen as f64) / 2.0 ;
    let n_pan_mutations = ((pan_size as f64 * pan_mu) / n_gen as f64) / 2.0;

    let mut core_weights : Vec<f64> = vec![1.0; core_size];
    normalize_vector(&mut core_weights);

    let mut pan_weights : Vec<f64> = vec![1.0; core_size];
    normalize_vector(&mut pan_weights);
    
    let mut core_genome = Population::new(pop_size, core_size, 4); // core genome alignment
    let mut pan_genome = Population::new(pop_size, pan_size, 2); // pangenome alignment

    let mut rng = rand::thread_rng();

    for _ in 0..n_gen { // Run for n_gen generations
        // mutate core genome
        core_genome.mutate_alleles(core_weights, n_core_mutations as i32, 4);
        pan_genome.mutate_alleles(pan_weights, n_pan_mutations as i32, 2);

        // sample new individuals
        let sampled_individuals: Vec<usize> = (0..pop_size).map(|_| rng.gen_range(0..pop_size)).collect();
        core_genome.next_generation(& sampled_individuals);
        pan_genome.next_generation(& sampled_individuals);
    }
}