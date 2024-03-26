extern crate rand;

use rand::distributions::WeightedIndex;
use rand::Rng;

struct Population {
    size: usize,
    allele_count: usize,
    pop: Vec<Vec<i32>>,
}

fn sample_with_weighted_replacement(weights: &[f64]) -> usize {
    let dist = WeightedIndex::new(weights).unwrap();
    let mut rng = rand::thread_rng();
    dist.sample(&mut rng)
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

    fn next_generation(&mut self, &Vec<usize> sample) {
        
        let mut next_pop: Vec<Vec<i32>> = sample.iter().map(|&i| self.pop[i].clone()).collect();

        self.pop = next_pop;
    }

    fn mutate_alleles(&self, weights : Vec<f64>) {
        // iterate over each genome
        for i in 0..self.size {
            // iterate for number of mutations required to reach mutation rate

        }

        let sampled_indices: Vec<usize> = (0..num_samples)
        .map(|_| sample_with_weighted_replacement(&weights))
        .collect();
    }

    fn sample_allele(&self) -> usize {
        let mut rng = rand::thread_rng();
        let mut cumulative_freq = 0.0;
        let r = rng.gen::<f64>();

        for i in 0..self.allele_count {
            cumulative_freq += self.allele_freq[i];
            if r <= cumulative_freq {
                return i;
            }
        }

        self.allele_count - 1 // Should not reach here
    }
}

fn main() {
    let pop_size = 100;
    let core_size = 100;
    let pan_size = 100;
    let n_gen = 10;
    
    let mut core_genome = Population::new(pop_size, core_size, 4); // core genome alignment
    let mut pan_genome = Population::new(pop_size, pan_size), 2; // pangenome alignment
    

    println!("Initial allele frequencies: {:?}", pop.allele_freq);

    let mut rng = rand::thread_rng();

    for _ in 0..n_gen { // Run for 10 generations

        
        let sampled_individuals: Vec<usize> = (0..pop_size).map(|_| rng.gen_range(0..pop_size)).collect();
        



        pop.next_generation();
        println!("Next generation allele frequencies: {:?}", pop.allele_freq);
    }

    // Sample some alleles from the population
    for _ in 0..10 {
        println!("Sampled allele: {}", pop.sample_allele());
    }
}