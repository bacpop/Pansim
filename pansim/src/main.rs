extern crate rand;
extern crate statrs;
extern crate rayon;
extern crate ndarray;

use rayon::prelude::*;

use statrs::distribution::Poisson;

use rand::rngs::StdRng;
use rand::distributions::WeightedIndex;
use rand::{Rng, SeedableRng};
use rand::seq::IteratorRandom;
use crate::rand::distributions::Distribution;

use ndarray::{Array1, Array2, Axis, s};
use ndarray::parallel::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

fn hamming_distance(x: &[u8], y: &[u8]) -> u64 {
    assert_eq!(x.len(), y.len(), "Vectors must have the same length");
    x.iter().zip(y).fold(0, |a, (b, c)| a + (*b ^ *c).count_ones() as u64)
}

fn jaccard_distance(row1: &[u8], row2:  &[u8]) -> (usize, usize) {
    assert_eq!(row1.len(), row2.len(), "Rows must have the same length");

    let intersection: usize = row1.iter().zip(row2.iter()).filter(|&(x, y)| *x == 1 && *y == 1).count();
    let union: usize = row1.iter().zip(row2.iter()).filter(|&(x, y)| *x == 1 || *y == 1).count();

    (intersection, union)
}


struct Population {
    pop: Array2<u8>,
    core : bool,
    core_vec : Vec<Vec<u8>>,
}

fn to_array2<T: Copy>(source: Vec<Array1<T>>) -> Result<Array2<T>, impl std::error::Error> {
    let width = source.len();
    let flattened: Array1<T> = source.into_iter().flat_map(|row| row.to_vec()).collect();
    let height = flattened.len() / width;
    flattened.into_shape((width, height))
}

impl Population {
    fn new(size: usize, allele_count: usize, max_variants: u8, core : bool) -> Self {
        //let mut pop = Array2::<u8>::zeros((size, allele_count));

        let mut start: Array1<u8> = Array1::zeros(allele_count);

        for j in 0..allele_count {
            start[j] = rand::thread_rng().gen_range(0..max_variants);
        }

        // Create a vector of Array1 filled with the same values
        let pop_vec: Vec<Array1<u8>> = std::iter::repeat(start)
            .take(size)
            .collect();

        // convert vector into 2D array
        let pop = to_array2(pop_vec).unwrap();

        let core_vec: Vec<Vec<u8>> = vec![vec![1, 2, 3],
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

        // Create a new Array2 by sampling rows from the original Array2
        let sampled_array = sample
            .iter()
            .map(|&i| self.pop.slice(s![i, ..]))
            .collect::<Vec<_>>();
        
        // Concatenate the sampled rows into a new Array2
        let next_pop = ndarray::stack(Axis(0), &sampled_array).unwrap();

        self.pop = next_pop;
    }

    fn mutate_alleles(&mut self, weights : &Vec<i32>, mutations : i32, rng : &StdRng) {

        let weighted_dist = WeightedIndex::new(weights).unwrap();
        let index = AtomicUsize::new(0);


        if self.core == false {
            self.pop.axis_iter_mut(Axis(0)).into_par_iter().for_each(|mut row| {
                // thread-specific random number generator
                let mut thread_rng = rng.clone();
                let current_index = index.fetch_add(1, Ordering::SeqCst);

                // Jump the state of the generator for this thread
                for _ in 0..current_index {
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
                    let value = row[mutant_site];
                    let new_allele : u8 = if value == 0 as u8 { 1 } else { 0 };


                    // set value in place
                    row[mutant_site] = new_allele;
                }
            }  
            );
        } else {
            self.pop.axis_iter_mut(Axis(0)).into_par_iter().for_each(|mut row| {
                    // thread-specific random number generator
                    let mut thread_rng = rng.clone();
                    let current_index = index.fetch_add(1, Ordering::SeqCst);

                    // Jump the state of the generator for this thread
                    for _ in 0..current_index {
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
                        let value = row[mutant_site];
                        let values = &self.core_vec[value as usize];

                        // sample new allele
                        let new_allele = values.iter().choose_multiple(&mut thread_rng, 1)[0];

                        // set value in place
                        row[mutant_site] = *new_allele;
                    }
                });
        }
    }

    fn pairwise_distances(&mut self, max_distances : usize, rng : &mut StdRng) -> Vec<f64> {
        let n_rows = self.pop.nrows();
    
        // determine which columns are all equal, ignore from distance calculations
        let array_f64 = self.pop.mapv(|x| x as f64);
        let column_variance = array_f64.var_axis(Axis(0), 0.0);

        // Determine which column indices have variance greater than 0
        let columns_to_iter: Vec<usize> = column_variance
            .iter()
            .enumerate()
            .filter_map(|(i, &variance)| if variance > 0.0 { Some(i) } else { None })
            .collect();
        
        let matches = column_variance.len() as f64 - columns_to_iter.len() as f64;

        let subset_array: Array2<u8> = self.pop.select(Axis(1), &columns_to_iter).to_owned().reversed_axes();
        //println!("{:?}", self.pop);
        //println!("{:?}", subset_array);
    
        //let mut distances : Vec<f64> = vec![0.0; n_rows * (n_rows - 1) / 2]; // Capacity for pairwise combinations
        
        // sample distances with replacement
        let mut distances : Vec<f64> = vec![0.0; max_distances as usize]; // Capacity for pairwise combinations

        //let mut idx = 0;
        for idx in 0..max_distances {
            let i = rng.gen_range(0..n_rows) as usize;
            let j = rng.gen_range(i + 1..n_rows + 1) as usize;
            
            let row1 = subset_array.index_axis(Axis(0), i);
            let row2 = subset_array.index_axis(Axis(0), j);

            if self.core == true {
                let distance = hamming_distance(row1.as_slice().unwrap(), &row2.as_slice().unwrap());
                let hamming_distance = distance as f64 / (column_variance.len() as f64);
                distances[idx] = hamming_distance;
            } else {
                let (intersection, union) = jaccard_distance(&row1.as_slice().unwrap(), &row2.as_slice().unwrap());
                let jaccard_distance = 1.0 - ((intersection as f64 + matches) / (union as f64 + matches));
                distances[idx] = jaccard_distance;
            }
        }
        distances
        }
}

fn main() {
    use std::time::Instant;
    let now = Instant::now();
    
    let n_threads = 4;
    rayon::ThreadPoolBuilder::new().num_threads(n_threads).build_global().unwrap();

    let pop_size = 1000;
    let core_size = 1200000;
    let pan_size = 6000;
    let n_gen = 100;

    let max_distances: usize = 100000;

    // core and pangenome mutation rates
    let core_mu = 0.05;
    let pan_mu = 0.05;

    // calculate number of mutations per genome per generation
    let n_core_mutations = (((core_size as f64 * core_mu) / n_gen as f64) / 2.0).ceil() ;
    let n_pan_mutations = (((pan_size as f64 * pan_mu) / n_gen as f64) / 2.0).ceil();

    let core_weights : Vec<i32> = vec![1; core_size];
    let pan_weights : Vec<i32> = vec![1; pan_size];

    let mut core_genome = Population::new(pop_size, core_size, 4, true); // core genome alignment
    let mut pan_genome = Population::new(pop_size, pan_size, 2, false); // pangenome alignment

    let seed: u64 = 0;
    let mut rng: StdRng = StdRng::seed_from_u64(seed);

    for j in 0..n_gen { // Run for n_gen generations
        let now_gen = Instant::now();
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

            core_genome.mutate_alleles(&core_weights, n_core_mutations as i32, &rng);

            // Jump the state of the generator
            for _ in 0..pop_size {
                rng.gen::<u64>(); // Discard some numbers to mimic jumping
            }

            //println!("finished mutating core genome {}", j);
            pan_genome.mutate_alleles(&pan_weights, n_pan_mutations as i32, &rng);
            //println!("finished mutating pangenome {}", j);

            // Jump the state of the generator
            for _ in 0..pop_size {
                rng.gen::<u64>(); // Discard some numbers to mimic jumping
            }
        } else {
            // else calculate hamming and jaccard distances
            let core_distances = core_genome.pairwise_distances(max_distances, &mut rng);
            let acc_distances = pan_genome.pairwise_distances(max_distances, &mut rng);
        }

        let elapsed = now_gen.elapsed();
        println!("Finished gen: {}", j);
        println!("Elapsed: {:.2?}", elapsed);
    }
    let elapsed = now.elapsed();
    println!("Total elapsed: {:.2?}", elapsed);


}