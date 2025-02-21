extern crate rand;
extern crate statrs;
extern crate rayon;
extern crate ndarray;
use rand_distr::{Beta};

use rayon::{prelude::*};

use statrs::distribution::Poisson;

use rand::rngs::StdRng;
use rand::distributions::WeightedIndex;
use rand::{Rng, SeedableRng};
use rand::seq::IteratorRandom;
use crate::rand::distributions::Distribution;
use rand::seq::SliceRandom;
use rand::distributions::Uniform;

use ndarray::{Array1, Array2, Axis, s};
use ndarray::Zip;
use std::f64::MIN_POSITIVE;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, RwLock};

use std::fs::File;
use std::io::{self, Write};
use std::usize;
use clap::{Arg, Command};

fn jaccard_distance(row1: &[u8], row2:  &[u8]) -> (usize, usize) {
    assert_eq!(row1.len(), row2.len(), "Rows must have the same length");

    let intersection: usize = row1.iter().zip(row2.iter()).filter(|&(x, y)| *x == 1 && *y == 1).count();
    let union: usize = row1.iter().zip(row2.iter()).filter(|&(x, y)| *x == 1 || *y == 1).count();

    (intersection, union)
}

fn average(numbers: &[f64]) -> f64 {
    numbers.iter().sum::<f64>() as f64 / numbers.len() as f64
}

fn standard_deviation(values: &[f64]) -> (f64, f64) {
   let mean = average(values);

   let sum_of_squares: f64 = values
        .iter()
        .map(|&x| (x - mean).powi(2))
        .sum();

   let variance = sum_of_squares / (values.len() as f64);
   (variance.sqrt(), mean)
}

fn sample_beta(num_samples: usize, rng : &mut StdRng) -> Vec<f64> {
    
    let mut return_vec : Vec<f64> = vec![0.0; num_samples];

    // U-shaped beta distribution
    let beta_dist = Beta::new(0.5, 0.5).unwrap();

    for j in 0..num_samples {
        let mut sample = beta_dist.sample(rng);

        // avoid sampling core with frequency 1.0
        while sample == 1.0 {
            sample = beta_dist.sample(rng);
        }
        return_vec[j] = sample;
    }
    return return_vec;
}

fn non_constant_columns(array: &Array2<u8>) -> Vec<usize> {
    (0..array.ncols())
        .filter(|&col| {
            let mut seen = [false; 256]; // Track encountered values
            let mut unique_count = 0;

            for &value in array.column(col).iter() {
                if !seen[value as usize] {
                    seen[value as usize] = true;
                    unique_count += 1;
                    if unique_count > 1 {
                        return true; // Early exit if more than one unique value
                    }
                }
            }
            false
        })
        .collect()
}

fn get_variable_loci (core: bool, pop: &Array2<u8>) -> (ndarray::ArrayBase<ndarray::OwnedRepr<u8>, ndarray::Dim<[usize; 2]>>, f64) {
    
    // Determine which column indices have variance greater than 0
    let columns_to_iter: Vec<usize> = non_constant_columns(&pop);

    // get matches for jaccard distance calculation
    let mut matches: f64 = 0.0;
    if core != true {
        matches = pop.axis_iter(Axis(1)).filter(|col| col.iter().all(|&x| x == 1)).count() as f64;
    }
    
    let subset_array: Array2<u8> = pop.select(Axis(1), &columns_to_iter);
    let mut contiguous_array: ndarray::ArrayBase<ndarray::OwnedRepr<u8>, ndarray::Dim<[usize; 2]>> = Array2::zeros((subset_array.dim().0, subset_array.dim().1));
    contiguous_array.assign(&subset_array);

    (contiguous_array, matches)
}

fn get_distance(i: usize, nrows: usize, core_genes: usize, matches: f64, core: bool,
    contiguous_array: &ndarray::Array2<u8>,
    ncols: usize
) -> Vec<f64> {
    let row1 = contiguous_array.index_axis(Axis(0), i);
    let row1_slice = row1.as_slice().unwrap().to_vec();  // Avoid multiple calls

    (0..nrows)
        .filter_map(|j| {
            if j == i {
                return None;  // Skip self-comparison
            }

            let row2 = contiguous_array.index_axis(Axis(0), j);
            let row2_slice = row2.as_slice().unwrap().to_vec();  // Single call

            let pair_distance = if core {
                let distance = hamming::distance_fast(&row1_slice, &row2_slice).unwrap();
                distance as f64 / (ncols as f64)
            } else {
                let (intersection, union) = jaccard_distance(&row1_slice, &row2_slice);
                1.0 - ((intersection as f64 + matches + core_genes as f64)
                    / (union as f64 + matches + core_genes as f64))
            };

            Some(pair_distance)  // Use `filter_map` instead of `map`
        })
        .collect()
}

struct Population {
    pop: Array2<u8>,
    core : bool,
    core_vec : Vec<Vec<u8>>,
    core_genes : usize,
    avg_gene_freq : f64,
}

// stacks vector of arrays into 2D array
fn to_array2<T: Copy>(source: Vec<Array1<T>>) -> Result<Array2<T>, impl std::error::Error> {
    let width = source.len();
    let flattened: Array1<T> = source.into_iter().flat_map(|row| row.to_vec()).collect();
    let height = flattened.len() / width;
    flattened.into_shape((width, height))
}

impl Population {
    fn new(size: usize, allele_count: usize, max_variants: u8, core : bool, avg_gene_freq: f64, rng : &mut StdRng, core_genes : usize, acc_sampling_vec: &Vec<f64>) -> Self {
        //let mut pop = Array2::<u8>::zeros((size, allele_count));

        // for multithreading
        let _index = AtomicUsize::new(0);
        let _update_rng = AtomicUsize::new(0);

        // generate vector of vectors to hold information in
        let start: Array1<u8> = Array1::zeros(allele_count);
        let pop_vec: Vec<Array1<u8>> = std::iter::repeat(start)
            .take(size)
            .collect();

        // convert vector into 2D array
        let mut pop = to_array2(pop_vec).unwrap();

        if core {
            // ensure core is same across all isolates
            let allele_vec: Vec<u8> = (0..allele_count)
            .map(|_| rng.gen_range(0..max_variants)).collect();

            pop.axis_iter_mut(Axis(0)).into_par_iter().for_each(|mut row| {
                for j in 0..allele_count
                {
                    row[j] = allele_vec[j];
                }
            }
            );
        } else {            
            let mut acc_array: Array1<u8> = Array1::zeros(allele_count);
            for j in 0..allele_count
            {
                //let sample_prop = acc_sampling_vec[j];
                let sampled_value: f64 = rng.gen();
                acc_array[j] = if sampled_value < avg_gene_freq { 1 } else { 0 };
            }

            pop.axis_iter_mut(Axis(0)).into_par_iter().for_each(|mut row| {
                
                // let mut thread_rng = rng.clone();
                // let current_index = _index.fetch_add(1, Ordering::SeqCst);
                // //let thread_index = rayon::current_thread_index();
                // //print!("{:?} ", thread_index);

                // // Jump the state of the generator for this thread
                // for _ in 0..current_index {
                //     thread_rng.gen::<u64>(); // Discard some numbers to mimic jumping
                // }

                // ensure all accessory genomes are identical at start
                for j in 0..allele_count
                {
                    //let sample_prop = acc_sampling_vec[j];
                    //et sampled_value: f64 = thread_rng.gen();
                    //row[j] = if sampled_value < avg_gene_freq { 1 } else { 0 };
                    row[j] = acc_array[j];
                    //_update_rng.fetch_add(1, Ordering::SeqCst);
                }
            }
            );
        }

        // // update rng in place
        // let rng_index: usize = _update_rng.load(Ordering::SeqCst);
        // //print!("{:?} ", rng_index);
        // for _ in 0..rng_index {
        //     rng.gen::<u64>(); // Discard some numbers to mimic jumping
        // }

        let core_vec: Vec<Vec<u8>> = vec![vec![1, 2, 3],
                                          vec![0, 2, 3],
                                          vec![0, 1, 3],
                                          vec![0, 1, 2]];

        Self {
            pop,
            core,
            core_vec,
            core_genes,
            avg_gene_freq,
        }
    }

    fn calc_gene_freq (&mut self) -> f64 {
        // Calculate the proportion of 1s for each row
        let proportions: Vec<f64> = self.pop.axis_iter(Axis(0))
        .map(|row| {
            let sum: usize = row.iter().map(|&x| x as usize).sum();
            let count = row.len();
            sum as f64 / count as f64
        })
        .collect();

        // Sum all the elements in the vector
        let sum: f64 = proportions.iter().sum();

        // Calculate the number of elements in the vector
        let count = proportions.len();

        // Calculate the average
        let average: f64 = sum as f64 / count as f64;

        average
    }

    fn sample_indices (&mut self, rng : &mut StdRng, avg_gene_num: i32, avg_pairwise_dists : Vec<f64>) -> Vec<usize> {
        // Calculate the proportion of 1s for each row
        let num_genes: Vec<i32> = self.pop.axis_iter(Axis(0))
        .map(|row| {
            let sum: i32 = row.iter().map(|&x| x as i32).sum();
            sum
            //let count = row.len();
            //sum as f64 / count as f64
        })
        .collect();

        //println!("proportions:\n{:?}", proportions);

        // Calculate the differences from avg_gene_freq
        let differences: Vec<i32> = num_genes.iter()
        .map(|&n_genes| (n_genes - avg_gene_num).abs())
        .collect();

        //println!("differences:\n{:?}", differences);

        // Convert differences to weights (lower difference should have higher weight)
        //let max_diff = differences.iter().cloned().fold(0./0., f64::max);
        let mut weights: Vec<f64> = differences.iter()
            .map(|&diff| 0.99_f64.powi(diff) ) // based on https://pmc.ncbi.nlm.nih.gov/articles/instance/5320679/bin/mgen-01-38-s001.pdf
            .collect();

        // update weights with average pairwise distance
        for i in 0..weights.len()
        {
            weights[i] *= avg_pairwise_dists[i];
        }

        //println!("max_diff:\n{:?}", max_diff);
        //println!("weights:\n{:?}", weights);

        // Create a WeightedIndex distribution based on weights
        let dist = WeightedIndex::new(&weights).unwrap();

        // Sample rows based on the distribution
        let sampled_indices: Vec<usize> = (0..self.pop.nrows())
        .map(|_| dist.sample(rng))
        .collect();
        
        //println!("sampled_indices:\n{:?}", sampled_indices);

        sampled_indices
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

    fn mutate_alleles(&mut self, mutations_vec : &Vec<i32>, rng : &mut StdRng, weighted_dist: &Vec<WeightedIndex<f32>>) {
        // index for random number generation
        let _index = AtomicUsize::new(0);
        let _update_rng = AtomicUsize::new(0);

        for site_idx in 0..mutations_vec.len() {

            let mutations = mutations_vec[site_idx];
            let poisson = Poisson::new(mutations as f64).unwrap();

            // avoid case where no mutations are allowed
            if mutations == 0 {
                continue;
            }

            if self.core == false {
                // generate Poisson sampler
                self.pop.axis_iter_mut(Axis(0)).into_par_iter().for_each(|mut row| {
                    // thread-specific random number generator
                    let mut thread_rng = rng.clone();
                    let current_index = _index.fetch_add(1, Ordering::SeqCst);
                    //let thread_index = rayon::current_thread_index();
                    //print!("{:?} ", thread_index);

                    // Jump the state of the generator for this thread
                    for _ in 0..current_index {
                        thread_rng.gen::<u64>(); // Discard some numbers to mimic jumping
                    }
                    
                    // sample from Poisson distribution for number of sites to mutate in this isolate
                    let n_sites = poisson.sample(&mut thread_rng) as usize;
                    _update_rng.fetch_add(1, Ordering::SeqCst);

                    // iterate for number of mutations required to reach mutation rate
                    for _ in 0..n_sites {
                        // sample new site to mutate
                        let mutant_site = weighted_dist[site_idx].sample(&mut thread_rng);
                        _update_rng.fetch_add(1, Ordering::SeqCst);
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
                        let current_index = _index.fetch_add(1, Ordering::SeqCst);
                        //let thread_index = rayon::current_thread_index();
                        //print!("{:?} ", thread_index);

                        // Jump the state of the generator for this thread
                        for _ in 0..current_index {
                            thread_rng.gen::<u64>(); // Discard some numbers to mimic jumping
                        }

                        // sample from Poisson distribution for number of sites to mutate in this isolate
                        let n_sites = thread_rng.sample(poisson) as usize;
                        _update_rng.fetch_add(1, Ordering::SeqCst);

                        // iterate for number of mutations required to reach mutation rate
                        for _ in 0..n_sites {
                            // sample new site to mutate
                            let mutant_site = weighted_dist[site_idx].sample(&mut thread_rng);
                            _update_rng.fetch_add(1, Ordering::SeqCst);

                            // get possible values to mutate to, must be different from current value
                            let value = row[mutant_site];
                            let values = &self.core_vec[value as usize];

                            // sample new allele
                            let new_allele = values.iter().choose_multiple(&mut thread_rng, 1)[0];
                            _update_rng.fetch_add(1, Ordering::SeqCst);

                            // set value in place
                            row[mutant_site] = *new_allele;
                        }
                    });
            }
        }
        // update rng in place
        let rng_index: usize = _update_rng.load(Ordering::SeqCst);
        //print!("{:?} ", rng_index);
        for _ in 0..rng_index {
            rng.gen::<u64>(); // Discard some numbers to mimic jumping
        }

    }

    fn recombine(&mut self, recombinations_vec : &Vec<f64>, rng : &mut StdRng, locus_weights: &Vec<Vec<f32>>) {
        // index for random number generation
        let _index = AtomicUsize::new(0);
        let _update_rng = AtomicUsize::new(0);
        
        for site_idx in 0..recombinations_vec.len() {
            
            let n_recombinations =  recombinations_vec[site_idx];

            let poisson_recomb = Poisson::new(n_recombinations as f64).unwrap();

            // Preallocate results vector with one entry per row
            //let mut loci: Vec<Vec<usize>> = vec![Vec::new(); self.pop.nrows()];
            let loci = Arc::new(RwLock::new(vec![Vec::new(); self.pop.nrows()]));
            let values = Arc::new(RwLock::new(vec![Vec::new(); self.pop.nrows()]));
            let recipients = Arc::new(RwLock::new(vec![Vec::new(); self.pop.nrows()]));
    
            // let contiguous_array:ndarray::ArrayBase<ndarray::OwnedRepr<u8>, ndarray::Dim<[usize; 2]>>;
            // let matches:f64; 
    
            // // get mutation matrix
            // match pangenome_matrix {
            //     Some(matrix) => {
            //         (contiguous_array, matches) =  get_variable_loci(false, &matrix);
            //     }
            //     None => {
            //         (contiguous_array, matches) =  get_variable_loci(false, &self.pop);
            //     }
            // }
    
            // recipient distribution, minus one to avoid comparison with self
            let dist: Uniform<usize> = Uniform::new(0, self.pop.nrows() - 1);
    
            // for each genome, determine which positions are being transferred
            self.pop.axis_iter(Axis(0)).into_par_iter().enumerate().for_each(|(row_idx , row)| {
                //use std::time::Instant;
                //let now = Instant::now();
                
                // thread-specific random number generator
                let mut thread_rng = rng.clone();
                let current_index = _index.fetch_add(1, Ordering::SeqCst);
                // Jump the state of the generator for this thread
                for _ in 0..current_index {
                    thread_rng.gen::<u64>(); // Discard some numbers to mimic jumping
                }
                
                // sample from Poisson distribution for number of sites to mutate in this isolate
                let n_sites = poisson_recomb.sample(&mut thread_rng) as usize;
                //let n_targets = poisson_recip.sample(&mut thread_rng) as usize;
                _update_rng.fetch_add(1, Ordering::SeqCst);
    
                // get sampling weights for each pairwise comparison
                // TODO remove this, jsut have same distance for all individuals
                //let binding = get_distance(row_idx, self.pop.nrows(), self.core_genes, matches, false, &contiguous_array, self.pop.ncols());
                //let mut elapsed = now.elapsed();
        
                //println!("finished distances: {}, {:.2?}", row_idx, elapsed);
                //let binding = vec![0.1; self.pop.nrows()];
                //let i_distances = binding
                //.iter().map(|i| {1.0 - i});
                //let sample_dist = WeightedIndex::new(i_distances).unwrap();
    
                //elapsed = now.elapsed();
                //println!("finished sampling dist: {}, {:.2?}", row_idx, elapsed);
    
                // Sample rows based on the distribution, adjusting as self comparison not conducted
                let sampled_recipients: Vec<usize> = (0..n_sites).map(|_| dist.sample(&mut thread_rng)).map(|value| value + (value >= row_idx) as usize).collect();
                let mut sampled_loci: Vec<usize> = Vec::with_capacity(n_sites);
    
                // for _ in 0..n_sites {
                //     let value = sample_dist.sample(&mut thread_rng);
                //     sampled_recipients.push(value + (value >= row_idx) as usize);
                // }
    
                //elapsed = now.elapsed();
                //println!("finished sampling total: {}, {:.2?}", row_idx, elapsed);
    
                //let sampled_recipients: Vec<usize> = vec![1, 7, 20, 705, 256];
    
                _update_rng.fetch_add(n_sites, Ordering::SeqCst);
    
                let mut sampled_values: Vec<u8> = vec![1; n_sites];
                // get non-zero indices
                if self.core == false {
                    //if accessory, set elements with no genes to 0
                    let mut non_zero_weights: Vec<f32> = locus_weights[site_idx].clone();
                    let mut total_0: usize = 0;
                    for (idx, &val) in row.indexed_iter() {
                        if val == 0  {
                            non_zero_weights[idx] = 0.0;
                            total_0 += 1
                        }
                    }
    
                    // // get all sites to be recombined
                    // sampled_loci = (0..n_sites)
                    // .map(|_| *non_zero_indices.choose(&mut thread_rng).unwrap()) // Sample with replacement
                    // .collect();
    
                    // ensure some non-zero values present
                    if total_0 < non_zero_weights.len() {
                        let locus_weighted_dist: WeightedIndex<f32> = WeightedIndex::new(non_zero_weights).unwrap();
    
                        // iterate for number of mutations required to reach mutation rate, include deletions and insertions
                        sampled_loci = thread_rng
                            .sample_iter(locus_weighted_dist)
                            .take(n_sites)
                            .collect();
                    }
                    
    
                } else {
                    
                    // sampled_loci = (0..n_sites)
                    // .map(|_| row.indexed_iter().map(|(idx, _)| idx).choose(&mut thread_rng).unwrap()) // Sample with replacement
                    // .collect();
    
                    sampled_loci = thread_rng
                        .sample_iter(rand::distributions::Uniform::new(0, self.pop.ncols()))
                        .take(n_sites)
                        .collect();
    
                    // assign site value from row
                    for site in 0..n_sites {
                        sampled_values[site] = row[sampled_loci[site]];
                    }
                }
    
                //elapsed = now.elapsed();
                //println!("finished getting sites total: {}, {:.2?}", row_idx, elapsed);
    
                // update the rng
                _update_rng.fetch_add(n_sites, Ordering::SeqCst);
    
                // assign values
                {
                    let mut mutex = loci.write().unwrap();  // Lock for writing
                    let entry = &mut mutex[row_idx];  // Now you can index safely
                    *entry = sampled_loci;
                }
                {
                    let mut mutex = values.write().unwrap();  // Lock for writing
                    let entry = &mut mutex[row_idx];  // Now you can index safely
                    *entry = sampled_values;
                }
                {
                    let mut mutex = recipients.write().unwrap();  // Lock for writing
                    let entry = &mut mutex[row_idx];  // Now you can index safely
                    *entry = sampled_recipients;
                }
    
                //elapsed = now.elapsed();
                //println!("finished entering data: {}, {:.2?}", row_idx, elapsed);
            }  
            );

            // go through entries in loci, values and recipients, mutating the rows in each case
            // randomise order in which rows are moved through 
            let mut row_indices: Vec<usize> = (0..self.pop.nrows()).collect();
            row_indices.shuffle(rng);

            for pop_idx in row_indices
            {
                // sample for given donor
                let sampled_loci: Vec<usize> = loci.write().unwrap()[pop_idx].to_vec();
                let sampled_recipients: Vec<usize> = recipients.write().unwrap()[pop_idx].to_vec();
                let sampled_values: Vec<u8> = values.write().unwrap()[pop_idx].to_vec();

                // println!("index: {}", pop_idx);
                // println!("sampled_loci: {:?}", sampled_loci);
                // println!("sampled_recipients: {:?}", sampled_recipients);
                // println!("sampled_values: {:?}", sampled_values);

                // update recipients in place if any recombinations allowed
                if sampled_loci.len() > 0 {
                    Zip::from(&sampled_recipients)
                    .and(&sampled_loci)
                    .and(&sampled_values)
                    .for_each(|&row_idx, &col_idx, &value| {
                        self.pop[[row_idx, col_idx]] = value;
                    });
                }
            }
            
            // update rng in place
            let rng_index: usize = _update_rng.load(Ordering::SeqCst);
            //print!("{:?} ", rng_index);
            for _ in 0..rng_index {
                rng.gen::<u64>(); // Discard some numbers to mimic jumping
            }
        }

    }

    fn average_distance(&mut self) -> Vec<f64> {
        let (contiguous_array, matches) =  get_variable_loci(self.core, &self.pop);
        
        let range = 0..self.pop.nrows();
        let distances: Vec<f64> = range.into_par_iter().map(|i| {

            let i_distances = get_distance(i, self.pop.nrows(), self.core_genes, matches, self.core, &contiguous_array, self.pop.ncols());
            
            let mut _final_distance = i_distances.iter().sum::<f64>() / i_distances.len() as f64;
            
            // ensure no zero distances that may cause no selection of isolates.
            if _final_distance == 0.0
            {
                _final_distance = MIN_POSITIVE;
            }

            _final_distance
        }).collect();

        //println!("new distances:\n{:?}", distances);
        distances
        }

    fn pairwise_distances(&mut self, max_distances : usize, range1: &Vec<usize>, range2: &Vec<usize>) -> Vec<f64> {
        let (contiguous_array, matches) =  get_variable_loci(self.core, &self.pop);

        //let mut idx = 0;
        let range = 0..max_distances;
        let distances: Vec<_> = range.into_par_iter().map(|current_index| {
            let i = range1[current_index];
            let j = range2[current_index];

            //println!("i:\n{:?}", i);
            //println!("j:\n{:?}", j);
            
            let row1 = contiguous_array.index_axis(Axis(0), i);
            let row2 = contiguous_array.index_axis(Axis(0), j);
            let row1_slice = row1.as_slice().unwrap().to_vec();
            let row2_slice = row2.as_slice().unwrap().to_vec(); 

            //println!("rowi:\n{:?}", row1);
            //println!("rowj:\n{:?}", row2);

            let mut _final_distance: f64 = 0.0;

            if self.core == true {
                let distance = hamming::distance_fast(&row1_slice, &row2_slice).unwrap();
                _final_distance = distance as f64 / (self.pop.ncols() as f64);
            } else {
                let (intersection, union) = jaccard_distance(&row1_slice, &row2_slice);
                _final_distance = 1.0 - ((intersection as f64 + matches + self.core_genes as f64) / (union as f64 + matches + self.core_genes as f64));
            }
            //println!("_final_distance:\n{:?}", _final_distance);
            _final_distance
        }).collect();
        distances
        }
}

fn main() -> io::Result<()> {

    // Define the command-line arguments using clap
    let matches = Command::new("pansim")
    .version("0.0.2")
    .author("Samuel Horsfield shorsfield@ebi.ac.uk")
    .about("Runs Wright-Fisher simulation, simulating neutral core genome evolution and two-speed accessory genome evolution.")
    .arg(Arg::new("pop_size")
        .long("pop_size")
        .help("Number of individuals in population.")
        .required(false)
        .default_value("1000"))
    .arg(Arg::new("core_size")
        .long("core_size")
        .help("Number of nucleotides in core genome.")
        .required(false)
        .default_value("1200000"))
    .arg(Arg::new("pan_genes")
        .long("pan_genes")
        .help("Total number of genes in pangenome (core + accessory).")
        .required(false)
        .default_value("6000"))
    .arg(Arg::new("core_genes")
        .long("core_genes")
        .help("Number of core genes in pangenome.")
        .required(false)
        .default_value("2000"))
    .arg(Arg::new("avg_gene_freq")
        .long("avg_gene_freq")
        .help("Average proportion of genes in pangenome present in an individual. Includes core and accessory genes")
        .required(false)
        .default_value("0.5"))
    .arg(Arg::new("n_gen")
        .long("n_gen")
        .help("Number of generations to simulate.")
        .required(false)
        .default_value("100"))
    .arg(Arg::new("max_distances")
        .long("max_distances")
        .help("Maximum number of pairwise distances to calculate.")
        .required(false)
        .default_value("100000"))
    .arg(Arg::new("core_mu")
        .long("core_mu")
        .help("Maximum average pairwise core distance to achieve by end of simulation.")
        .required(false)
        .default_value("0.05"))
    .arg(Arg::new("HR_rate")
        .long("HR_rate")
        .help("Homologous recombination rate, as number of core sites transferred per core genome mutation.")
        .required(false)
        .default_value("0.05"))
    .arg(Arg::new("HGT_rate")
        .long("HGT_rate")
        .help("HGT rate, as number of accessory sites transferred per core genome mutation.")
        .required(false)
        .default_value("0.05"))
    .arg(Arg::new("competition")
        .long("competition")
        .help("Adds competition based on average pairwise genome distance.")
        .required(false)
        .takes_value(false))
    .arg(Arg::new("rate_genes1")
        .long("rate_genes1")
        .help("Proportion of accessory pangenome that mutates per generation in gene compartment 1. Must be >= 0.0")
        .required(false)
        .default_value("1.0"))
    .arg(Arg::new("rate_genes2")
        .long("rate_genes2")
        .help("Proportion of accessory pangenome that mutates per generation in gene compartment 2. Must be >= 0.0")
        .required(false)
        .default_value("prop_genes2"))
    .arg(Arg::new("prop_genes2")
        .long("prop_comp2")
        .help("Proportion of pangenome made up of compartment 2 genes. Must be 0.0 <= X <= 0.5")
        .required(false)
        .default_value("2.0"))
    .arg(Arg::new("seed")
        .long("seed")
        .help("Seed for random number generation.")
        .required(false)
        .default_value("0"))
    .arg(Arg::new("outpref")
        .long("outpref")
        .help("Output prefix path.")
        .required(false)
        .default_value("distances"))
    .arg(Arg::new("print_dist")
        .long("print_dist")
        .required(false)
        .takes_value(false))
        .help("Print per-generation average pairwise distances.")
    .arg(Arg::new("threads")
        .long("threads")
        .help("Number of threads.")
        .required(false)
        .default_value("1"))
    .arg(Arg::new("verbose")
        .long("verbose")
        .help("Prints generation and time to completion")
        .required(false)
        .takes_value(false))
    .get_matches();

    // Set the argument to a variable
    let pop_size: usize = matches.value_of_t("pop_size").unwrap();
    let core_size: usize = matches.value_of_t("core_size").unwrap();
    let pan_genes: usize = matches.value_of_t("pan_genes").unwrap();
    let core_genes: usize = matches.value_of_t("core_genes").unwrap();
    let mut avg_gene_freq: f64 = matches.value_of_t("avg_gene_freq").unwrap();
    let HR_rate: f64 = matches.value_of_t("HR_rate").unwrap();
    let HGT_rate: f64 = matches.value_of_t("HGT_rate").unwrap();
    let n_gen: i32 = matches.value_of_t("n_gen").unwrap();
    let outpref = matches.value_of("outpref").unwrap_or("distances");
    let max_distances: usize = matches.value_of_t("max_distances").unwrap();
    let core_mu: f64 = matches.value_of_t("core_mu").unwrap();
    let rate_genes1: f64 = matches.value_of_t("rate_genes1").unwrap();
    let rate_genes2: f64 = matches.value_of_t("rate_genes2").unwrap();
    let prop_comp2: f32 = matches.value_of_t("prop_comp2").unwrap();
    let mut n_threads: usize = matches.value_of_t("threads").unwrap();
    let verbose = matches.is_present("verbose");
    let competition = matches.is_present("competition");
    let seed: u64 = matches.value_of_t("seed").unwrap();
    let print_dist: bool = matches.is_present("print_dist");

    //let verbose = true;

    // time testing
    //use std::time::Instant;
    //let now = Instant::now();

    // validate all variables
    if core_genes > pan_genes {
        println!("core_genes must be less than or equal to pan_size");
        return Ok(())
    }

    if (HR_rate < 0.0 || HGT_rate < 0.0) {
        println!("HR_rate and HGT_rate must be above 0.0");
        println!("HR_rate: {}", HR_rate);
        println!("HGT_rate: {}", HGT_rate);
        return Ok(())
    }

    if rate_genes1 < 0.0 || rate_genes2 < 0.0 {
        println!("rate_genes1 and rate_genes2 must be >= 0.0");
        println!("rate_genes1: {}", rate_genes1);
        println!("rate_genes2: {}", rate_genes2);
        return Ok(())
    }

    if prop_comp2 < 0.0 || prop_comp2 > 0.5 {
        println!("prop_comp2 must be 0.0 <= prop_comp2 <= 0.5");
        println!("prop_comp2: {}", prop_comp2);
        return Ok(())
    }

    if (pop_size < 1) || (core_size < 1) || (pan_genes < 1) || (n_gen < 1) || (max_distances < 1) {
        println!("pop_size, core_size, pan_genes, n_gen and max_distances must all be above 1");
        println!("pop_size: {}", pop_size);
        println!("core_size: {}", core_size);
        println!("pan_genes: {}", pan_genes);
        println!("n_gen: {}", n_gen);
        println!("max_distances: {}", max_distances);
        return Ok(())
    }

    if (core_mu < 0.0) || (core_mu > 1.0) {
        println!("core_mu must be between 0.0 and 1.0");
        println!("core_mu: {}", core_mu);
        return Ok(())
    }

    if (avg_gene_freq <= 0.0) || (avg_gene_freq > 1.0) {
        println!("avg_gene_freq must be above 0.0 and below or equal to 1.0");
        println!("avg_gene_freq: {}", avg_gene_freq);
        return Ok(())
    }


    if n_threads < 1 {
        n_threads = 1;
    }

    // enable multithreading
    rayon::ThreadPoolBuilder::new().num_threads(n_threads).build_global().unwrap();

    let pan_size = pan_genes - core_genes;

    // adjust avg_gene_freq by size of core genome
    // determine proportion of accessory genome that should be present on average
    let core_prop : f64 = core_genes as f64 / pan_genes as f64;
    let acc_prop : f64 = 1.0 - core_prop;
    avg_gene_freq = (avg_gene_freq - core_prop) / acc_prop;
    if avg_gene_freq < 0.0 {
        avg_gene_freq = 0.0;
    }
    if verbose {
        println!("avg_gene_freq adjusted to {}", avg_gene_freq);
    }
    let avg_gene_num: i32 = (avg_gene_freq * pan_size as f64).round() as i32;
    
    // calculate number of mutations per genome per generation, should this be whole pangenome or just accessory genes?
    let n_core_mutations = vec![(((core_size as f64 * core_mu) / n_gen as f64) / 2.0).ceil() as i32];
    

    // calculate average recombinations per genome
    let n_recombinations_core: Vec<f64> = vec![((n_core_mutations[0] as f64 * HR_rate)).round()];
    let n_recombinations_pan_total = ((n_core_mutations[0] as f64 * HGT_rate)).round();
    let n_recombinations_pan_gene1 = n_recombinations_pan_total * (rate_genes1 / (rate_genes1 + rate_genes2));
    let n_recombinations_pan_gene2 = n_recombinations_pan_total * (rate_genes2 / (rate_genes1 + rate_genes2));
    let n_recombinations_pan: Vec<f64> = vec![n_recombinations_pan_gene1, n_recombinations_pan_gene2];

    // set weights for sampling of sites
    let core_weights : Vec<Vec<f32>> = vec![vec![1.0; core_size]; 1];
    let mut pan_weights : Vec<Vec<f32>> = vec![vec![0.0; pan_size]; 2];

    // calculate sites for fast accessory genome
    let num_gene1_sites = (pan_size as f32 * (1.0 - prop_comp2)).round() as usize;
    let n_pan_mutations_gene1 = ((pan_size as f64 * rate_genes1)).ceil() as i32;
    let mut n_pan_mutations_gene2 = ((pan_size as f64 * rate_genes2)).ceil() as i32;

    // set number of gene 1 mutations to 0 if there are no gene 2 present
    if num_gene1_sites == pan_size {
        n_pan_mutations_gene2 = 0;
    }

    let n_pan_mutations = vec![n_pan_mutations_gene1, n_pan_mutations_gene2];
    
    // create weights for either rate compartments
    // for gene rate 1
    for i in 0..num_gene1_sites {
        pan_weights[0][i] = 1.0;
    }
    // gene rate 2
    for i in num_gene1_sites..pan_size {
        pan_weights[1][i] = 1.0;
    }

    let mut rng: StdRng = StdRng::seed_from_u64(seed);

    // generate sampling distribution for genes in accessory genome
    let acc_sampling_vec = sample_beta(pan_size, &mut rng);

    let mut core_genome = Population::new(pop_size, core_size, 4, true, avg_gene_freq, &mut rng, core_genes, & acc_sampling_vec); // core genome alignment
    let mut pan_genome = Population::new(pop_size, pan_size, 2, false, avg_gene_freq, &mut rng, core_genes, & acc_sampling_vec); // pangenome alignment

    // weighted distribution samplers
    let core_weighted_dist: Vec<WeightedIndex<f32>> = core_weights.clone().into_iter().map(|dist| WeightedIndex::new(dist.clone()).unwrap()).collect();
    let pan_weighted_dist: Vec<WeightedIndex<f32>> = pan_weights.clone().into_iter().map(|dist| WeightedIndex::new(dist.clone()).unwrap()).collect();

    // hold pairwise core and accessory distances per generation
    let mut avg_acc_dist = vec![0.0; n_gen as usize];
    let mut avg_core_dist = vec![0.0; n_gen as usize];
    let mut std_acc_dist = vec![0.0; n_gen as usize];
    let mut std_core_dist = vec![0.0; n_gen as usize];


    // generate random numbers to sample indices
    // TODO make it so that equivalent distances aren't sample, sample with replacement from one?
    let range1: Vec<usize> = (0..max_distances).map(|_| rng.gen_range(0..pop_size)).collect();
    let mut range2: Vec<usize> = vec![0; max_distances];

    let mut i2 = 0;
    // sample same range, ensure self-comparisons not included
    for i1 in range1.clone() {
        let mut entry = rng.gen_range(0..pop_size - 1);
        if entry >= i1 {entry += 1};
        range2[i2] = entry;
        i2 += 1;
    }
    
    for j in 0..n_gen { // Run for n_gen generations
        //let now_gen = Instant::now();
        
        // sample new individuals if not at first generation
        if j > 0 {
            //let sampled_individuals: Vec<usize> = (0..pop_size).map(|_| rng.gen_range(0..pop_size)).collect();
            let mut avg_pairwise_dists = vec![1.0; pop_size];
            
            // include competition
            if competition == true
            {
                avg_pairwise_dists = pan_genome.average_distance();
            }
            
            let sampled_individuals = pan_genome.sample_indices(&mut rng, avg_gene_num, avg_pairwise_dists);
            core_genome.next_generation(& sampled_individuals);
            //println!("finished copying core genome {}", j);
            pan_genome.next_generation(& sampled_individuals);
            //println!("finished copying pangenome {}", j);
        }
        
        // if at final generation, just sample, otherwise mutate
        if j < (n_gen - 1) {
            // mutate core genome
            //println!("started {}", j);
            core_genome.mutate_alleles(&n_core_mutations, &mut rng, &core_weighted_dist);

            //println!("finished mutating core genome {}", j);
            pan_genome.mutate_alleles(&n_pan_mutations, &mut rng, &pan_weighted_dist);
            //println!("finished mutating pangenome {}", j);

            // recombine populations
            if HR_rate > 0.0 {
                core_genome.recombine(&n_recombinations_core, &mut rng, &core_weights);
            }
            if HGT_rate > 0.0 {
                pan_genome.recombine(&n_recombinations_pan, &mut rng, &pan_weights);
            }

        } else {
            let final_avg_gene_freq = pan_genome.calc_gene_freq();
            if verbose {
                println!("final avg_gene_freq: {}", final_avg_gene_freq);
            }
            
            // else calculate hamming and jaccard distances
            let core_distances = core_genome.pairwise_distances(max_distances, &range1, &range2);
            let acc_distances = pan_genome.pairwise_distances(max_distances, &range1, &range2);

            let mut output_file = outpref.to_owned();
            let extension: &str = ".tsv";
            output_file.push_str(extension);
            let mut file = File::create(output_file)?;

            // Iterate through the vectors and write each pair to the file
            for (core, acc) in core_distances.iter().zip(acc_distances.iter()) {
                writeln!(file, "{}\t{}", core, acc);
            }
        }

        // get average distances
        if print_dist {
            let core_distances = core_genome.pairwise_distances(max_distances, &range1, &range2);
            let acc_distances = pan_genome.pairwise_distances(max_distances, &range1, &range2);

            let mut std_core = 0.0;
            let mut avg_core = 0.0;
            (std_core, avg_core) = standard_deviation(&core_distances);
            
            let mut std_acc = 0.0;
            let mut avg_acc = 0.0;
            (std_acc, avg_acc) = standard_deviation(&acc_distances);
            
            avg_core_dist[j as usize] = avg_core;
            avg_acc_dist[j as usize] = avg_acc;

            std_core_dist[j as usize] = std_core;
            std_acc_dist[j as usize] = std_acc;
        }

        //let elapsed = now_gen.elapsed();
        if verbose {
            println!("Finished gen: {}", j + 1);    
        }
        //println!("Elapsed: {:.2?}", elapsed);
    }

    // print per generation distances
    if print_dist {
        let mut output_file = outpref.to_owned();
        let extension: &str = "_per_gen.tsv";
        output_file.push_str(extension);

        let mut file = File::create(output_file)?;

        // Iterate through the vectors and write each pair to the file
        for (avg_core, avg_acc, std_core, std_acc) in avg_core_dist.iter().zip(avg_acc_dist.iter()).zip(std_core_dist.iter()).zip(std_acc_dist.iter()).map(|(((w, x), y), z)| (w, x, y, z)) {
            writeln!(file, "{}\t{}\t{}\t{}", avg_core, std_core, avg_acc, std_acc);
        }
    }

    //let elapsed = now.elapsed();
    
    // if verbose {
    //     println!("Total elapsed: {:.2?}", elapsed);
    // }

    // println!("Total elapsed: {:.2?}", elapsed);

    return Ok(())
}
