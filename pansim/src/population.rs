use rand_distr::Beta;

//use rand_distr::Distribution;
use rayon::prelude::*;

use statrs::distribution::Poisson;

use rand::distributions::Distribution;
use rand::distributions::Uniform;
use rand::distributions::WeightedIndex;
use rand::rngs::StdRng;
use rand::seq::IteratorRandom;
use rand::seq::SliceRandom;
use rand::{Rng};

use ndarray::Zip;
use ndarray::{s, Array1, Array2, Axis};
use std::f64::MIN_POSITIVE;
use std::f64::MAX;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, RwLock};

use logsumexp::LogSumExp;
use crate::distances::*;
use safe_arch::*;

use std::fs::File;
use std::io::{self, Write};

use std::usize;

// fn jaccard_distance_test(v1: &[u8], v2: &[u8]) -> (u32, u32) {
//     assert_eq!(v1.len(), v2.len(), "Vectors must have the same length.");

//     let mut intersection: u32 = 0;
//     let mut union: u32 = 0;


//     for (&x, &y) in v1.iter().zip(v2.iter()) {
//         if x == 1 || y == 1 {
//             union += 1;
//             if x == 1 && y == 1 {
//                 intersection += 1;
//             }
//         }
//     }
//     (intersection, union)
// }

fn safe_pow(base: f64, exp: f64) -> f64 {
    let result = base.powf(exp);
    if result.is_infinite() {
        if result.is_sign_positive() {
            MAX
        } else {
            MIN_POSITIVE
        }
    } else {
        result
    }
}

pub fn normalize_avg_dists(numbers: &[f64]) -> Vec<f64> {
    let sum = numbers.iter().sum::<f64>() as f64;

    //let minvalue = numbers.iter().fold(f64::INFINITY, |a, &b| a.min(b));

    let normalised: Vec<f64> = numbers.iter().map(|&x| (x / sum)).collect();
    normalised
}

fn average(numbers: &[f64]) -> f64 {
    numbers.iter().sum::<f64>() as f64 / numbers.len() as f64
}

pub fn standard_deviation(values: &[f64]) -> (f64, f64) {
    let mean = average(values);

    let sum_of_squares: f64 = values.iter().map(|&x| (x - mean).powi(2)).sum();

    let variance = sum_of_squares / (values.len() as f64);
    (variance.sqrt(), mean)
}

pub fn sample_beta(num_samples: usize, rng: &mut StdRng) -> Vec<f64> {
    let mut return_vec: Vec<f64> = vec![0.0; num_samples];

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

fn get_distance(
    i: usize,
    nrows: usize,
    core_genes: usize,
    matches: f64,
    core: bool,
    contiguous_array: &ndarray::Array2<u8>,
    ncols: usize,
) -> Vec<f64> {
    let row1 = contiguous_array.index_axis(Axis(0), i);
    let row1_slice = row1.as_slice().unwrap(); // Avoid multiple calls

    (0..nrows)
        .filter_map(|j| {
            if j == i {
                return None; // Skip self-comparison
            }

            let row2 = contiguous_array.index_axis(Axis(0), j);
            let row2_slice = row2.as_slice().unwrap(); // Single call

            let pair_distance = if core {
                let distance = hamming_bitwise_fast(row1_slice, row2_slice) / 2;
                distance as f64 / (ncols as f64)
            } else {
                
                let (intersection, union) = jaccard_distance_fast(row1_slice, row2_slice);
                // let (intersection_test, union_test) = jaccard_distance_test(row1_slice, row2_slice);
                // println!("intersection_test: {:?} intersection: {:?}", intersection_test, intersection);
                // println!("union_test: {:?} union: {:?}", union_test, union);
                1.0 - ((intersection as f64 + matches + core_genes as f64)
                    / (union as f64 + matches + core_genes as f64))
            };

            Some(pair_distance) // Use `filter_map` instead of `map`
        })
        .collect()
}

// Mapping function
pub fn int_to_base(n: u8) -> char {
    match n {
        1 => 'A',
        2 => 'C',
        4 => 'G',
        8 => 'T',
        _ => 'N', // Unknown
    }
}

pub struct Population {
    pop: Array2<u8>,
    core: bool,
    core_vec: Vec<Vec<u8>>,
    core_genes: usize,
    avg_gene_freq: f64,
}

// stacks vector of arrays into 2D array
fn to_array2<T: Copy>(source: Vec<Array1<T>>) -> Result<Array2<T>, impl std::error::Error> {
    let width = source.len();
    let flattened: Array1<T> = source.into_iter().flat_map(|row| row.to_vec()).collect();
    let height = flattened.len() / width;
    flattened.into_shape((width, height))
}

impl Population {
    pub fn new(
        size: usize,
        allele_count: usize,
        max_variants: u8,
        core: bool,
        avg_gene_freq: f64,
        rng: &mut StdRng,
        core_genes: usize,
        acc_sampling_vec: &Vec<f64>,
    ) -> Self {

        // generate vector of vectors to hold information in
        let start: Array1<u8> = Array1::zeros(allele_count);
        let pop_vec: Vec<Array1<u8>> = std::iter::repeat(start).take(size).collect();

        // convert vector into 2D array
        let mut pop = to_array2(pop_vec).unwrap();

        if core {
            // ensure core is same across all isolates
            let allele_vec: Vec<u8> = (0..allele_count)
                .map(|_| rng.gen_range(0..max_variants))
                .map(|i| {1 << i})
                .collect();

            pop.axis_iter_mut(Axis(0))
                .into_par_iter()
                .for_each(|mut row| {
                    for j in 0..allele_count {
                        row[j] = allele_vec[j];
                    }
                });
        } else {
            let mut acc_array: Array1<u8> = Array1::zeros(allele_count);
            for j in 0..allele_count {
                //let sample_prop = acc_sampling_vec[j];
                let sampled_value: f64 = rng.gen();
                acc_array[j] = if sampled_value < avg_gene_freq { 1 } else { 0 };
            }

            pop.axis_iter_mut(Axis(0))
                .into_par_iter()
                .for_each(|mut row| {

                    // ensure all accessory genomes are identical at start
                    for j in 0..allele_count {
                        row[j] = acc_array[j];
                    }
                });
        }

        let core_vec: Vec<Vec<u8>> =
            vec![vec![2, 4, 8], vec![1, 4, 8], vec![1, 2, 8], vec![1, 2, 4]];

        Self {
            pop,
            core,
            core_vec,
            core_genes,
            avg_gene_freq,
        }
    }

    pub fn calc_gene_freq(&mut self) -> f64 {
        // Calculate the proportion of 1s for each row
        let proportions: Vec<f64> = self
            .pop
            .axis_iter(Axis(0))
            .map(|row| {
                let sum: usize = row.iter().map(|&x| x as usize).sum();
                let count = row.len();
                sum as f64 / count as f64
            })
            .collect();

        //println!("proportions: {:?}", proportions);
        
        // Sum all the elements in the vector
        let sum: f64 = proportions.iter().sum();

        // Calculate the number of elements in the vector
        let count = proportions.len();

        // Calculate the average
        let average: f64 = sum as f64 / count as f64;

        average
    }

    pub fn sample_indices(
        &mut self,
        rng: &mut StdRng,
        avg_gene_num: i32,
        avg_pairwise_dists: Vec<f64>,
        selection_coefficients: &Vec<f64>,
        verbose: bool,
        no_control_genome_size: bool,
        genome_size_penalty: f64,
        competition_strength: f64
    ) -> Vec<usize> {
        // Calculate the proportion of 1s for each row
        let num_genes: Vec<i32> = self
            .pop
            .axis_iter(Axis(0))
            .map(|row| {
                let sum: i32 = row.iter().map(|&x| x as i32).sum();
                sum
                //let count = row.len();
                //sum as f64 / count as f64
            })
            .collect();

        let mut selection_weights: Vec<f64> = vec![1.0; self.pop.nrows()];

        // ensure accessory genome present
        if self.pop.ncols() > 0 {
            
            // TODO generate log sum value for entire row, then do logsumexp across whole selection weights array
            selection_weights = self
                .pop
                .axis_iter(Axis(0))
                .map(|row| {
                    let log_values: Vec<f64> = row
                        .iter()
                        .enumerate()
                        .map(|(col_idx, &col_val)| (1.0 + selection_coefficients[col_idx] * col_val as f64).ln())
                        .collect();

                    //println!("log_values: {:?}", log_values);
                    
                    //println!("log_values: {:?}", log_values);
                    let neg_inf = log_values.contains(&std::f64::NEG_INFINITY);

                    //let log_mean = log_values.into_iter().map(|x| x).ln_sum_exp() - (row.len() as f64).ln();
                    let mut log_sum = 0.0;
                    if neg_inf == false {
                        log_sum = log_values.iter().sum();
                    }

                    log_sum
                })
                .collect();

            // TODO: work out why when prop_positive = 0, genome size still increases (should favour reduction in genome size)
            let logsumexp_value = selection_weights.iter().ln_sum_exp();

            //println!("logsumexp_value: {:?}", logsumexp_value);

            //println!("raw_selection_weights: {:?}", selection_weights);

            // Exponentiate and normalize
            selection_weights = selection_weights.into_iter()
                .map(|x| (x - logsumexp_value).exp()) // exp(log(w) - logsumexp)
                .collect();

            //println!("pre_norm_selection_weights: {:?}", selection_weights);

            let sum_weights: f64 = selection_weights.iter().sum();
            //println!("sum_weights: {:?}", sum_weights);
            selection_weights = selection_weights.iter().map(|&w| if w != std::f64::NEG_INFINITY {w / sum_weights} else {0.0}).collect();
        }

        // Convert differences to weights (lower difference should have higher weight)
        //println!("raw_weights: {:?}", selection_weights);
        let mut weights : Vec<f64>;
        if no_control_genome_size == false {
            // Calculate the differences from avg_gene_freq
            let differences: Vec<i32> = num_genes
                .iter()
                .map(|&n_genes| (n_genes - avg_gene_num).abs())
                .collect();

            //println!("differences: {:?}", differences);

            weights = differences
                .iter()
                .enumerate()
                .map(|(row_idx, &diff)| genome_size_penalty.powi(diff) * selection_weights[row_idx]) // based on https://pmc.ncbi.nlm.nih.gov/articles/instance/5320679/bin/mgen-01-38-s001.pdf
                .collect();
        } else {
            weights = selection_weights.clone();
        }

        // normalise pairwise dists to minimum value
        // let norm_avg_pairwise_dists = normalize_avg_dists(&avg_pairwise_dists);

        // println!("post_genome_size_weights: {:?}", weights);
        // update weights with average pairwise distance
        for i in 0..weights.len() {
            //let scaled_distance = safe_pow(norm_avg_pairwise_dists[i], 1.0 / competition_strength);
            let scaled_distance = weights[i] * avg_pairwise_dists[i];
            weights[i] = scaled_distance;
        }

        // println!("norm_avg_pairwise_dists: {:?}", norm_avg_pairwise_dists);
        // println!("post_pairwise_weights: {:?}", weights);
        // let mean_avg_pairwise_dists = average(&avg_pairwise_dists);
        // println!("mean_avg_pairwise_dists: {:?}", mean_avg_pairwise_dists);

        // determine whether weights is only 0s
        let max_final_weights = weights.iter().cloned().fold(-1./0. /* -inf */, f64::max);

        // if verbose {
        //     // printing selection weights pre-size selection
        //     // println!("selection_weights: {:?}", selection_weights);
        //     // //println!("selection_coefficients: {:?}", selection_coefficients);
        //     // let max_selection_coefficients = selection_coefficients.iter().cloned().fold(-1./0. /* -inf */, f64::max);
        //     // println!("max_selection_coefficients: {:?}", max_selection_coefficients);
        //     // let min_selection_coefficients = selection_coefficients.iter().copied().fold(f64::INFINITY, f64::min);
        //     // println!("min_selection_coefficients: {:?}", min_selection_coefficients);
        //     // let max_selection_weights = selection_weights.iter().cloned().fold(-1./0. /* -inf */, f64::max);
        //     // println!("max_selection_weights: {:?}", max_selection_weights);
        //     // let min_selection_weights = selection_weights.iter().copied().fold(f64::INFINITY, f64::min);
        //     // println!("min_selection_weights: {:?}", min_selection_weights);
        //     // let mean_selection_weights = average(&selection_weights);
        //     // println!("mean_selection_weights: {:?}", mean_selection_weights);

        //     //println!("differences: {:?}", differences);

        //     // printing selection weights post-size selection
        //     //println!("final_weights: {:?}", weights);
        //     let mean_final_weights = average(&weights);
        //     println!("mean_final_weights: {:?}", mean_final_weights);
        //     println!("max_final_weights: {:?}", max_final_weights);
        //     let min_final_weights = weights.iter().copied().fold(f64::INFINITY, f64::min);
        //     println!("min_final_weights: {:?}", min_final_weights);

        //     //println!("max_diff: {:?}", max_diff);
        //     //println!("weights: {:?}", weights);
        // }

        // account for only zeros
        if max_final_weights == 0.0 {
            weights = vec![1.0; self.pop.nrows()];
        }

        // Create a WeightedIndex distribution based on weights
        let dist = WeightedIndex::new(&weights).unwrap();

        // Sample rows based on the distribution
        let sampled_indices: Vec<usize> = (0..self.pop.nrows()).map(|_| dist.sample(rng)).collect();

        //println!("sampled_indices: {:?}", sampled_indices);

        sampled_indices
    }

    pub fn next_generation(&mut self, sample: &Vec<usize>) {
        let nrows = sample.len();
        let ncols = self.pop.ncols();
    
        // Pre-allocate the new population array
        let mut next_pop = Array2::zeros((nrows, ncols));
    
        // Fill in each row directly
        for (row_idx, &sample_idx) in sample.iter().enumerate() {
            let source_row = self.pop.slice(s![sample_idx, ..]);
            let mut target_row = next_pop.slice_mut(s![row_idx, ..]);
            target_row.assign(&source_row);
        }
    
        self.pop = next_pop;
    }

    pub fn mutate_alleles(
        &mut self,
        mutations_vec: &Vec<f64>,
        weighted_dist: &Vec<WeightedIndex<f32>>,
    ) {
        // index for random number generation
        let _index = AtomicUsize::new(0);
        let _update_rng = AtomicUsize::new(0);

        for site_idx in 0..mutations_vec.len() {
            let mutations = mutations_vec[site_idx];
            
            // avoid rate parameter of 0
            if mutations == 0.0 {
                continue;
            }
            
            let poisson = Poisson::new(mutations).unwrap();

            if self.core == false {
                // generate Poisson sampler
                self.pop
                    .axis_iter_mut(Axis(0))
                    .into_par_iter()
                    .for_each(|mut row| {
                        // thread-specific random number generator
                        let mut thread_rng = rand::thread_rng();
                        //let thread_index = rayon::current_thread_index();
                        //print!("{:?} ", thread_index);

                        // sample from Poisson distribution for number of sites to mutate in this isolate
                        let n_sites = poisson.sample(&mut thread_rng) as usize;

                        // iterate for number of mutations required to reach mutation rate
                        for _ in 0..n_sites {
                            // sample new site to mutate
                            let mutant_site = weighted_dist[site_idx].sample(&mut thread_rng);
                            let value = row[mutant_site];
                            let new_allele: u8 = if value == 0 as u8 { 1 } else { 0 };

                            // set value in place
                            row[mutant_site] = new_allele;
                        }
                    });
            } else {
                self.pop
                    .axis_iter_mut(Axis(0))
                    .into_par_iter()
                    .for_each(|mut row| {
                        // thread-specific random number generator
                        let mut thread_rng = rand::thread_rng();
                        //let thread_index = rayon::current_thread_index();
                        //print!("{:?} ", thread_index);

                        // sample from Poisson distribution for number of sites to mutate in this isolate
                        let n_sites = thread_rng.sample(poisson) as usize;

                        // iterate for number of mutations required to reach mutation rate
                        for _ in 0..n_sites {
                            // sample new site to mutate
                            let mutant_site = weighted_dist[site_idx].sample(&mut thread_rng);

                            // get possible values to mutate to, must be different from current value
                            let value = row[mutant_site];
                            let values = &self.core_vec[1 >> value];

                            // sample new allele
                            let new_allele = values.iter().choose_multiple(&mut thread_rng, 1)[0];

                            // set value in place
                            row[mutant_site] = *new_allele;
                        }
                    });
            }
        }
    }

    pub fn recombine(
        &mut self,
        recombinations_vec: &Vec<f64>,
        rng: &mut StdRng,
        locus_weights: &Vec<Vec<f32>>,
    ) {
        // index for random number generation
        let _index = AtomicUsize::new(0);
        let _update_rng = AtomicUsize::new(0);

        for site_idx in 0..recombinations_vec.len() {
            let n_recombinations = recombinations_vec[site_idx];

            // avoid rate parameter of 0
            if n_recombinations == 0.0 {
                continue;
            }

            let poisson_recomb = Poisson::new(n_recombinations).unwrap();

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
            self.pop
                .axis_iter(Axis(0))
                .into_par_iter()
                .enumerate()
                .for_each(|(row_idx, row)| {
                    //use std::time::Instant;
                    //let now = Instant::now();

                    // thread-specific random number generator
                    let mut thread_rng = rand::thread_rng();

                    // sample from Poisson distribution for number of sites to mutate in this isolate
                    let n_sites = poisson_recomb.sample(&mut thread_rng) as usize;

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
                    let sampled_recipients: Vec<usize> = (0..n_sites)
                        .map(|_| dist.sample(&mut thread_rng))
                        .map(|value| value + (value >= row_idx) as usize)
                        .collect();
                    let mut sampled_loci: Vec<usize> = Vec::with_capacity(n_sites);

                    // for _ in 0..n_sites {
                    //     let value = sample_dist.sample(&mut thread_rng);
                    //     sampled_recipients.push(value + (value >= row_idx) as usize);
                    // }

                    //elapsed = now.elapsed();
                    //println!("finished sampling total: {}, {:.2?}", row_idx, elapsed);

                    //let sampled_recipients: Vec<usize> = vec![1, 7, 20, 705, 256];

                    let mut sampled_values: Vec<u8> = vec![1; n_sites];
                    // get non-zero indices
                    if self.core == false {
                        //if accessory, set elements with no genes to 0
                        let mut non_zero_weights: Vec<f32> = locus_weights[site_idx].clone();
                        let mut total_0: usize = 0;
                        for (idx, &val) in row.indexed_iter() {
                            let mut update: bool = false;
                            
                            // check for any zeroes in either vector
                            if non_zero_weights[idx] == 0.0 {
                                update = true;
                            }

                            if val == 0 {
                                non_zero_weights[idx] = 0.0;
                                update = true;
                            }

                            if update == true 
                            {
                                total_0 += 1;
                            }
                        }

                        // // get all sites to be recombined
                        // sampled_loci = (0..n_sites)
                        // .map(|_| *non_zero_indices.choose(&mut thread_rng).unwrap()) // Sample with replacement
                        // .collect(); 
                        
                        //let sum : f32 = non_zero_weights.clone().iter().sum();
                        // if sum == 0.0 {
                        //     println!("total_0: {:?}", total_0);
                        //     println!("non_zero_weights.len(): {:?}", non_zero_weights.len());
                            
                        //     println!("non_zero_weights.sum(): {:?}", sum);
                        // }


                        // ensure some non-zero values present
                        if total_0 < non_zero_weights.len() {
                            let locus_weighted_dist: WeightedIndex<f32> =
                                WeightedIndex::new(non_zero_weights).unwrap();

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


                    // assign values
                    {
                        let mut mutex = loci.write().unwrap(); // Lock for writing
                        let entry = &mut mutex[row_idx]; // Now you can index safely
                        *entry = sampled_loci;
                    }
                    {
                        let mut mutex = values.write().unwrap(); // Lock for writing
                        let entry = &mut mutex[row_idx]; // Now you can index safely
                        *entry = sampled_values;
                    }
                    {
                        let mut mutex = recipients.write().unwrap(); // Lock for writing
                        let entry = &mut mutex[row_idx]; // Now you can index safely
                        *entry = sampled_recipients;
                    }

                    //elapsed = now.elapsed();
                    //println!("finished entering data: {}, {:.2?}", row_idx, elapsed);
                });

            // go through entries in loci, values and recipients, mutating the rows in each case
            // randomise order in which rows are moved through
            let mut row_indices: Vec<usize> = (0..self.pop.nrows()).collect();
            row_indices.shuffle(rng);

            for pop_idx in row_indices {
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

        }
    }

    pub fn average_distance(&mut self) -> Vec<f64> {
        //let (contiguous_array, matches) = get_variable_loci(self.core, &self.pop);

        let range = 0..self.pop.nrows();
        let distances: Vec<f64> = range
            .into_par_iter()
            .map(|i| {
                let i_distances = get_distance(
                    i,
                    self.pop.nrows(),
                    self.core_genes,
                    0.0,
                    self.core,
                    &self.pop,
                    self.pop.ncols(),
                );

                let (sum, count) = i_distances.iter().fold((0.0, 0), |(s, c), &x| (s + x, c + 1));
                let mut _final_distance = sum / count as f64;

                // ensure no zero distances that may cause no selection of isolates.
                if _final_distance == 0.0 {
                    _final_distance = MIN_POSITIVE;
                }

                _final_distance
            })
            .collect();

        //println!("new distances: {:?}", distances);
        distances
    }

    // TODO update vector in place
    pub fn pairwise_distances(
        &mut self,
        max_distances: usize,
        range1: &Vec<usize>,
        range2: &Vec<usize>,
    ) -> Vec<f64> {
        //let (contiguous_array, matches) = get_variable_loci(self.core, &self.pop);

        //let mut idx = 0;
        let range = 0..max_distances;
        let distances: Vec<_> = range
            .into_par_iter()
            .map(|current_index| {
                let i = range1[current_index];
                let j = range2[current_index];

                let row1 = self.pop.index_axis(Axis(0), i);
                let row2 = self.pop.index_axis(Axis(0), j);
                let row1_slice = row1.as_slice().unwrap();
                let row2_slice = row2.as_slice().unwrap();
                //println!("i: {:?}", i);
                //println!("j: {:?}", j);

                //println!("rowi: {:?}", row1);
                //println!("rowj: {:?}", row2);

                let mut _final_distance: f64 = 0.0;

                if self.core == true {
                    
                    let distance = hamming_bitwise_fast(row1_slice, row2_slice) / 2;
           
                    // let test_distance = row1_slice.to_vec().iter().zip(&row2_slice.to_vec()).filter(|&(a, b)| a != b).count();
                    // println!("distance: {:?} test_distance: {:?}", distance, test_distance);

                    _final_distance = distance as f64 / (self.pop.ncols() as f64);
                } else {
                    let (intersection, union) = jaccard_distance_fast(row1_slice, row2_slice);
                    // let (intersection_test, union_test) = jaccard_distance_test(row1_slice, row2_slice);
                    // println!("intersection_test: {:?} intersection: {:?}", intersection_test, intersection);
                    // println!("union_test: {:?} union: {:?}", union_test, union);
                    _final_distance = 1.0
                        - ((intersection as f64 + self.core_genes as f64)
                            / (union as f64 + self.core_genes as f64));
                }
                //println!("_final_distance: {:?}", _final_distance);
                _final_distance
            })
            .collect();
        distances
    }

    pub fn write(&mut self, outpref: &str) -> io::Result<()>
    {
        // core genome
        if self.core == true {
            // Open output file
            let mut output_file = outpref.to_owned();
            let extension: &str = "_core_genome.csv";
            output_file.push_str(extension);

            let mut file = File::create(output_file)?;

            // Iterate rows and write
            for row in self.pop.outer_iter() {
                let line: Vec<String> = row.iter().map(|&x| int_to_base(x).to_string()).collect();
                writeln!(file, "{}", line.join(","))?;
            }
        } else {
            // Open output file
            let mut output_file = outpref.to_owned();
            let extension: &str = "_pangenome.csv";
            output_file.push_str(extension);

            let mut file = File::create(output_file)?;

            // Iterate rows and write
            for row in self.pop.outer_iter() {
                let mut line: Vec<String> = vec![1_u8.to_string(); self.core_genes];
                line.extend(row.iter().map(|&x| x.to_string()));
                writeln!(file, "{}", line.join(","))?;
            }
        }
        Ok(())
    }
}