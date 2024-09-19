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
use std::sync::atomic::{AtomicUsize, Ordering};

use std::fs::File;
use std::io::{self, Write};
use std::usize;
use clap::{Arg, Command};

fn hamming_distance(x: &[u8], y: &[u8]) -> u64 {
    assert_eq!(x.len(), y.len(), "Vectors must have the same length");
    x.iter().zip(y).filter(|&(a, b)| a != b).count() as u64
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
    fn new(size: usize, allele_count: usize, max_variants: u8, core : bool, avg_gene_freq: f64, rng : &mut StdRng, core_genes : usize) -> Self {
        //let mut pop = Array2::<u8>::zeros((size, allele_count));

        let mut start: Array1<u8> = Array1::zeros(allele_count);

        if core {
            for j in 0..allele_count {
                start[j] = rng.gen_range(0..max_variants);
            }
        } else {
            // let sample_item: [(u8, f64); 2] = [(0, 1.0 - avg_gene_freq), (1, avg_gene_freq)];
            // //Create the WeightedIndex distribution using the weights
            // let weights = sample_item.iter().map(|&(_, weight)| weight).collect::<Vec<_>>();
            // println!("weights:\n{:?}", weights);
            // let weighted_dist = WeightedIndex::new(&weights).unwrap();
            // //let weighted_dist = WeightedIndex::new(sample_item.iter().map(|(_, weight)| weight)).unwrap();
            // for j in 0..allele_count {
            //     start[j] = sample_item[weighted_dist.sample(rng)].0;
            // }

            // set determined number of genes to 1
            let gene_num : usize = (allele_count as f64 * avg_gene_freq).round() as usize;
            
            for j in 0..gene_num {
                start[j] = 1;
            }
        }

        // Create a vector of Array1 filled with the same values
        let pop_vec: Vec<Array1<u8>> = std::iter::repeat(start)
            .take(size)
            .collect();

        // convert vector into 2D array
        let pop = to_array2(pop_vec).unwrap();

        // if !core {
        //     //println!("pop initial:\n{:?}", pop);

        //     let proportions: Vec<f64> = pop.axis_iter(Axis(0))
        //     .map(|row| {
        //         let sum: usize = row.iter().map(|&x| x as usize).sum();
        //         let count = row.len();
        //         sum as f64 / count as f64
        //     })
        //     .collect();
    
        //     println!("proportions initial:\n{:?}", proportions);
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

    fn sample_indices (&mut self, rng : &mut StdRng) -> Vec<usize> {
        // Calculate the proportion of 1s for each row
        let proportions: Vec<f64> = self.pop.axis_iter(Axis(0))
        .map(|row| {
            let sum: usize = row.iter().map(|&x| x as usize).sum();
            let count = row.len();
            sum as f64 / count as f64
        })
        .collect();

        //println!("proportions:\n{:?}", proportions);

        // Calculate the differences from avg_gene_freq
        let differences: Vec<f64> = proportions.iter()
        .map(|&prop| (prop - self.avg_gene_freq).abs())
        .collect();

        //println!("differences:\n{:?}", differences);

        // Convert differences to weights (lower difference should have higher weight)
        let max_diff = differences.iter().cloned().fold(0./0., f64::max);
        let weights: Vec<f64> = differences.iter()
            .map(|&diff| max_diff - diff) // Inverting so smaller differences give larger weights
            .collect();

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

    fn mutate_alleles(&mut self, mutations : i32, rng : &mut StdRng, weighted_dist: &WeightedIndex<f32>) {
        // index for random number generation
        let _index = AtomicUsize::new(0);
        let _update_rng = AtomicUsize::new(0);

        let poisson = Poisson::new(mutations as f64).unwrap();

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
                    let mutant_site = weighted_dist.sample(&mut thread_rng);
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
                    _update_rng.fetch_add(1, Ordering::SeqCst);

                    // sample from Poisson distribution for number of sites to mutate in this isolate
                    let n_sites = thread_rng.sample(poisson) as usize;

                    // iterate for number of mutations required to reach mutation rate
                    for _ in 0..n_sites {
                        // sample new site to mutate
                        let mutant_site = weighted_dist.sample(&mut thread_rng);
                        _update_rng.fetch_add(1, Ordering::SeqCst);

                        // get possible values to mutate to, must be different from current value
                        let value = row[mutant_site];
                        let values = &self.core_vec[value as usize];

                        // sample new allele
                        let new_allele = values.iter().choose_multiple(&mut thread_rng, 1)[0];
                        _update_rng.fetch_add(1, Ordering::SeqCst);

                        // TODO update rng for all times thread_rng is sampled!

                        // set value in place
                        row[mutant_site] = *new_allele;
                    }
                });
        }
        // update rng in place
        let rng_index: usize = _update_rng.load(Ordering::SeqCst);
        //print!("{:?} ", rng_index);
        for _ in 0..rng_index {
            rng.gen::<u64>(); // Discard some numbers to mimic jumping
        }

    }

    fn pairwise_distances(&mut self, max_distances : usize, range1: &Vec<usize>, range2: &Vec<usize>) -> Vec<f64> {
        // determine which columns are all equal, ignore from distance calculations
        let array_f64 = self.pop.mapv(|x| x as f64);
        let column_variance = array_f64.var_axis(Axis(0), 0.0);

        // Determine which column indices have variance greater than 0
        let columns_to_iter: Vec<usize> = column_variance
            .iter()
            .enumerate()
            .filter_map(|(i, &variance)| if variance > 0.0 { Some(i) } else { None })
            .collect();

        // get matches for jaccard distance calculation
        let mut matches: f64 = 0.0;
        if self.core != true {
            matches = self.pop.axis_iter(Axis(1)).filter(|col| col.iter().all(|&x| x == 1)).count() as f64;
        }
        
        let subset_array: Array2<u8> = self.pop.select(Axis(1), &columns_to_iter);
        let mut contiguous_array = Array2::zeros((subset_array.dim().0, subset_array.dim().1));
        contiguous_array.assign(&subset_array);
        //println!("{:?}", self.pop);
        //println!("{:?}", contiguous_array);

        //let mut idx = 0;
        let range = 0..max_distances;
        let distances: Vec<_> = range.into_par_iter().map(|current_index| {
            let i = range1[current_index];
            let j = range2[current_index];

            //println!("i:\n{:?}", i);
            //println!("j:\n{:?}", j);
            
            let row1 = contiguous_array.index_axis(Axis(0), i);
            let row2 = contiguous_array.index_axis(Axis(0), j);

            //println!("rowi:\n{:?}", row1);
            //println!("rowj:\n{:?}", row2);

            let mut _final_distance: f64 = 0.0;

            if self.core == true {
                let distance = hamming_distance(row1.as_slice().unwrap(), &row2.as_slice().unwrap());
                _final_distance = distance as f64 / (column_variance.len() as f64);
            } else {
                let (intersection, union) = jaccard_distance(&row1.as_slice().unwrap(), &row2.as_slice().unwrap());
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
    .version("0.0.1")
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
    .arg(Arg::new("pan_mu")
        .long("pan_mu")
        .help("Maximum average pairwise pangenome distance to achieve by end of simulation.")
        .required(false)
        .default_value("0.2"))
    .arg(Arg::new("proportion_fast")
        .long("proportion_fast")
        .help("Proportion of genes in pangenome in fast compartment.")
        .required(false)
        .default_value("0.5"))
    .arg(Arg::new("speed_fast")
        .long("speed_fast")
        .help("Proportional difference in speed of mutation between slow and fast genes. Must be >=1.0.")
        .required(false)
        .default_value("2.0"))
    .arg(Arg::new("seed")
        .long("seed")
        .help("Seed for random number generation.")
        .required(false)
        .default_value("0"))
    .arg(Arg::new("output")
        .long("output")
        .help("Output file path.")
        .required(false)
        .default_value("distances.tsv"))
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
    let n_gen: i32 = matches.value_of_t("n_gen").unwrap();
    let output = matches.value_of("output").unwrap_or("distances.tsv");
    let max_distances: usize = matches.value_of_t("max_distances").unwrap();
    let core_mu: f64 = matches.value_of_t("core_mu").unwrap();
    let pan_mu: f64 = matches.value_of_t("pan_mu").unwrap();
    let proportion_fast: f32 = matches.value_of_t("proportion_fast").unwrap();
    let speed_fast: f32 = matches.value_of_t("speed_fast").unwrap();
    let mut n_threads: usize = matches.value_of_t("threads").unwrap();
    let verbose = matches.is_present("verbose");
    let seed: u64 = matches.value_of_t("seed").unwrap();

    //let verbose = true;

    // time testing
    //use std::time::Instant;
    //let now = Instant::now();

    // validate all variables
    if core_genes > pan_genes {
        println!("core_genes must be less than or equal to pan_size");
        return Ok(())
    }

    if (proportion_fast < 0.0) || (proportion_fast > 1.0) {
        println!("proportion_fast must be between 0.0 and 1.0");
        return Ok(())
    }

    if speed_fast < 1.0 {
        println!("speed_fast must be above 1.0");
        return Ok(())
    }

    if (pop_size < 1) || (core_size < 1) || (pan_genes < 1) || (n_gen < 1) || (max_distances < 1) {
        println!("pop_size, core_size, pan_genes, n_gen and max_distances must all be above 1");
        return Ok(())
    }

    if (core_mu < 0.0) || (core_mu > 1.0) {
        println!("core_mu must be between 0.0 and 1.0");
        return Ok(())
    }

    if (pan_mu < 0.0) || (pan_mu > 1.0) {
        println!("pan_mu must be between 0.0 and 1.0");
        return Ok(())
    }

    if (avg_gene_freq <= 0.0) || (avg_gene_freq > 1.0) {
        println!("avg_gene_freq must be above 0.0 and below or equal to 1.0");
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
    
    // calculate number of mutations per genome per generation
    let n_core_mutations = (((core_size as f64 * core_mu) / n_gen as f64) / 2.0).ceil() ;
    let n_pan_mutations = (((pan_size as f64 * pan_mu) / n_gen as f64) / 2.0).ceil();

    // set weights for sampling of sites
    let core_weights : Vec<f32> = vec![1.0; core_size];
    let mut pan_weights : Vec<f32> = vec![1.0; pan_size];

    // calculate sites for fast accessory genome
    let num_fast_sites = (pan_size as f32 * proportion_fast).round() as usize;
    for i in 0..num_fast_sites {
        pan_weights[i] = speed_fast;
    }

    let mut rng: StdRng = StdRng::seed_from_u64(seed);

    let mut core_genome = Population::new(pop_size, core_size, 4, true, avg_gene_freq, &mut rng, core_genes); // core genome alignment
    let mut pan_genome = Population::new(pop_size, pan_size, 2, false, avg_gene_freq, &mut rng, core_genes); // pangenome alignment

    // weighted distribution samplers
    let core_weighted_dist = WeightedIndex::new(core_weights).unwrap();
    let pan_weighted_dist = WeightedIndex::new(pan_weights).unwrap();

    for j in 0..n_gen { // Run for n_gen generations
        //let now_gen = Instant::now();
        
        // sample new individuals if not at first generation
        if j > 0 {
            //let sampled_individuals: Vec<usize> = (0..pop_size).map(|_| rng.gen_range(0..pop_size)).collect();
            let sampled_individuals = pan_genome.sample_indices(&mut rng);
            core_genome.next_generation(& sampled_individuals);
            //println!("finished copying core genome {}", j);
            pan_genome.next_generation(& sampled_individuals);
            //println!("finished copying pangenome {}", j);
        }
        
        // if at final generation, just sample, otherwise mutate
        if j < (n_gen - 1) {
            // mutate core genome
            //println!("started {}", j);
            core_genome.mutate_alleles(n_core_mutations as i32, &mut rng, &core_weighted_dist);

            //println!("finished mutating core genome {}", j);
            pan_genome.mutate_alleles(n_pan_mutations as i32, &mut rng, &pan_weighted_dist);
            //println!("finished mutating pangenome {}", j);
        } else {
            let final_avg_gene_freq = pan_genome.calc_gene_freq();
            if verbose {
                println!("final avg_gene_freq: {}", final_avg_gene_freq);
            }
            
            // else calculate hamming and jaccard distances
            // generate random numbers to sample indices
            let range1: Vec<usize> = (0..max_distances).map(|_| rng.gen_range(0..pop_size)).collect();
            let range2: Vec<usize> = (0..max_distances).map(|_| rng.gen_range(0..pop_size)).collect();

            let core_distances = core_genome.pairwise_distances(max_distances, &range1, &range2);
            let acc_distances = pan_genome.pairwise_distances(max_distances, &range1, &range2);

            let mut file = File::create(output)?;

            // Iterate through the vectors and write each pair to the file
            for (core, acc) in core_distances.iter().zip(acc_distances.iter()) {
                writeln!(file, "{}\t{}", core, acc);
            }
        }

        //let elapsed = now_gen.elapsed();
        if verbose {
            println!("Finished gen: {}", j + 1);    
        }
        //println!("Elapsed: {:.2?}", elapsed);
    }
    //let elapsed = now.elapsed();
    
    // if verbose {
    //     println!("Total elapsed: {:.2?}", elapsed);
    // }

    // println!("Total elapsed: {:.2?}", elapsed);

    return Ok(())
}
