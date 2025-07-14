use rand::prelude::*;
use rand::rngs::StdRng;
use rand::distributions::Uniform;
use statrs::distribution::Exp;
use rand::distributions::WeightedIndex;
use rand::{SeedableRng};

use pansim::population::*;

use clap::{Arg, Command};

use std::fs::File;
use std::io::{self, Write};

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
        .help("Average number of accessory pangenome that mutates per generation in gene compartment 1. Must be >= 0.0")
        .required(false)
        .default_value("1.0"))
    .arg(Arg::new("rate_genes2")
        .long("rate_genes2")
        .help("Average number of accessory pangenome that mutates per generation in gene compartment 2. Must be >= 0.0")
        .required(false)
        .default_value("1000.0"))
    .arg(Arg::new("prop_genes2")
        .long("prop_genes2")
        .help("Proportion of pangenome made up of compartment 2 genes. Must be 0.0 <= X <= 1.0")
        .required(false)
        .default_value("0.1"))
    .arg(Arg::new("prop_positive")
        .long("prop_positive")
        .help("Proportion of pangenome made up of positively selected genes. Must be 0.0 <= X <= 1.0. If negative, neutral selection is simulated.")
        .required(false)
        .default_value("-0.1"))
        .allow_hyphen_values(true)
    .arg(Arg::new("pos_lambda")
        .long("pos_lambda")
        .help("Lambda value for exponential distribution of positively selected genes. Must be > 0.0")
        .required(false)
        .default_value("0.000001"))
    .arg(Arg::new("neg_lambda")
        .long("neg_lambda")
        .help("Lambda value for exponential distribution of negatively selected genes. Must be > 0.0.")
        .required(false)
        .default_value("5.0"))
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
    .arg(Arg::new("print_matrices")
        .long("print_matrices")
        .required(false)
        .takes_value(false))
        .help("Prints core and accessory matrices.")
    .arg(Arg::new("print_selection")
        .long("print_selection")
        .required(false)
        .takes_value(false))
        .help("Prints selection coefficients.")
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
    .arg(Arg::new("no_control_genome_size")
        .long("no_control_genome_size")
        .help("Removes penalisation of genome sizes deviating from average.")
        .required(false)
        .takes_value(false))
    .arg(Arg::new("genome_size_penalty")
        .long("genome_size_penalty")
        .help("Multiplier for each gene difference between avg_gene_freq and observed value. Default = 0.99")
        .required(false)
        .default_value("0.99"))
    .arg(Arg::new("competition_strength")
        .long("competition_strength")
        .help("Strength of competition felt by strain to all others. Default = 1.0")
        .required(false)
        .default_value("1.0"))
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
    let prop_genes2: f64 = matches.value_of_t("prop_genes2").unwrap();
    let prop_positive: f64 = matches.value_of_t("prop_positive").unwrap();
    let pos_lambda: f64 = matches.value_of_t("pos_lambda").unwrap();
    let neg_lambda: f64 = matches.value_of_t("neg_lambda").unwrap();
    let mut n_threads: usize = matches.value_of_t("threads").unwrap();
    let verbose = matches.is_present("verbose");
    let competition = matches.is_present("competition");
    let seed: u64 = matches.value_of_t("seed").unwrap();
    let print_dist: bool = matches.is_present("print_dist");
    let print_matrices: bool = matches.is_present("print_matrices");
    let print_selection: bool = matches.is_present("print_selection");
    let no_control_genome_size: bool = matches.is_present("no_control_genome_size");
    let genome_size_penalty: f64 = matches.value_of_t("genome_size_penalty").unwrap();
    let competition_strength: f64 = matches.value_of_t("competition_strength").unwrap();

    //let verbose = true;

    // time testing
    //use std::time::Instant;
    //let now = Instant::now();

    // validate all variables
    if core_genes > pan_genes {
        println!("core_genes must be less than or equal to pan_size");
        return Ok(());
    }

    if (HR_rate < 0.0 || HGT_rate < 0.0) {
        println!("HR_rate and HGT_rate must be above 0.0");
        println!("HR_rate: {}", HR_rate);
        println!("HGT_rate: {}", HGT_rate);
        return Ok(());
    }

    if (pos_lambda <= 0.0 || neg_lambda <= 0.0) {
        println!("pos_lambda and neg_lambda must be above 0.0");
        println!("pos_lambda: {}", pos_lambda);
        println!("neg_lambda: {}", neg_lambda);
        return Ok(());
    }

    if rate_genes1 < 0.0 || rate_genes2 < 0.0 {
        println!("rate_genes1 and rate_genes2 must be >= 0");
        println!("rate_genes1: {}", rate_genes1);
        println!("rate_genes2: {}", rate_genes2);
        return Ok(());
    }

    if prop_genes2 < 0.0 || prop_genes2 > 1.0 {
        println!("prop_genes2 must be 0.0 <= prop_genes2 <= 1.0");
        println!("prop_genes2: {}", prop_genes2);
        return Ok(());
    }

    if (pop_size < 1) || (core_size < 1) || (pan_genes < 1) || (n_gen < 1) || (max_distances < 1) {
        println!("pop_size, core_size, pan_genes, n_gen and max_distances must all be above 1");
        println!("pop_size: {}", pop_size);
        println!("core_size: {}", core_size);
        println!("pan_genes: {}", pan_genes);
        println!("n_gen: {}", n_gen);
        println!("max_distances: {}", max_distances);
        return Ok(());
    }

    if (core_mu < 0.0) || (core_mu > 1.0) {
        println!("core_mu must be between 0.0 and 1.0");
        println!("core_mu: {}", core_mu);
        return Ok(());
    }

    if (avg_gene_freq <= 0.0) || (avg_gene_freq > 1.0) {
        println!("avg_gene_freq must be above 0.0 and below or equal to 1.0");
        println!("avg_gene_freq: {}", avg_gene_freq);
        return Ok(());
    }

    if n_threads < 1 {
        n_threads = 1;
    }

    // enable multithreading
    rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build_global()
        .unwrap();

    let pan_size = pan_genes - core_genes;

    // adjust avg_gene_freq by size of core genome
    // determine proportion of accessory genome that should be present on average
    let core_prop: f64 = core_genes as f64 / pan_genes as f64;
    let acc_prop: f64 = 1.0 - core_prop;
    avg_gene_freq = (avg_gene_freq - core_prop) / acc_prop;
    if avg_gene_freq < 0.0 {
        avg_gene_freq = 0.0;
    }
    if verbose {
        println!("avg_gene_freq adjusted to {}", avg_gene_freq);
    }
    let avg_gene_num: i32 = (avg_gene_freq * pan_size as f64).round() as i32;

    // calculate number of mutations per genome per generation, should this be whole pangenome or just accessory genes?
    let n_core_mutations =
        vec![(((core_size as f64 * core_mu) / n_gen as f64) / 2.0).ceil() as f64];

    // calculate average recombinations per genome
    let n_recombinations_core: Vec<f64> = vec![(n_core_mutations[0] as f64 * HR_rate).round()];
    let n_recombinations_pan_total = (n_core_mutations[0] as f64 * HGT_rate).round();
    let mut n_recombinations_pan: Vec<f64> = vec![];

    // set weights for sampling of sites
    let core_weights: Vec<Vec<f32>> = vec![vec![1.0; core_size]; 1];
    let mut pan_weights: Vec<Vec<f32>> = vec![];
    let mut n_pan_mutations: Vec<f64> = vec![];
    let mut selection_weights: Vec<f64> = vec![0.0; pan_size];

    let mut rng: StdRng = StdRng::seed_from_u64(seed);
    
    // calculate selection weights for each gene
    if prop_positive >= 0.0 {
        let uniform: Uniform<f64> = Uniform::new(0.0, 1.0);
        let exponential_pos: Exp = Exp::new(pos_lambda).unwrap();
        let exponential_neg: Exp = Exp::new(neg_lambda).unwrap();
        
        for i in 0..pan_size {
            let weight: f64 = uniform.sample(&mut rng) as f64;

            let mut selection_coeffient: f64 = 0.0;

            // positively selected gene
            if weight <= prop_positive {
                selection_coeffient = exponential_pos.sample(&mut rng);
                //selection_coeffient = 100.0;
            } else {
                selection_coeffient = exponential_neg.sample(&mut rng);

                while selection_coeffient > 1.0 {
                    selection_coeffient = exponential_neg.sample(&mut rng);
                }
                // if selection_coeffient > 1.0 {
                //     selection_coeffient = 1.0;
                // }
                selection_coeffient = -1.0 * selection_coeffient;
            }
            selection_weights[i] = selection_coeffient;
        }
    }

    if print_selection {

        let mut output_file = outpref.to_owned();
        let extension: &str = "_selection.tsv";
        output_file.push_str(extension);

        let mut file = File::create(output_file)?;
        let line: Vec<String> = selection_weights.iter().map(|&x| x.to_string()).collect();
                writeln!(file, "{}", line.join("\n"))?;

    }

    // calculate sites for fast accessory genome
    let num_gene1_sites = (pan_size as f64 * (1.0 - prop_genes2)).round() as usize;
    let prop_gene1_sites: f64 = num_gene1_sites as f64 / pan_size as f64;
    let prop_gene2_sites: f64 = 1.0 - prop_gene1_sites;

    // create weights for either rate compartments
    // for gene rate 1 if any exist, otherwise don't add
    if num_gene1_sites > 0 {
        let mut pan_weights_1: Vec<f32> = vec![0.0; pan_size];
        for i in 0..num_gene1_sites {
            pan_weights_1[i] = 1.0;
        }
        pan_weights.push(pan_weights_1);
        n_pan_mutations.push(rate_genes1);
        let n_recombinations_pan_gene1 =
            n_recombinations_pan_total * prop_gene1_sites; //(rate_genes1 / (rate_genes1 + rate_genes2));
        n_recombinations_pan.push(n_recombinations_pan_gene1);
    } 
    
    // gene 2 rates of mutation and HGT, if any genes exist, otherwise don't add
    if num_gene1_sites < pan_size {
        let mut pan_weights_2: Vec<f32> = vec![0.0; pan_size];
        for i in num_gene1_sites..pan_size {
            pan_weights_2[i] = 1.0;
        }
        n_pan_mutations.push(rate_genes2);
        pan_weights.push(pan_weights_2);

        let n_recombinations_pan_gene2 =
            n_recombinations_pan_total * prop_gene2_sites; //(rate_genes2 / (rate_genes1 + rate_genes2));
        n_recombinations_pan.push(n_recombinations_pan_gene2);
    }

    // generate sampling distribution for genes in accessory genome
    let acc_sampling_vec = sample_beta(pan_size, &mut rng);

    let mut core_genome = Population::new(
        pop_size,
        core_size,
        4,
        true,
        avg_gene_freq,
        &mut rng,
        core_genes,
        &acc_sampling_vec,
    ); // core genome alignment
    let mut pan_genome = Population::new(
        pop_size,
        pan_size,
        2,
        false,
        avg_gene_freq,
        &mut rng,
        core_genes,
        &acc_sampling_vec,
    ); // pangenome alignment

    // weighted distribution samplers
    let core_weighted_dist: Vec<WeightedIndex<f32>> = core_weights
        .clone()
        .into_iter()
        .map(|dist| WeightedIndex::new(dist.clone()).unwrap())
        .collect();
    let pan_weighted_dist: Vec<WeightedIndex<f32>> = pan_weights
        .clone()
        .into_iter()
        .map(|dist| WeightedIndex::new(dist.clone()).unwrap())
        .collect();

    // hold pairwise core and accessory distances per generation
    let mut avg_acc_dist = vec![0.0; n_gen as usize];
    let mut avg_core_dist = vec![0.0; n_gen as usize];
    let mut std_acc_dist = vec![0.0; n_gen as usize];
    let mut std_core_dist = vec![0.0; n_gen as usize];

    // generate random numbers to sample indices
    // TODO make it so that equivalent distances aren't sample, sample with replacement from one?
    let range1: Vec<usize> = (0..max_distances)
        .map(|_| rng.gen_range(0..pop_size))
        .collect();
    let mut range2: Vec<usize> = vec![0; max_distances];

    let mut i2 = 0;
    // sample same range, ensure self-comparisons not included
    for i1 in range1.clone() {
        let mut entry = rng.gen_range(0..pop_size - 1);
        if entry >= i1 {
            entry += 1
        };
        range2[i2] = entry;
        i2 += 1;
    }

    for j in 0..n_gen {
        // Run for n_gen generations
        //let now_gen = Instant::now();

        // sample new individuals
        //let sampled_individuals: Vec<usize> = (0..pop_size).map(|_| rng.gen_range(0..pop_size)).collect();
        let mut avg_pairwise_dists = vec![1.0; pop_size];

        // include competition
        if competition == true {
            avg_pairwise_dists = pan_genome.average_distance();
        }

        let sampled_individuals =
            pan_genome.sample_indices(&mut rng, avg_gene_num, avg_pairwise_dists, &selection_weights, verbose, no_control_genome_size, genome_size_penalty, competition_strength);
        
        core_genome.next_generation(&sampled_individuals);
        //println!("finished copying core genome {}", j);
        pan_genome.next_generation(&sampled_individuals);
        //println!("finished copying pangenome {}", j);

        // mutate core genome
        //println!("started {}", j);
        core_genome.mutate_alleles(&n_core_mutations, &core_weighted_dist);

        //println!("finished mutating core genome {}", j);
        pan_genome.mutate_alleles(&n_pan_mutations, &pan_weighted_dist);
        //println!("finished mutating pangenome {}", j);

        // recombine populations
        if HR_rate > 0.0 {
            core_genome.recombine(&n_recombinations_core, &mut rng, &core_weights);
        }
        if HGT_rate > 0.0 {
            pan_genome.recombine(&n_recombinations_pan, &mut rng, &pan_weights);
        }
        
        // if at final generation, sample
        if j == n_gen -1 {

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
            let avg_gene_freq = pan_genome.calc_gene_freq();
            println!("avg_gene_freq: {}", avg_gene_freq);
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
        for (avg_core, avg_acc, std_core, std_acc) in avg_core_dist
            .iter()
            .zip(avg_acc_dist.iter())
            .zip(std_core_dist.iter())
            .zip(std_acc_dist.iter())
            .map(|(((w, x), y), z)| (w, x, y, z))
        {
            writeln!(file, "{}\t{}\t{}\t{}", avg_core, std_core, avg_acc, std_acc);
        }
    }

    if print_matrices {
        core_genome.write(outpref);
        pan_genome.write(outpref);
    }

    //let elapsed = now.elapsed();

    // if verbose {
    //     println!("Total elapsed: {:.2?}", elapsed);
    // }

    // println!("Total elapsed: {:.2?}", elapsed);

    return Ok(());
}
