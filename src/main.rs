use std::process;

use clap::Parser;
use extract_sites::SnpPos;
use fasta_io::FastaFile;
use validate_inputs::Inputs;

mod extract_sites;
mod fasta_io;
mod validate_inputs;

#[derive(Parser)]
#[clap(
    version = "0.1",
    author = "Sudaraka Mallawaarachchi <smallawaarachchi@gmail.com>",
    about = "This is a draft Rust implementation of a SNP extractor. This version loads the whole multi-fasta alignment into memory, which is not efficient!"
)]
pub struct Args {
    #[arg(short = 'i', long = "input")]
    input_filepath: String,
    #[arg(short = 'o', long = "output")]
    output_folder_name: String,
    #[arg(short = 'm', long = "maf", default_value = "0.01")]
    minor_allele_frequency: f32,
}

fn main() {
    let args = cli_args();

    // Parse the inputs
    let input = match Inputs::parse_inputs(
        args.input_filepath,
        args.output_folder_name,
        args.minor_allele_frequency,
    ) {
        Ok(input) => input,
        Err(err) => {
            eprintln!("Error parsing inputs: {}", err);
            process::exit(1);
        }
    };

    println!("{:?}", input);

    // Read the fasta file
    let fasta_file = match FastaFile::read_fasta(input.read_path().to_path_buf()) {
        Ok(fasta_file) => fasta_file,
        Err(err) => {
            eprintln!("Error reading fasta file: {}", err);
            process::exit(1);
        }
    };

    // Extract variants that pass the maf
    let snp_pos = match SnpPos::get_snp_positions(&fasta_file, input.maf()) {
        Ok(snp_pos) => snp_pos,
        Err(err) => {
            eprintln!("Error extracting SNPs: {}", err);
            process::exit(1);
        }
    };
    
    // Write a snp-only fasta file with positions
    let _ = FastaFile::write_fasta(&fasta_file, &snp_pos, &input);


    // println!("<<{:?}>>", fasta_file);

    // let snp_pos = check_invariance(fasta_file.fasta_sites(), fasta_file.seqlen() - 1, fasta_file.num_seqs(), input.maf);
    // let snp_pos = check_invariance(&fasta_file, input.maf);
    // println!("{:?}", snp_pos);

    // write_fasta(args.output_folder_name, &fasta_vec, &snp_pos, num_seqs, &headers);
}

pub fn cli_args() -> Args {
    Args::parse()
}
